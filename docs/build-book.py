#!/usr/bin/env python3

import subprocess
import sys
from glob import glob
from html.parser import HTMLParser
from pathlib import Path
from shutil import copyfile

from jinja2 import Environment, FileSystemLoader, select_autoescape

CONTENT_PATH = Path("docs")
THEME_PATH = Path("docs") / "theme"
OUT_PATH = Path("_site")

THEME_FILES = [
    "ayu-highlight.css",
    "book.js",
    "clipboard.min.js",
    "favicon.png",
    "highlight.css",
    "highlight.js",
    "tomorrow-night.css",
    "css/chrome.css",
    "css/general.css",
    "css/print.css",
    "css/variables.css",
    "FontAwesome/css/font-awesome.min.css",
    "FontAwesome/fonts/fontawesome-webfont.eot",
    "FontAwesome/fonts/fontawesome-webfont.svg",
    "FontAwesome/fonts/fontawesome-webfont.ttf",
    "FontAwesome/fonts/fontawesome-webfont.woff",
    "FontAwesome/fonts/fontawesome-webfont.woff2",
]


class StaticFileFinder:
    def __init__(self):
        self.static_files = []

    def handle_starttag(self, tag, attrs):
        if tag == "img":
            attrs = dict(attrs)
            if "src" in attrs:
                self.static_files.append(attrs["src"])


class HeadingParser(HTMLParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.current_tag = None
        self.toc = []
        self.static_file_finder = StaticFileFinder()

    def add_content(self, content):
        if self.toc:
            if self.current_tag is not None:
                key = "name"
            else:
                key = "content"
            self.toc[-1][key] += content

    def handle_starttag(self, tag, attrs):
        self.static_file_finder.handle_starttag(tag, attrs)
        if tag in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_tag is None
            self.current_tag = tag
            ident = dict(attrs).get("id")
            self.toc.append(
                {
                    "level": int(tag[1]),
                    "name": "",
                    "link": "#" + ident,
                    "content": self.get_starttag_text(),
                }
            )
        else:
            self.add_content(self.get_starttag_text())

    def handle_data(self, data):
        self.add_content(data)

    def handle_endtag(self, tag):
        if tag in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_tag == tag
            self.current_tag = None
        self.add_content("</" + tag + ">")


class Book(object):
    def __init__(self, chapters):
        super().__init__()
        self.chapters = chapters
        self.toc = []
        self.static_files = set()
        self.build_toc()

    def build_toc(self):
        for chap in self.chapters:
            content = get_chapter_content(chap)
            heading_parser = HeadingParser(convert_charrefs=False)
            heading_parser.feed(content)
            toc = heading_parser.toc
            toc[0]["link"] = get_chapter_name(chap)
            self.toc.extend(heading_parser.toc)
            self.static_files |= set(heading_parser.static_file_finder.static_files)


def get_chapters():
    chapters = sorted(glob("docs/[0-9]*-*.html"))
    chapters.insert(0, "docs/index.html")
    return Book(chapters)


def get_chapter_name(chapter):
    return Path(chapter).name


def get_chapter_content(chapter):
    with open(chapter, "r") as file:
        return file.read()


def render(book):
    env = Environment(
        loader=FileSystemLoader(str(THEME_PATH)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("index.html.jinja")

    output = []
    for i, chapter in enumerate(book.chapters):
        content = get_chapter_content(chapter)
        output.append(
            template.render(
                content=content,
                path_to_root="",
                mathjax_support=True,
                default_theme="Light",
                language="en-us",
                title="Mutual Information between Trajectories",
                book_title="Mutual Information between Trajectories",
                path=get_chapter_name(chapter),
                next={"link": get_chapter_name(book.chapters[i + 1])}
                if len(book.chapters) > i + 1
                else None,
                previous={"link": get_chapter_name(book.chapters[i - 1])}
                if i > 0
                else None,
                chapters=book.toc,
            )
        )
    return output


def copy_static_files(book):
    filelist = THEME_FILES
    for file in filelist:
        (OUT_PATH / file).parent.mkdir(parents=True, exist_ok=True)
        dest = copyfile(THEME_PATH / file, OUT_PATH / file)
        print(f"copied static file to {dest}", file=sys.stderr)

    filelist = book.static_files
    for file in filelist:
        (OUT_PATH / file).parent.mkdir(parents=True, exist_ok=True)
        dest = copyfile(CONTENT_PATH / file, OUT_PATH / file)
        print(f"copied static file to {dest}", file=sys.stderr)


def invoke_pandoc():
    chapters = sorted(glob("docs/[0-9]*-*.md"))
    chapters.insert(0, "docs/index.md")
    for chapter in chapters:
        outfile = Path(chapter).with_suffix(".html")
        args = [
            "pandoc",
            "--filter",
            "pandoc-crossref",
            "--filter",
            "pandoc-citeproc",
            "--csl",
            "docs/elsevier-with-titles.csl",
            "--mathjax",
            "--bibliography",
            "docs/library.bib",
            "-o",
            str(outfile),
            chapter,
        ]
        print(" ".join(args), file=sys.stderr)
        subprocess.run(args, capture_output=True, check=True)


def main():
    invoke_pandoc()

    book = get_chapters()
    path = OUT_PATH
    path.mkdir(exist_ok=True)

    copy_static_files(book)

    output = render(book)
    for chapter, result in zip(book.chapters, output):
        chaptername = get_chapter_name(chapter)
        with (path / chaptername).open("w") as file:
            file.write(result)
            print(f"wrote {file.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
