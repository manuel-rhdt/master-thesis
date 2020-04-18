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
            ident = dict(attrs).get("id", "")
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


from io import StringIO


class ChapterSeparator(HTMLParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.chapters = []
        self.depth = 0

    @property
    def current_chapter(self):
        self.chapters.get(-1)

    def add_subsection_to(self, section, name):
        new_section = {
            "name": name,
            "parent": section,
            "path": section["path"],
            "depth": self.depth + 1,
        }
        if "subsections" in section:
            section["subsections"].append(new_section)
        else:
            section["subsections"] = [new_section]
        return section["subsections"][-1]

    def start_section(self, name):
        if self.depth == 0:
            self.chapters.append(
                {"name": name, "path": name + ".html", "content": StringIO()}
            )
            self.current_section = self.current_chapter
        else:
            self.current_section = self.add_subsection_to(self.current_section, name)
        self.depth += 1

    def end_section(self):
        self.current_section = self.current_section["parent"]

        if self.depth > 0:
            self.depth -= 1

        if self.depth == 0:
            self.current_chapter["content"] = self.current_chapter["content"].getvalue()

    def add_to_chapter(self, content):
        if self.current_chapter is not None:
            try:
                self.current_chapter["content"].write(content)
            except AttributeError:
                pass

    def add_to_heading(self, content):
        if self.current_heading is not None:
            self.current_heading.write(content)

    def handle_html_content(self, content):
        self.add_to_chapter(content)
        self.add_to_heading(content)

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "section":
            self.start_section(dict(attrs)["id"])
        elif tag.lower() in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_heading is None, "Nested Headings are not allowed"
            self.current_heading = StringIO()

        self.handle_html_content(self.get_starttag_text())

    def handle_data(self, data):
        self.handle_html_content(data)

    def handle_endtag(self, tag):
        self.handle_html_content(f"</{tag}>")
        if tag.lower() == "section":
            self.end_section()
        elif tag.lower() in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_heading is not None, "Spurious " + tag + " end tag."
            self.current_section["title"] = self.current_heading.getvalue()
            self.current_heading = None


class Book(object):
    def __init__(self, content):
        super().__init__()
        self.static_files = set()
        self.content = content
        self.build_toc()
        self.build_chapters()

    def build_toc(self):
        heading_parser = HeadingParser(convert_charrefs=False)
        heading_parser.feed(self.content)
        self.toc = heading_parser.toc
        self.static_files |= set(heading_parser.static_file_finder.static_files)

    def build_chapters(self):
        chapter_separator = ChapterSeparator(convert_charrefs=False)
        print(self.content)
        chapter_separator.feed(self.content)
        self.chapters = chapter_separator.chapters
        self.chapter_names = chapter_separator.chapter_names


def render(book):
    env = Environment(
        loader=FileSystemLoader(str(THEME_PATH)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("index.html.jinja")

    output = []
    for i, chapter in enumerate(book.chapters):
        chapter_name = book.chapter_names[i]
        output.append(
            template.render(
                content=chapter,
                path_to_root="",
                mathjax_support=True,
                default_theme="Light",
                language="en-us",
                title="Mutual Information between Trajectories",
                book_title="Mutual Information between Trajectories",
                path=chapter_name,
                next={"link": book.chapter_names[i + 1]}
                if len(book.chapters) > i + 1
                else None,
                previous={"link": book.chapter_names[i - 1]} if i > 0 else None,
                chapters=book.toc,
            )
        )

    print_html = template.render(
        content=book.content,
        path_to_root="",
        mathjax_support=True,
        default_theme="Light",
        language="en-us",
        title="Mutual Information between Trajectories",
        book_title="Mutual Information between Trajectories",
        path=chapter_name,
        next={"link": book.chapter_names[i + 1]}
        if len(book.chapters) > i + 1
        else None,
        previous={"link": book.chapter_names[i - 1]} if i > 0 else None,
        chapters=book.toc,
        is_print=True,
    )

    book.print_html = print_html

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


def invoke_pandoc() -> str:
    chapters = sorted(glob("docs/[0-9]*-*.md"))
    chapters.insert(0, "docs/index.md")

    chapter_html = []
    for chapter_number, chapter in enumerate(chapters):
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
            "--number-sections",
            f"--number-offset={chapter_number}",
            "--from",
            "markdown+smart+auto_identifiers",
            "--section-divs",
            chapter,
        ]
        print(" ".join(args), file=sys.stderr)
        html = subprocess.run(args, stdout=subprocess.PIPE, check=True, encoding="utf8")
        chapter_html.append(html.stdout)
    return "\n".join(chapter_html)


def main():
    html = invoke_pandoc()

    book = Book(html)
    path = OUT_PATH
    path.mkdir(exist_ok=True)

    copy_static_files(book)

    output = render(book)
    for chapter_name, result in zip(book.chapter_names, output):
        with (path / chapter_name).open("w") as file:
            file.write(result)
            print(f"wrote {file.name}", file=sys.stderr)

    with (path / "print.html").open("w") as file:
        file.write(book.print_html)
        print(f"wrote {file.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
