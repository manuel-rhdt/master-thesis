#!/usr/bin/env python3

import subprocess
import sys
from glob import glob
from html.parser import HTMLParser
from io import StringIO
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


class ChapterSeparator(HTMLParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.chapters = []
        self.static_file_finder = StaticFileFinder()
        self.depth = 0
        self.current_chapter = None
        self.current_section = None
        self.current_heading = None

    @property
    def static_files(self):
        return self.static_file_finder.static_files

    def add_subsection_to(self, section, name):
        new_section = {
            "name": name,
            "parent": section,
            "path": self.current_chapter["path"] + f"#{name}",
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
            self.current_chapter = self.chapters[-1]
            self.current_section = self.current_chapter
        else:
            self.current_section = self.add_subsection_to(self.current_section, name)
        self.depth += 1

    def end_section(self):
        self.current_section = self.current_section.get("parent")

        if self.depth > 0:
            self.depth -= 1

        if self.depth == 0:
            self.current_chapter["content"] = self.current_chapter["content"].getvalue()
            self.current_chapter = None

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
        self.static_file_finder.handle_starttag(tag, attrs)
        if tag.lower() == "section":
            self.start_section(dict(attrs)["id"])
        self.handle_html_content(self.get_starttag_text())
        if tag.lower() in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_heading is None, "Nested Headings are not allowed"
            self.current_heading = StringIO()

    def handle_data(self, data):
        self.handle_html_content(data)

    def handle_endtag(self, tag):
        if tag.lower() in ["h1", "h2", "h3", "h4", "h5", "h6"]:
            assert self.current_heading is not None, "Spurious " + tag + " end tag."
            self.current_section["title"] = self.current_heading.getvalue()
            self.current_heading = None
        self.handle_html_content(f"</{tag}>")
        if tag.lower() == "section":
            self.end_section()


class Book(object):
    def __init__(self, content):
        super().__init__()
        self.static_files = set()
        self.content = content
        self.build_chapters()

    def build_chapters(self):
        chapter_separator = ChapterSeparator(convert_charrefs=False)
        chapter_separator.feed(self.content)
        self.chapters = chapter_separator.chapters
        self.static_files = set(chapter_separator.static_files)
        self.content = "\n".join(ch["content"] for ch in self.chapters)


def render(book):
    env = Environment(
        loader=FileSystemLoader(str(THEME_PATH)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("index.html.jinja")

    output = {}
    for i, chapter in enumerate(book.chapters):
        chapter_path = chapter["path"]
        output[chapter_path] = template.render(
            content=chapter["content"],
            path_to_root="",
            mathjax_support=True,
            default_theme="Light",
            language="en",
            title="Mutual Information between Trajectories",
            book_title="Mutual Information between Trajectories",
            path=chapter_path,
            next={"link": book.chapters[i + 1]["path"]}
            if len(book.chapters) > i + 1
            else None,
            previous={"link": book.chapters[i - 1]["path"]} if i > 0 else None,
            chapters=book.chapters,
            git_repository_url="https://github.com/manuel-rhdt/master-thesis",
            git_repository_icon="fa-github",
        )

    print_html = template.render(
        content=book.content,
        path_to_root="",
        mathjax_support=True,
        default_theme="Light",
        language="en-us",
        title="Mutual Information between Trajectories",
        book_title="Mutual Information between Trajectories",
        path="print.html",
        chapters=book.chapters,
        is_print=True,
        git_repository_url="https://github.com/manuel-rhdt/master-thesis",
        git_repository_icon="fa-github",
    )

    output["print.html"] = print_html

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
    for chapter_name, result in output.items():
        with (path / chapter_name).open("w") as file:
            file.write(result)
            print(f"wrote {file.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
