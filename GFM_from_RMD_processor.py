#!/usr/bin/env python3
from pathlib import Path
from datetime import date
import sys

preamble = '''---
layout: post
title: {title}
tags: [{tags}]
author: {author}
excerpt_separator: <!--more-->
---

'''

def get_contents(filename):
    try:
        with open(filename, 'r') as f:
            return f.read()
    except IOError as e:
        print("File not found", filename)
        raise RuntimeException("Bad filename") from e


def get_image_folder(filename):
    img_folder = Path(filename.stem + '_files')
    return img_folder


def get_values():
    vals = dict()
    vals['title'] = input("What is the title of the post? ")
    vals['tags'] = input("What tags should be added to the post (comma separated list) ")
    vals['author'] = input("What is the Author name? ")
    return vals


def edit_image_links(current, folder_name):
    lines = current.splitlines()
    output = []
    old_path, new_path = str(folder_name), str('/media' / folder_name)
    for line in lines:
        if old_path in line:
            output.append(line.replace(old_path, new_path))
        else:
            output.append(line)
    return '\n'.join(output)            


def update_contents(filename, new_contents):
    with open(filename, 'w') as updated_file:
        updated_file.write(new_contents)


def process_file(filename):
    filename = Path(filename)
    contents = get_contents(filename)
    img_folder = get_image_folder(filename)
    new_contents = preamble.format(**get_values()) + contents
    if img_folder.exists():
        new_contents = edit_image_links(new_contents, img_folder)
        img_folder.rename(Path('media') / img_folder)
    update_contents(filename, new_contents)
    filename.rename(Path('_posts') / filename.with_stem(date.today().strftime("%Y-%m-%d-") + filename.stem))


if __name__ == '__main__':
    print("Attempting to process the filename passed to the script")
    process_file(sys.argv[1])
    print("Add the marker <!--more--> to the line where you want to break the summary view in the new _posts file")
