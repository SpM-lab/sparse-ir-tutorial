#!/usr/bin/env python3

import sys
import shutil
import re

header = {}
header["python"] = \
"""\
---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
"""

header["julia"] = \
"""\
---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Julia 1.7
  language: julia
  name: julia-1.7
---
"""

def update_meta_data(filename):
    f = open(filename, "r")
    f_out = open(filename+"_update_work", "w")

    lang = ""

    line = f.readline()
    if line[:-1] != '---':
        return
    line = f.readline()
    p = re.compile(r'\s+language:\s+(\S+)')
    while line[:-1] != '---':
        line = f.readline()
        m = p.match(line)
        if m is not None:
            lang = m.groups()[0]
    if lang not in ["python", "julia"]:
        raise RuntimeError("Unknown language: " + lang)
    print(header[lang], end='', file=f_out)

    line = f.readline()
    while line:
        print(line, end='', file=f_out)
        line = f.readline()
    f.close()
    f_out.close()
    shutil.move(filename+"_update_work", filename)


if len(sys.argv) < 2:
    print("Invalid number of arguments!")

for filename in sys.argv[1:]:
    update_meta_data(filename)