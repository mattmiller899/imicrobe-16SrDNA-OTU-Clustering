#!/bin/bash

if [ -d "~/venv/juITS16S" ]; then
  rm -rf "~/venv/juITS16S"
fi

# there is a problem with VirtualBox shared filesystems
# that prevents a Python virtual environment from being created

/usr/bin/python3 -m venv ~/venv/juITS16S
source ~/venv/juITS16S/bin/activate
pip install --upgrade pip
pip install -e .[test]
