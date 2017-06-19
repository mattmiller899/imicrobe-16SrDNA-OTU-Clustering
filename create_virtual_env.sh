#!/bin/bash

if [ -d "~/venv/16SrDNA" ]; then
  rm -rf "~/venv/16SrDNA"
fi

# there is a problem with VirtualBox shared filesystems
# that prevents a Python virtual environment from being created
# in a host directory

/usr/bin/python3 -m venv ~/venv/16SrDNA
source ~/venv/16SrDNA/bin/activate
pip install --upgrade pip
pip install -e .[test]
