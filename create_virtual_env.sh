#!/bin/bash -i

VENV_NAME=16SrDNA

if [ -d "~/venv/${VENV_NAME}" ]; then
  rm -rf "~/venv/${VENV_NAME}"
fi

# there is a problem with VirtualBox shared filesystems
# that prevents a Python virtual environment from being created
# in a host directory

python3 -m venv ~/venv/${VENV_NAME}
source ~/venv/${VENV_NAME}/bin/activate
pip install --upgrade pip

# install numpy explicitly or biom-format will not install
pip install numpy

pip install -e .[test]
