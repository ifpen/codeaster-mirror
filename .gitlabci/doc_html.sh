#!/bin/bash -e

# waiting for sphinx + pydata-sphinx-theme from container...
printf "\nSetup sphinx with pydata theme...\n"
source /opt/public/debian-12_mpi.sh
python3 -m venv .venvdoc
source .venvdoc/bin/activate
python3 -m pip install -U 'pip>=25.1'
python3 -m pip install .
python3 -m pip install --group doc
sed -i -e "s%/.*/sphinx-build%$(pwd)/.venvdoc/bin/sphinx-build%g" build/mpi/c4che/*_cache.py

printf "\nGenerate html documentation...\n"
# do not use the makefile that wipes the venv environment
./waf doc

getstatus="git status -uno --porcelain"
if [ $(${getstatus} | wc -l) != 0 ]; then
    printf "\nChanges must be committed:\n"
    ${getstatus}
    exit 1
fi
