#!/bin/bash -e

echo "+ compiling..."
export DEVTOOLS_COMPUTER_ID=none

# only add mpi4py
. env.d/version.sh
export PREREQ_PATH=/opt/public/${VERSION}/gcc-openblas-ompi
export PYPATH_MPI4PY="$(find ${PREREQ_PATH}/mpi4py-*/lib/python* -name site-packages)"
export PYTHONPATH="${PYPATH_MPI4PY}:${PYTHONPATH}"

# ensure to call the right python interpreter
PYDIR=$(. ${PREREQ_PATH}/*_mpi.sh ; echo $(dirname $(which python3)))
export PATH="${PYDIR}:${PATH}"

jobs=$(( ${NPROC_MAX:-8} * 3 / 4 ))
export ASTER_BUILD=debug

# mpi build
./configure --prefix=./mini --without-repo --no-enable-all
make install -j ${jobs}
make distclean

# sequential build
./configure --prefix=./mini --without-repo --no-enable-all --disable-mpi
make install -j ${jobs}
