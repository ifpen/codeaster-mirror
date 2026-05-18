# This file set the environment for code_aster.
# Configuration for cronos mpi
. $(readlink -n -f $(dirname ${BASH_SOURCE:-${(%):-%x}}))/version.sh
. /software/restricted/simumeca/aster/prerequisites/${VERSION}/gcc-mkl-ompi-py311/cronos_mpi.sh
