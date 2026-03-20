#!/bin/bash -e

echo "+ compiling..."

opts=( "--prefix=./install" )
[ -f data-src/README ] && opts+=( "--with-data=data-src" )

if [ "${BUILDTYPE}" = "nightly-coverage" ]; then
    if [ "${ASTER_BUILD}" != "debug" ]; then
        echo "ERROR: ASTER_BUILD must be set as 'debug'"
        exit 4
    fi
    opts+=( "--coverage" )
fi
if [ "${BUILDTYPE}" = "nightly-sanitize" ]; then
    if [ "${ASTER_BUILD}" != "debug" ]; then
        echo "ERROR: ASTER_BUILD must be set as 'debug'"
        exit 4
    fi
    opts+=( "--enable-asan" )
fi

jobs=$(( ${NPROC_MAX:-8} * 3 / 4 ))

if [ "${OSNAME}" != "win" ]; then
    ./configure "${opts[@]}"
    make install -j ${jobs}
else
    ./waf_std configure --mingw-cross-compilation "${opts[@]}"
    ./waf_std install -j ${jobs}

    ./data/post_build_win.sh ./install
fi
