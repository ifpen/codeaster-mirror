#!/bin/bash

prefix=$(dirname $(dirname $(readlink -n -f ${0})))

export ASTER_REQS_USE_DEBUG=${ASTER_REQS_USE_DEBUG:-0}

defaultdir=${SCRATCHDIR:-${HOME}}/.cache/codeaster-prerequisites
export ASTER_REQS_ROOT="${ASTER_REQS_ROOT:-${defaultdir}}"

_test() {
    printf "\nchecking for ${1}...\n"
    ${1} --version 2> /dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: please install ${1}."
        return 1
    fi
    return 0
}

is_intranet() {
    printf "checking for EDF network... "
    curl -m 2 https://gitlab.pleiade.edf.fr > /dev/null 2>&1
    okintranet=$?
    test ${okintranet} -eq 0 && echo ok || echo no
    return ${okintranet}
}

do_install() {
    if [ ${ASTER_REQS_CONFIRM:-n} != "y" ]; then
        printf "✨ Do you want to download the requirements ([y]/n, c to cancel)? "
        read answer
        answer=$(tr '[:upper:]' '[:lower:]' <<< "${answer}")
        [ "${answer}" = "n" ] && return 0
        [ "${answer}" = "c" ] && echo "exiting..." && return 1
    fi

    ASTER_REQS_PACKAGE=${ASTER_REQS_PACKAGE:-"gcc-ompi"}
    if [ "${ASTER_REQS_PACKAGE}" = "ask" ]; then
        echo
        echo "Select the package you want to download:"
        echo "  1. full (embedding gcc & openmpi, recommended)"
        echo "  2. with gcc, without openmpi (you must have the same version of openmpi on the host)"
        echo "  3. without gcc, without openmpi (you must have the same version of openmpi on the host)"
        printf "✨ your choice ([1]/2/3, 0 to cancel)? "
        read answer
        [ "${answer}" = "2" ] && ASTER_REQS_PACKAGE="gcc-noompi"
        [ "${answer}" = "3" ] && ASTER_REQS_PACKAGE="nogcc-noompi"
        [ "${answer}" = "0" ] && echo "exiting..." && return 1
    fi
    [ ${ASTER_REQS_USE_DEBUG} -eq 1 ] && VARIANT="-debug"
    arch="codeaster-prerequisites-${VERSION}-${ASTER_REQS_PACKAGE}${VARIANT}.sh"

    _test curl || return 1

    is_intranet
    okintranet=$?

    mkdir -p "${ASTER_REQS_ROOT}"
    echo
    if [ ${okintranet} -eq 0 ]; then
        echo "⏳ downloading requirements archive ${arch} (it may take a few minutes)..."
        curl -fsSL https://minio.retd.edf.fr/codeaster/prereq/${arch} \
            -o "${ASTER_REQS_ROOT}/${arch}"
    else
        echo "❌ downloading requirements from internet is not yet supported, sorry"
        return 1
    fi
    if [ ! -f "${ASTER_REQS_ROOT}/${arch}" ]; then
        echo "❌ download failed"
        return 1
    fi

    (
        cd "${ASTER_REQS_ROOT}"
        chmod 755 ${arch}
        ./${arch} && rm -f ${arch}
    )

    printf "✅ code_aster requirements installed\n\n"

    return 0
}

check_requirements_main() {
    args=( "-v" )
    [ ${ASTER_REQS_USE_DEBUG} -eq 1 ] && args=()

    found=$(find "${ASTER_REQS_ROOT}" -name "codeaster-prerequisites-${VERSION}-*" -type d 2> /dev/null | grep ${args[@]} debug)
    if [ -z "${found}" ]; then
        echo
        echo "code_aster requirements not found in '${ASTER_REQS_ROOT}/'."
        echo "If they are already installed, define ASTER_CONFIG environment variable" \
            "to the environment file."
        do_install || return 4
    fi
    found=$(find "${ASTER_REQS_ROOT}" -name "codeaster-prerequisites-${VERSION}-*" -type d 2> /dev/null | grep ${args[@]} debug)
    if [ -z "${found}" ]; then
        return 1
    fi
    if [ "${DEBUG}" = "1" ] && [ "${found}/VERSION" ]; then
        echo "package details:"
        cat "${found}/VERSION"
    fi
    export ASTER_CONFIG=$(ls ${found}/*.sh 2> /dev/null | head -1)
    return 0
}

check_requirements_main "$@"
