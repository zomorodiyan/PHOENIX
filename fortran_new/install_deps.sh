#!/bin/bash
# Install dependencies for PHOENIX CFD: gfortran, python3, numpy, matplotlib, OpenMP
# Skips anything that is already installed.
set -e

# Detect package manager
if command -v apt-get >/dev/null 2>&1; then
    PM="apt"
elif command -v dnf >/dev/null 2>&1; then
    PM="dnf"
elif command -v yum >/dev/null 2>&1; then
    PM="yum"
elif command -v pacman >/dev/null 2>&1; then
    PM="pacman"
elif command -v zypper >/dev/null 2>&1; then
    PM="zypper"
else
    echo "Error: no supported package manager found (apt/dnf/yum/pacman/zypper)" >&2
    exit 1
fi

echo "Detected package manager: $PM"

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
    if command -v sudo >/dev/null 2>&1; then
        SUDO="sudo"
    else
        echo "Error: must run as root or have sudo installed" >&2
        exit 1
    fi
fi

# Check if a system package is installed
pkg_installed() {
    case "$PM" in
        apt)     dpkg -s "$1" >/dev/null 2>&1 ;;
        dnf|yum) rpm -q "$1" >/dev/null 2>&1 ;;
        pacman)  pacman -Qi "$1" >/dev/null 2>&1 ;;
        zypper)  rpm -q "$1" >/dev/null 2>&1 ;;
    esac
}

# Check if a python module is importable
py_has() {
    python3 -c "import $1" >/dev/null 2>&1
}

MISSING=()

add_if_missing() {
    if pkg_installed "$1"; then
        echo "  [skip] $1 already installed"
    else
        echo "  [add ] $1"
        MISSING+=("$1")
    fi
}

# Pick package names by distro
case "$PM" in
    apt)
        BASE_PKGS=(build-essential gfortran libgomp1 python3 python3-pip)
        NUMPY_PKG=python3-numpy
        MPL_PKG=python3-matplotlib
        ;;
    dnf|yum)
        BASE_PKGS=(gcc-gfortran libgomp python3 python3-pip)
        NUMPY_PKG=python3-numpy
        MPL_PKG=python3-matplotlib
        ;;
    pacman)
        BASE_PKGS=(gcc-fortran python python-pip)
        NUMPY_PKG=python-numpy
        MPL_PKG=python-matplotlib
        ;;
    zypper)
        BASE_PKGS=(gcc-fortran libgomp1 python3 python3-pip)
        NUMPY_PKG=python3-numpy
        MPL_PKG=python3-matplotlib
        ;;
esac

echo "Checking packages..."
for p in "${BASE_PKGS[@]}"; do
    add_if_missing "$p"
done

if py_has numpy; then
    echo "  [skip] $NUMPY_PKG (numpy already importable)"
else
    add_if_missing "$NUMPY_PKG"
fi

if py_has matplotlib; then
    echo "  [skip] $MPL_PKG (matplotlib already importable)"
else
    add_if_missing "$MPL_PKG"
fi

if [ ${#MISSING[@]} -eq 0 ]; then
    echo "Nothing to install — all dependencies already present."
else
    echo "Installing: ${MISSING[*]}"
    case "$PM" in
        apt)
            $SUDO apt-get update
            $SUDO apt-get install -y "${MISSING[@]}"
            ;;
        dnf)
            $SUDO dnf install -y "${MISSING[@]}"
            ;;
        yum)
            $SUDO yum install -y "${MISSING[@]}"
            ;;
        pacman)
            $SUDO pacman -Sy --noconfirm "${MISSING[@]}"
            ;;
        zypper)
            $SUDO zypper install -y "${MISSING[@]}"
            ;;
    esac
fi

echo
echo "Verifying installation..."
gfortran --version | head -n1
python3 --version
python3 -c "import numpy; print('numpy', numpy.__version__)"
python3 -c "import matplotlib; print('matplotlib', matplotlib.__version__)"
echo "OpenMP support:"
echo 'program t; use omp_lib; print *, "threads=", omp_get_max_threads(); end program' \
    | gfortran -fopenmp -x f95 - -o /tmp/_omp_test && /tmp/_omp_test && rm -f /tmp/_omp_test

echo
echo "All dependencies ready."
