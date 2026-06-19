#!/usr/bin/env bash

set -euo pipefail

script_dir() {
	local source="${BASH_SOURCE[0]}"
	while [ -h "$source" ]; do
		local dir
		dir="$(cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd)"
		source="$(readlink "$source")"
		[[ "$source" != /* ]] && source="$dir/$source"
	done
	cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd
}

LDIR="$(script_dir)"

echo "Installing CLARK from $LDIR"
echo "Compiler: ${CXX:-c++}"

if ! command -v "${CXX:-c++}" >/dev/null 2>&1; then
	echo "Error: C++ compiler '${CXX:-c++}' was not found." >&2
	echo "Install g++ or set CXX to a working compiler, then rerun ./install.sh." >&2
	exit 1
fi

if ! command -v make >/dev/null 2>&1; then
	echo "Error: make was not found." >&2
	exit 1
fi

make -C "$LDIR" all

echo
echo "CLARK executables are available in: $LDIR/exe"
echo "Run 'make test' to verify the installation."
