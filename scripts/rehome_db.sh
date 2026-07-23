#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/rehome_db.sh <db-directory> [--dry-run]

Repairs the absolute-path metadata scripts/set_targets.sh cached under
<db-directory> (e.g. .protozoa, .protozoa.fileToAccssnTaxID,
.protozoa.fileToTaxIDs) after that directory has been copied or moved to
a different machine or path. Run this once right after copying the
database directory, then re-run scripts/set_targets.sh as usual (no
network access needed -- it will regenerate a clean targets.txt from the
repaired cache).

Options:
  --dry-run   report what would change without writing anything.
USAGE
}

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

die() {
	echo "Error: $*" >&2
	exit 1
}

if [ "$#" -lt 1 ]; then
	usage
	exit 1
fi

SCRIPT_DIR="$(script_dir)"
PYTHON_CMD="${CLARK_PYTHON:-python3}"

DBDR=""
EXTRA_ARGS=()

while [ "$#" -gt 0 ]; do
	case "$1" in
		--dry-run)
			EXTRA_ARGS+=(--dry-run)
			shift
			;;
		--help|-h)
			usage
			exit 0
			;;
		*)
			[ -z "$DBDR" ] || die "unrecognized extra argument: $1"
			DBDR="$1"
			shift
			;;
	esac
done

[ -n "$DBDR" ] || die "<db-directory> is required"
[ -d "$DBDR" ] || die "database directory '$DBDR' does not exist"
command -v "$PYTHON_CMD" >/dev/null 2>&1 || die "Python 3 is required"

"$PYTHON_CMD" "$SCRIPT_DIR/rehome_db.py" "$DBDR" "${EXTRA_ARGS[@]}"
