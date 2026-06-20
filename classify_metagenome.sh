#!/bin/sh
SCRIPT_DIR=$(CDPATH= cd "$(dirname "$0")" && pwd -P)
CLARK_HOME="${CLARK_HOME:-$SCRIPT_DIR}"
export CLARK_HOME
exec "$SCRIPT_DIR/scripts/classify_metagenome.sh" "$@"
