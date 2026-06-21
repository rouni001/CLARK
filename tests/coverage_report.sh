#!/usr/bin/env bash

set -euo pipefail

REPO_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
COVERAGE_DIR="$REPO_DIR/build/coverage"
GCOV_OUT="$COVERAGE_DIR/gcov.out"
MIN_COVERAGE="${COVERAGE_MIN:-10}"

if [ -z "${COVERAGE_CXXFLAGS+x}" ]; then
	COVERAGE_CXXFLAGS="-O0 -g --coverage"
fi

if [ -z "${COVERAGE_LDFLAGS+x}" ]; then
	COVERAGE_LDFLAGS="--coverage"
fi

command -v gcov >/dev/null 2>&1 || {
	echo "gcov is required to generate coverage reports." >&2
	exit 1
}

make -C "$REPO_DIR" clean
make -C "$REPO_DIR" all unit-tests CXXFLAGS="$COVERAGE_CXXFLAGS" LDFLAGS="$COVERAGE_LDFLAGS"
"$REPO_DIR/tests/run_tests.sh"

mkdir -p "$COVERAGE_DIR"
: > "$GCOV_OUT"

while IFS= read -r -d '' gcno; do
	gcov -t "$gcno" >> "$GCOV_OUT"
done < <(find "$REPO_DIR/exe" -name '*.gcno' -print0 | sort -z)

awk -v repo="$REPO_DIR" -v min="$MIN_COVERAGE" '
function normalize(path) {
	sub(/^\.\//, "", path)
	if (index(path, repo "/") == 1) {
		path = substr(path, length(repo) + 2)
	}
	if (path ~ /^build\/(default|light|spaced)\//) {
		sub(/^build\/(default|light|spaced)\//, "src/", path)
	}
	return path
}

function wanted(path) {
	return path ~ /^src\/.*\.(cc|hh)$/
}

/^[[:space:]]*-:[[:space:]]*0:Source:/ {
	src = $0
	sub(/^[[:space:]]*-:[[:space:]]*0:Source:/, "", src)
	src = normalize(src)
	include = wanted(src)
	next
}

include && /^[[:space:]]*([0-9]+|#####):[[:space:]]*[0-9]+:/ {
	count = $1
	line = $2
	gsub(/:/, "", count)
	gsub(/:/, "", line)
	key = src ":" line

	if (!(key in seen)) {
		seen[key] = 1
		total++
		file_total[src]++
	}
	if (count != "#####" && count + 0 > 0 && !(key in covered)) {
		covered[key] = 1
		covered_total++
		file_covered[src]++
	}
}

END {
	percent = total ? covered_total * 100 / total : 0
	printf("Coverage summary for compiled CLARK C++ sources\n")
	printf("Covered lines: %d\n", covered_total)
	printf("Executable lines: %d\n", total)
	printf("Line coverage: %.2f%%\n", percent)
	printf("Required minimum: %.2f%%\n", min + 0)
	printf("\nPer-file coverage:\n")
	for (file in file_total) {
		file_percent = file_total[file] ? (file_covered[file] + 0) * 100 / file_total[file] : 0
		printf("%7.2f%% %5d/%-5d %s\n", file_percent, file_covered[file] + 0, file_total[file], file)
	}
	if (percent + 0.0001 < min + 0) {
		exit 2
	}
}
' "$GCOV_OUT"
