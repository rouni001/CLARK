#!/usr/bin/env bash
# One-time fix so scripts/rapl_energy.py (and batch-classify/run_all.sh -e)
# can read Intel RAPL's energy counters without sudo on every run.

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: sudo scripts/setup_rapl_permissions.sh [--dry-run]

Since a kernel fix for CVE-2020-8694 (a RAPL-based power side-channel
attack), most distro kernels ship
/sys/class/powercap/intel-rapl:*/energy_uj readable by root only, so
reading it (directly, or via scripts/rapl_energy.py) needs sudo every
time. This script installs a udev rule that makes those files
world-readable again every time the RAPL devices are (re-)enumerated
(boot, or a module reload), and applies the same permission immediately
to whatever RAPL devices already exist right now -- so after running this
once (with sudo), no future run needs sudo.

Options:
  --dry-run          print the udev rule and what would change, without
                      writing anything or requiring root.
  --rule-file <path> where to write the udev rule (default:
                      /etc/udev/rules.d/51-rapl-permissions.rules).
  -h, --help         show this message.
USAGE
}

die() {
	echo "Error: $*" >&2
	exit 1
}

DRY_RUN=0
RULE_FILE="/etc/udev/rules.d/51-rapl-permissions.rules"
SYSFS_GLOB="${RAPL_SETUP_SYSFS_GLOB:-/sys/class/powercap/intel-rapl:*/energy_uj}"

while [ "$#" -gt 0 ]; do
	case "$1" in
		--dry-run)
			DRY_RUN=1
			shift
			;;
		--rule-file)
			[ "$#" -ge 2 ] || die "--rule-file requires a path"
			RULE_FILE="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			die "unrecognized option: $1 (see -h)"
			;;
	esac
done

RULE_CONTENT='SUBSYSTEM=="powercap", KERNEL=="intel-rapl:*", ACTION=="add", RUN+="/bin/chmod -R a+r /sys%p"'

if [ "$DRY_RUN" -eq 1 ]; then
	echo "Would write $RULE_FILE with:"
	echo "  $RULE_CONTENT"
	echo "Would then chmod a+r any existing $SYSFS_GLOB files and reload udev rules."
	exit 0
fi

[ "$(id -u)" -eq 0 ] || die "this needs root to write a udev rule and change /sys permissions -- re-run with sudo: sudo $0"

mkdir -p "$(dirname "$RULE_FILE")"
printf '%s\n' "$RULE_CONTENT" > "$RULE_FILE"
echo "Wrote $RULE_FILE"

shopt -s nullglob
matched=0
for f in $SYSFS_GLOB; do
	chmod a+r "$f" 2>/dev/null && matched=$((matched + 1))
done
shopt -u nullglob
echo "Applied a+r to $matched existing RAPL energy_uj file(s) for this boot."

if command -v udevadm >/dev/null 2>&1; then
	udevadm control --reload-rules 2>/dev/null || true
	udevadm trigger --subsystem-match=powercap 2>/dev/null || true
	echo "Reloaded udev rules."
else
	echo "Note: udevadm not found here; the rule will still take effect on the next boot on a normal system with udev running."
fi

echo "Done. RAPL energy_uj should now be readable without sudo (verify as a normal, non-root user: cat /sys/class/powercap/intel-rapl:0/energy_uj)."
