#!/usr/bin/env python3
"""Measure CPU package energy (Intel RAPL) around a wrapped command.

Reads the "package" RAPL domains under /sys/class/powercap/intel-rapl:*
(one per CPU socket) before and after running the given command, sums the
energy deltas across sockets, and prints one ENERGY_PROFILE line with the
result. This is the same underlying MSR data Intel's PCM tool reports
(MSR_PKG_ENERGY_STATUS) via a different interface -- reading it straight
from sysfs avoids needing root + the msr kernel module that PCM requires.

RAPL's package-energy counter measures the CPU package silicon only, not
whole-system ("wall") power -- it excludes PSU losses, fans, disks, and
(on server/Xeon platforms) DRAM, which RAPL reports as a separate domain
not summed here. There is no universal correction factor for this gap; as
a rough, clearly-labeled lower/upper bracket, this script also reports the
raw value scaled down by a fixed fraction (default 30%) alongside the raw
package figure -- treat both as component-level estimates, not a
measurement of true wall power.

Only top-level "package-*" domains are summed (e.g. intel-rapl:0,
intel-rapl:1 for a 2-socket machine); their "core"/"uncore" sub-domains
(intel-rapl:0:0, intel-rapl:0:1, ...) are already included in the package
total and are skipped to avoid double-counting.

Usage:
  scripts/rapl_energy.py [--sysfs-dir DIR] [--scale FRACTION] -- <command> [args...]

Exits with the wrapped command's exit code. If RAPL is unavailable or
unreadable (no /sys/class/powercap, permission denied, no package
domains), prints a warning to stderr and still runs the command, just
without an ENERGY_PROFILE line.
"""

import argparse
import os
import re
import subprocess
import sys
import time

PACKAGE_DIR_RE = re.compile(r"^intel-rapl:\d+$")


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sysfs-dir", default="/sys/class/powercap", help="powercap sysfs root (default: /sys/class/powercap)")
    parser.add_argument("--scale", type=float, default=0.70, help="fraction of the raw package energy to also report (default: 0.70, i.e. a 30%% scale-down)")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command to run, after --")
    args = parser.parse_args(argv)
    if not args.command:
        parser.error("no command given (usage: rapl_energy.py [options] -- <command> [args...])")
    if args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("no command given (usage: rapl_energy.py [options] -- <command> [args...])")
    return args


def read_int(path):
    with open(path, "r", encoding="utf-8") as handle:
        return int(handle.read().strip())


def discover_package_domains(sysfs_dir):
    """Return [(name, energy_uj_path, max_energy_range_uj_or_None), ...]
    for each top-level RAPL "package-*" domain."""
    domains = []
    try:
        entries = sorted(os.listdir(sysfs_dir))
    except OSError as exc:
        print("Warning: RAPL energy accounting unavailable: cannot list %s (%s)" % (sysfs_dir, exc), file=sys.stderr)
        return domains

    for entry in entries:
        if not PACKAGE_DIR_RE.match(entry):
            continue
        domain_dir = sysfs_dir + "/" + entry
        name_path = domain_dir + "/name"
        energy_path = domain_dir + "/energy_uj"
        range_path = domain_dir + "/max_energy_range_uj"
        try:
            with open(name_path, "r", encoding="utf-8") as handle:
                name = handle.read().strip()
        except OSError as exc:
            print("Warning: RAPL domain %s: cannot read name (%s), skipping" % (entry, exc), file=sys.stderr)
            continue
        if not name.startswith("package"):
            continue
        try:
            read_int(energy_path)
        except OSError as exc:
            print("Warning: RAPL domain %s (%s): cannot read energy_uj (%s), skipping" % (entry, name, exc), file=sys.stderr)
            continue
        max_range = None
        try:
            max_range = read_int(range_path)
        except OSError:
            pass
        domains.append((name, energy_path, max_range))
    return domains


def read_total_uj(domains):
    total = 0
    readings = {}
    for name, energy_path, _max_range in domains:
        value = read_int(energy_path)
        readings[energy_path] = value
        total += value
    return readings


def delta_with_wraparound(before, after, domains):
    total_delta = 0
    for name, energy_path, max_range in domains:
        b = before[energy_path]
        a = after[energy_path]
        if a >= b:
            total_delta += a - b
        elif max_range is not None:
            total_delta += (a - b) + max_range
        else:
            print("Warning: RAPL domain %s: counter appears to have wrapped and max_energy_range_uj is unavailable; skipping its contribution" % name, file=sys.stderr)
    return total_delta


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)

    domains = discover_package_domains(args.sysfs_dir)
    if not domains:
        print("Warning: no readable RAPL package domains found under %s; running without energy measurement." % args.sysfs_dir, file=sys.stderr)
        return subprocess.call(args.command)

    before = read_total_uj(domains)
    start = time.monotonic()
    exit_code = subprocess.call(args.command)
    elapsed = time.monotonic() - start
    after = read_total_uj(domains)

    delta_uj = delta_with_wraparound(before, after, domains)
    joules = delta_uj / 1e6
    scaled_joules = joules * args.scale

    print(
        "ENERGY_PROFILE unit=joules pkg_joules=%.6f pkg_joules_scaled=%.6f scale=%.2f elapsed_s=%.6f domains=%d"
        % (joules, scaled_joules, args.scale, elapsed, len(domains))
    )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
