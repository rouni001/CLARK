#!/usr/bin/env python3
"""Measure CPU package (and, if exposed, DRAM) energy via Intel RAPL around a wrapped command.

Reads RAPL domains under /sys/class/powercap/intel-rapl:* before and after
running the given command, sums the energy deltas separately for the
"package" domains (one per CPU socket) and the "dram" domain(s) if the
platform exposes them, and prints one ENERGY_PROFILE line with the
results. This is the same underlying MSR data Intel's PCM tool reports
(MSR_PKG_ENERGY_STATUS / MSR_DRAM_ENERGY_STATUS) via a different interface
-- reading it straight from sysfs avoids needing root + the msr kernel
module that PCM requires.

DRAM availability varies by platform and kernel: most server/Xeon
platforms expose it as its own RAPL domain (named exactly "dram"), either
nested under a package (e.g. intel-rapl:0:2) or as a sibling top-level
zone, depending on CPU generation and kernel version; many client/desktop
CPUs don't expose it via RAPL at all. This script looks in both places and
reports dram_domains=0 / dram_joules=0.0 plainly when it isn't available,
rather than silently omitting the field.

Domains are classified purely by their "name" file, not by nesting depth,
so both layouts are handled the same way:
  - name starting with "package" -> summed into pkg_joules (one entry per
    socket; this is the whole package including cores+uncore, so its own
    "core"/"uncore" sub-domains are NOT separately summed here to avoid
    double-counting).
  - name exactly "dram" -> summed into dram_joules.
Anything else (e.g. "core", "uncore", "psys") is ignored.

Neither figure is whole-system ("wall") power -- package excludes PSU
losses, fans, disks, and (where dram is a separate domain) DRAM itself;
dram excludes DIMM voltage-regulator losses. There is no universal
correction factor for either gap; as a rough, clearly-labeled second
estimate, this script also reports each raw figure (and their sum) scaled
down by a fixed fraction (default 30%) -- treat all of these as
component-level estimates, not a measurement of true wall power.

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

DOMAIN_DIR_RE = re.compile(r"^intel-rapl:[\d:]+$")


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sysfs-dir", default="/sys/class/powercap", help="powercap sysfs root (default: /sys/class/powercap)")
    parser.add_argument("--scale", type=float, default=0.70, help="fraction of each raw energy figure to also report (default: 0.70, i.e. a 30%% scale-down)")
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


def discover_domains(sysfs_dir):
    """Return {"package": [...], "dram": [...]}, each a list of
    (name, energy_uj_path, max_energy_range_uj_or_None) tuples, by
    scanning every intel-rapl:* entry (top-level or nested subdomain)
    under sysfs_dir and classifying it by its "name" file."""
    domains = {"package": [], "dram": []}
    try:
        entries = sorted(os.listdir(sysfs_dir))
    except OSError as exc:
        print("Warning: RAPL energy accounting unavailable: cannot list %s (%s)" % (sysfs_dir, exc), file=sys.stderr)
        return domains

    for entry in entries:
        if not DOMAIN_DIR_RE.match(entry):
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

        if name.startswith("package"):
            group = "package"
        elif name == "dram":
            group = "dram"
        else:
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
        domains[group].append((name, energy_path, max_range))
    return domains


def read_total_uj(domains):
    readings = {}
    for name, energy_path, _max_range in domains:
        readings[energy_path] = read_int(energy_path)
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

    domains = discover_domains(args.sysfs_dir)
    pkg_domains = domains["package"]
    dram_domains = domains["dram"]
    all_domains = pkg_domains + dram_domains

    if not pkg_domains:
        print("Warning: no readable RAPL package domains found under %s; running without energy measurement." % args.sysfs_dir, file=sys.stderr)
        return subprocess.call(args.command)

    before = read_total_uj(all_domains)
    start = time.monotonic()
    exit_code = subprocess.call(args.command)
    elapsed = time.monotonic() - start
    after = read_total_uj(all_domains)

    pkg_uj = delta_with_wraparound(before, after, pkg_domains)
    dram_uj = delta_with_wraparound(before, after, dram_domains)

    pkg_joules = pkg_uj / 1e6
    dram_joules = dram_uj / 1e6
    total_joules = pkg_joules + dram_joules

    print(
        "ENERGY_PROFILE unit=joules pkg_joules=%.6f pkg_joules_scaled=%.6f "
        "dram_joules=%.6f dram_joules_scaled=%.6f "
        "total_joules=%.6f total_joules_scaled=%.6f "
        "scale=%.2f elapsed_s=%.6f pkg_domains=%d dram_domains=%d"
        % (
            pkg_joules, pkg_joules * args.scale,
            dram_joules, dram_joules * args.scale,
            total_joules, total_joules * args.scale,
            args.scale, elapsed, len(pkg_domains), len(dram_domains),
        )
    )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
