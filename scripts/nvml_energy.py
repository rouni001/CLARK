#!/usr/bin/env python3
"""Measure NVIDIA GPU energy (NVML) around a wrapped command.

Reads NVML's per-device energy counter before and after running the given
command, sums the deltas across the selected GPU(s), and prints one
ENERGY_PROFILE_GPU line. No changes to cuCLARK (or any other GPU program)
are needed -- this wraps the process from the outside exactly like
rapl_energy.py does for CPU package energy, just reading from NVML
instead of RAPL sysfs.

Primary method: nvmlDeviceGetTotalEnergyConsumption(), a cumulative
millijoule counter most modern (Volta+) NVIDIA GPUs support directly --
no polling, no integration error. If a GPU/driver doesn't support it,
this falls back to sampling nvmlDeviceGetPowerUsage() (instantaneous
milliwatts) on a background thread at --sample-interval and integrating
power over time (trapezoidal). The fallback is inherently approximate --
prefer the counter method wherever it's available (`method=counter` vs.
`method=sampling` in the output tells you which one ran).

By default every GPU NVML can see is included; pass --devices 0,1 to
restrict to specific NVML device indices (e.g. to match CUDA_VISIBLE_DEVICES
when only some GPUs on the box are actually used by the wrapped command).

NVML GPU energy is closer to true board power than CPU package RAPL is to
system wall power (on most server GPUs it already includes onboard HBM/
GDDR), but it still excludes the host machine's own consumption and PSU
AC-DC conversion losses. For consistency with rapl_energy.py's reporting,
this script also prints the raw figure scaled down by the same default
30% (--scale) as a rough second estimate.

Usage:
  scripts/nvml_energy.py [--devices 0,1] [--scale FRACTION] [--sample-interval SECONDS] -- <command> [args...]

Exits with the wrapped command's exit code. If pynvml isn't installed, NVML
can't be initialized, or no matching GPU is found, prints a warning to
stderr and still runs the command, just without an ENERGY_PROFILE_GPU line.
"""

import argparse
import subprocess
import sys
import threading
import time


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--devices", default="", help="comma-separated NVML device indices to measure (default: all devices NVML finds)")
    parser.add_argument("--scale", type=float, default=0.70, help="fraction of the raw GPU energy to also report (default: 0.70, i.e. a 30%% scale-down)")
    parser.add_argument("--sample-interval", type=float, default=0.05, help="seconds between power samples for the sampling fallback (default: 0.05)")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="command to run, after --")
    args = parser.parse_args(argv)
    if not args.command:
        parser.error("no command given (usage: nvml_energy.py [options] -- <command> [args...])")
    if args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("no command given (usage: nvml_energy.py [options] -- <command> [args...])")
    return args


def discover_handles(pynvml, requested):
    count = pynvml.nvmlDeviceGetCount()
    indices = requested if requested else list(range(count))
    handles = []
    for i in indices:
        if i >= count:
            print("Warning: NVML device index %d does not exist (found %d device(s)); skipping" % (i, count), file=sys.stderr)
            continue
        handles.append(pynvml.nvmlDeviceGetHandleByIndex(i))
    return handles


def supports_total_energy(pynvml, handle):
    try:
        pynvml.nvmlDeviceGetTotalEnergyConsumption(handle)
        return True
    except Exception:
        return False


class PowerSampler(threading.Thread):
    """Fallback for GPUs without a total-energy counter: poll
    instantaneous power draw and trapezoidally integrate it over time."""

    def __init__(self, pynvml, handles, interval):
        super().__init__()
        self._pynvml = pynvml
        self._handles = handles
        self._interval = interval
        self._stop_event = threading.Event()
        self.energy_mj = 0.0

    def _read_power_mw(self, handle):
        try:
            return self._pynvml.nvmlDeviceGetPowerUsage(handle)
        except Exception:
            return 0.0

    def run(self):
        last_time = time.monotonic()
        last_power = [self._read_power_mw(h) for h in self._handles]
        while not self._stop_event.wait(self._interval):
            now = time.monotonic()
            dt = now - last_time
            power = [self._read_power_mw(h) for h in self._handles]
            for lp, p in zip(last_power, power):
                # milliwatts * seconds = millijoules
                self.energy_mj += (lp + p) / 2.0 * dt
            last_power = power
            last_time = now

    def stop(self):
        self._stop_event.set()


def measure_with_counter(pynvml, handles, command):
    before = [pynvml.nvmlDeviceGetTotalEnergyConsumption(h) for h in handles]
    start = time.monotonic()
    exit_code = subprocess.call(command)
    elapsed = time.monotonic() - start
    after = [pynvml.nvmlDeviceGetTotalEnergyConsumption(h) for h in handles]
    delta_mj = sum(a - b for a, b in zip(after, before))
    return exit_code, delta_mj, elapsed


def measure_with_sampling(pynvml, handles, command, interval):
    sampler = PowerSampler(pynvml, handles, interval)
    sampler.start()
    start = time.monotonic()
    exit_code = subprocess.call(command)
    elapsed = time.monotonic() - start
    sampler.stop()
    sampler.join()
    return exit_code, sampler.energy_mj, elapsed


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    requested = [int(x) for x in args.devices.split(",") if x.strip() != ""] if args.devices else []

    try:
        import pynvml
    except ImportError:
        print("Warning: pynvml is not installed (pip install nvidia-ml-py); running without GPU energy measurement.", file=sys.stderr)
        return subprocess.call(args.command)

    try:
        pynvml.nvmlInit()
    except Exception as exc:
        print("Warning: NVML could not be initialized (%s); running without GPU energy measurement." % exc, file=sys.stderr)
        return subprocess.call(args.command)

    try:
        handles = discover_handles(pynvml, requested)
        if not handles:
            print("Warning: no NVIDIA GPU found via NVML; running without GPU energy measurement.", file=sys.stderr)
            return subprocess.call(args.command)

        if supports_total_energy(pynvml, handles[0]):
            method = "counter"
            exit_code, delta_mj, elapsed = measure_with_counter(pynvml, handles, args.command)
        else:
            method = "sampling"
            exit_code, delta_mj, elapsed = measure_with_sampling(pynvml, handles, args.command, args.sample_interval)

        joules = delta_mj / 1000.0
        print(
            "ENERGY_PROFILE_GPU unit=joules gpu_joules=%.6f gpu_joules_scaled=%.6f scale=%.2f elapsed_s=%.6f devices=%d method=%s"
            % (joules, joules * args.scale, args.scale, elapsed, len(handles), method)
        )
        return exit_code
    finally:
        pynvml.nvmlShutdown()


if __name__ == "__main__":
    sys.exit(main())
