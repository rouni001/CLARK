#!/usr/bin/env python3
"""Manifest-driven RefSeq downloader for CLARK."""

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import random
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


USER_AGENT = "CLARK-refseq-downloader/1.4.5"
GZIP_ERRORS = (EOFError, OSError) + ((gzip.BadGzipFile,) if hasattr(gzip, "BadGzipFile") else ())


def parse_args():
    parser = argparse.ArgumentParser(description="Download CLARK RefSeq files from a manifest.")
    parser.add_argument("--download-list", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--timestamp", required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--failed-list", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--resume", type=int, choices=(0, 1), default=1)
    parser.add_argument("--attempts", type=int, default=5)
    parser.add_argument("--first-pass-attempts", type=int, default=3)
    parser.add_argument("--retry-delay", type=int, default=2)
    return parser.parse_args()


def setup_logging(path):
    logging.basicConfig(
        filename=path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )


def read_download_list(path):
    entries = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError("malformed download entry at line %d" % line_number)
            checksum = fields[2] if len(fields) > 2 and fields[2] else ""
            entries.append(
                {
                    "index": len(entries) + 1,
                    "source": fields[0],
                    "url": fields[1],
                    "checksum": checksum,
                }
            )
    return entries


def filename_from_url(url):
    path = urllib.parse.urlsplit(url).path
    name = os.path.basename(path)
    if not name:
        raise ValueError("download URL has no filename: %s" % url)
    return name


def gzip_is_valid(path):
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except GZIP_ERRORS:
        return False


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checksum_is_valid(path, expected):
    if not expected:
        return True
    return md5sum(path).lower() == expected.lower()


def manifest_row(timestamp, database, action, source, url, local_path, status):
    return "\t".join([timestamp, database, action, source, url, local_path, status])


def state_record(entry, output_path, status, attempt_count, error=""):
    return {
        "index": entry["index"],
        "source": entry["source"],
        "url": entry["url"],
        "local_path": str(output_path),
        "status": status,
        "attempts": attempt_count,
        "error": error,
        "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def append_jsonl(path, records):
    with open(path, "a", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")


def download_once(url, partial_path, resume):
    headers = {"User-Agent": USER_AGENT}
    mode = "wb"
    existing_size = partial_path.stat().st_size if partial_path.exists() else 0
    scheme = urllib.parse.urlsplit(url).scheme

    if resume and existing_size > 0 and scheme in ("http", "https"):
        headers["Range"] = "bytes=%d-" % existing_size
        mode = "ab"

    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        status = getattr(response, "status", None)
        if mode == "ab" and status != 206:
            mode = "wb"
        with open(partial_path, mode) as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)


def fetch_entry(entry, data_dir, max_attempts, failure_status, retry_delay, resume):
    url = entry["url"]
    source = entry["source"]
    output = data_dir / filename_from_url(url)
    decompressed = data_dir / output.name[:-3] if output.name.endswith(".gz") else None

    if output.exists():
        if gzip_is_valid(output) and checksum_is_valid(output, entry["checksum"]):
            row = manifest_row("", "", "download", source, url, str(output), "skipped-existing")
            return {"entry": entry, "status": "skipped-existing", "row": row, "path": output, "attempts": 0}
        output.unlink()

    if decompressed is not None and decompressed.exists() and decompressed.stat().st_size > 0:
        row = manifest_row("", "", "download", source, url, str(decompressed), "skipped-existing")
        return {"entry": entry, "status": "skipped-existing", "row": row, "path": decompressed, "attempts": 0}

    partial = data_dir / (output.name + ".part")
    if not resume and partial.exists():
        partial.unlink()

    last_error = ""
    for attempt in range(1, max_attempts + 1):
        try:
            download_once(url, partial, resume)
            if gzip_is_valid(partial) and checksum_is_valid(partial, entry["checksum"]):
                partial.replace(output)
                row = manifest_row("", "", "download", source, url, str(output), "downloaded")
                return {"entry": entry, "status": "downloaded", "row": row, "path": output, "attempts": attempt}
            last_error = "downloaded file failed gzip or checksum validation"
        except (OSError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last_error = str(exc)
        if partial.exists():
            partial.unlink()
        if attempt < max_attempts and retry_delay > 0:
            time.sleep(retry_delay * attempt + random.uniform(0, 0.25))

    row = manifest_row("", "", "download", source, url, str(output), failure_status)
    return {
        "entry": entry,
        "status": failure_status,
        "row": row,
        "path": output,
        "attempts": max_attempts,
        "error": last_error,
    }


def print_progress(label, completed, total):
    percent = 100 if total == 0 else int(completed * 100 / total)
    message = "%s: %d/%d files complete (%d%%)." % (label, completed, total, percent)
    if sys.stdout.isatty():
        sys.stdout.write("\r%s%s" % (message, "\n" if completed >= total else ""))
        sys.stdout.flush()
    else:
        print(message, flush=True)


def run_pass(entries, label, threads, max_attempts, failure_status, retry_delay, resume, data_dir):
    if not entries:
        return []

    results = []
    completed = 0
    total = len(entries)
    progress_step = max(1, min(total // 100, 100))
    next_report = progress_step
    print_progress(label, 0, total)

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(fetch_entry, entry, data_dir, max_attempts, failure_status, retry_delay, resume)
            for entry in entries
        ]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            completed += 1
            if completed >= next_report or completed == total:
                print_progress(label, completed, total)
                while next_report <= completed:
                    next_report += progress_step

    results.sort(key=lambda item: item["entry"]["index"])
    return results


def write_manifest_rows(path, rows, timestamp, database):
    with open(path, "a", encoding="utf-8") as handle:
        for row in rows:
            fields = row.split("\t")
            fields[0] = timestamp
            fields[1] = database
            handle.write("\t".join(fields) + "\n")


def write_failed_list(path, failures):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("source\turl\tlocal_path\tattempts\terror\n")
        for item in failures:
            handle.write(
                "%s\t%s\t%s\t%s\t%s\n"
                % (
                    item["entry"]["source"],
                    item["entry"]["url"],
                    item["path"],
                    item["attempts"],
                    item.get("error", ""),
                )
            )


def main():
    args = parse_args()
    if args.threads <= 0:
        raise SystemExit("--threads must be a positive integer")
    if args.attempts <= 0 or args.first_pass_attempts <= 0:
        raise SystemExit("attempt counts must be positive integers")

    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.log)
    entries = read_download_list(args.download_list)
    threads = min(args.threads, max(1, len(entries)))

    logging.info("Starting RefSeq download: database=%s files=%d threads=%d", args.database, len(entries), threads)
    first_results = run_pass(
        entries,
        "RefSeq download progress",
        threads,
        args.first_pass_attempts,
        "deferred",
        args.retry_delay,
        bool(args.resume),
        data_dir,
    )
    deferred = [item["entry"] for item in first_results if item["status"] == "deferred"]
    retry_results = []
    if deferred:
        print("Retrying %d deferred RefSeq download(s) after the first pass." % len(deferred), flush=True)
        retry_results = run_pass(
            deferred,
            "RefSeq deferred retry progress",
            threads,
            args.attempts,
            "failed",
            args.retry_delay,
            bool(args.resume),
            data_dir,
        )

    all_results = first_results + retry_results
    append_jsonl(
        args.state,
        [state_record(item["entry"], item["path"], item["status"], item["attempts"], item.get("error", "")) for item in all_results],
    )
    write_manifest_rows(args.manifest, [item["row"] for item in all_results], args.timestamp, args.database)

    failures = [item for item in all_results if item["status"] == "failed"]
    if failures:
        write_failed_list(args.failed_list, failures)
        print("Failed to download %d RefSeq genome file(s) after deferred retries." % len(failures), file=sys.stderr)
        print("First failed URLs:", file=sys.stderr)
        for item in failures[:20]:
            print("  [%s] %s" % (item["entry"]["source"], item["entry"]["url"]), file=sys.stderr)
        return 1

    logging.info("Finished RefSeq download: database=%s files=%d", args.database, len(entries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
