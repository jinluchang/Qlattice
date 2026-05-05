#!/usr/bin/env python3

# Usage:
# monitor-slurm-jobs.py [--check] [--timeout MINUTES] [--interval MINUTES] [--dir DIR]
#
# Monitor running SLURM jobs and kill those whose output files (slurm-$id.out)
# have not been updated for a specified time.
#
# Options:
#   --check           Only report stale jobs, do not kill them.
#   --timeout MINUTES Minutes of inactivity before killing a job (default: 15).
#   --interval MINUTES Minutes between checks in continuous mode (default: 1).
#   --once            Run one check and exit (do not loop).
#   --dir DIR         Directory to search for slurm-*.out files (default: cwd).

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print(
        "You are using not supported Python {}.{}.".format(
            sys.version_info.major, sys.version_info.minor
        )
    )
    sys.exit(1)

import os
import time
import glob
import subprocess
import argparse

def get_running_jobs():
    """Return a dict of {job_id: (user, partition, name, state, time_str)} for the current user's running jobs."""
    username = os.environ.get("USER") or os.getlogin()
    try:
        result = subprocess.run(
            ["squeue", "-u", username, "-h", "-o", "%i %u %P %j %T %M"],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        print("ERROR: squeue command not found. Is SLURM installed?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: squeue failed: {e.stderr.strip()}")
        sys.exit(1)
    jobs = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split()
        if len(parts) < 6:
            continue
        job_id, user, partition, name, state, time_str = (
            parts[0],
            parts[1],
            parts[2],
            parts[3],
            parts[4],
            parts[5],
        )
        if state != "PENDING":
            jobs[job_id] = (user, partition, name, state, time_str)
    return jobs

def find_output_file(job_id, search_dir):
    """Find the output file for a given job_id in search_dir."""
    pattern = os.path.join(search_dir, f"slurm-{job_id}.out")
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    return None

def file_age_seconds(filepath):
    """Return the number of seconds since the file was last modified."""
    try:
        mtime = os.path.getmtime(filepath)
    except OSError:
        return None
    return time.time() - mtime

def cancel_job(job_id):
    """Cancel a SLURM job."""
    try:
        subprocess.run(
            ["scancel", job_id],
            capture_output=True,
            text=True,
            check=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: failed to cancel job {job_id}: {e.stderr.strip()}")
        return False

def format_age(seconds):
    """Format seconds into a human-readable string in minutes."""
    if seconds is None:
        return "unknown"
    return f"{seconds / 60:.1f}m"

def monitor_once(search_dir, timeout_seconds, check_only):
    """Check all running jobs and kill stale ones. Returns list of stale jobs."""
    jobs = get_running_jobs()
    if not jobs:
        print("No running SLURM jobs found.")
        return []
    stale_jobs = []
    now_str = time.strftime("%Y-%m-%d %H:%M:%S")
    for job_id, (user, partition, name, state, time_str) in jobs.items():
        outfile = find_output_file(job_id, search_dir)
        if outfile is None:
            print(
                f"[{now_str}] Job {job_id} ({name}): no output file found (slurm-{job_id}.out)"
            )
            continue
        age = file_age_seconds(outfile)
        if age is not None and age > timeout_seconds:
            age_str = format_age(age)
            stale_jobs.append((job_id, name, outfile, age_str))
            if check_only:
                print(
                    f"[{now_str}] STALE Job {job_id} ({name}): {outfile} inactive for {age_str}"
                )
            else:
                print(
                    f"[{now_str}] KILLING Job {job_id} ({name}): {outfile} inactive for {age_str}"
                )
                if cancel_job(job_id):
                    print(f"  Job {job_id} cancelled.")
        else:
            age_str = format_age(age)
            print(
                f"[{now_str}] OK Job {job_id} ({name}): {outfile} last updated {age_str} ago"
            )
    return stale_jobs

def main():
    parser = argparse.ArgumentParser(
        description="Monitor running SLURM jobs and kill those with stale output files.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only report stale jobs, do not kill them.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=15,
        help="Minutes of inactivity before killing a job (default: 15).",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=1,
        help="Minutes between checks in continuous mode (default: 1).",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Run one check and exit (do not loop).",
    )
    parser.add_argument(
        "--dir",
        type=str,
        default=".",
        help="Directory to search for slurm-*.out files (default: cwd).",
    )
    args = parser.parse_args()
    search_dir = os.path.abspath(args.dir)
    timeout_seconds = args.timeout * 60
    interval_seconds = args.interval * 60
    print(f"Monitoring SLURM jobs in {search_dir}")
    print(f"Timeout: {args.timeout} minutes")
    print(f"Mode: {'check only' if args.check else 'kill stale jobs'}")
    if not args.once:
        print(f"Check interval: {args.interval:.1f} minutes")
    print()
    if args.once:
        stale = monitor_once(search_dir, timeout_seconds, args.check)
        if stale:
            print(f"\n{len(stale)} stale job(s) found.")
        else:
            print("\nNo stale jobs found.")
        sys.exit(1 if stale and not args.check else 0)
    try:
        while True:
            stale = monitor_once(search_dir, timeout_seconds, args.check)
            if stale:
                print(f"\n{len(stale)} stale job(s) found.")
            else:
                print("\nNo stale jobs found.")
            print(f"Next check in {args.interval:.1f} minutes...\n")
            time.sleep(interval_seconds)
    except KeyboardInterrupt:
        print("\nStopped.")

if __name__ == "__main__":
    main()
