#!/usr/bin/env bash
set -eu

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
binary="${RTK_BINARY:-$repo_root/build/RTK}"
duration="${1:-4h}"
if [ "$#" -gt 0 ]; then shift; fi
if [ ! -x "$binary" ]; then
    printf 'RTK executable not found: %s\n' "$binary" >&2
    exit 1
fi
binary="$(realpath "$binary")"
mkdir -p "$repo_root/runs"
run_dir="$(mktemp -d "$repo_root/runs/rooftop_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")"
printf 'Run directory: %s\n' "$run_dir"
cd "$run_dir"

{
    printf 'Started UTC: %s\n' "$(date -u +%FT%TZ)"
    printf 'Duration: %s\nCommand: ' "$duration"
    printf '%q ' "$binary" --stream-rooftop "$@"
    printf '\n'
} >run.log
# 保留 RTK 的退出码；正常响应计时器的 SIGTERM 时返回 0，超时强杀仍为非零。
if timeout --preserve-status --signal=TERM --kill-after=30s "$duration" \
    "$binary" --stream-rooftop "$@" >>run.log 2>&1; then
    status=0
else
    status=$?
fi
printf '%s\n' "$status" >exit_status.txt
printf 'Ended UTC: %s\nExit status: %s\n' "$(date -u +%FT%TZ)" "$status" >>run.log
printf 'Exit status: %s; log: %s/run.log\n' "$status" "$run_dir"
exit "$status"
