#!/bin/bash
# The memory check of the hypothesis tests. The program test/memory.goose
# runs one call on the largest input of the package. This script measures
# the peak size of that program and compares it with the size of its input.
# A call must stay under 8 times its input.
#
# Run it from the package directory:
#   bash test/memory.sh
# The script needs /usr/bin/time. It skips the check when the tool is
# missing.

set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
root="$(dirname "$here")"
# A private directory keeps another user from planting a link at a fixed name.
work="$(mktemp -d "${TMPDIR:-/tmp}/statistics-memory.XXXXXX")" || exit 1
trap 'rm -rf "$work"' EXIT
bin="$work/memory"
max_input=4194304
input_bytes=$((max_input * 8))
limit=$((input_bytes * 8))

cd "$root" || exit 1
goose build -o "$bin" test/memory.goose || exit 1

peak_of() {
    # The peak size of the process in bytes.
    if [ "$(uname -s)" = "Darwin" ]; then
        /usr/bin/time -l "$@" 2>&1 >/dev/null | awk '/maximum resident set size/ {print $1}'
    else
        /usr/bin/time -f "%M" "$@" 2>&1 >/dev/null | tail -1 | awk '{print $1 * 1024}'
    fi
}

full_base=$(peak_of env N=$max_input OP=0 "$bin")
if [ -z "$full_base" ]; then
    echo "memory check skipped, the tool gives no peak size"
    exit 0
fi
half_base=$(peak_of env N=$((max_input / 2)) OP=0 "$bin")
echo "baseline $full_base bytes, half $half_base bytes"

bad=0
for op in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27; do
    # The calls that read two samples hold half of the input in each one.
    n=$max_input
    base=$full_base
    case "$op" in
        2|5|6|7|8|10|11|13|14|17|18|19|20|27) n=$((max_input / 2)); base=$half_base ;;
    esac
    peak=$(peak_of env N=$n OP=$op "$bin")
    extra=$((peak - base))
    echo "call $op peak $peak bytes, extra $extra bytes"
    if [ "$extra" -gt "$limit" ]; then
        echo "FAIL call $op holds $extra bytes, the limit is $limit bytes"
        bad=$((bad + 1))
    fi
done

if [ "$bad" -gt 0 ]; then
    echo "$bad calls hold too much memory"
    exit 1
fi
echo "every call stays under 8 times its input"
