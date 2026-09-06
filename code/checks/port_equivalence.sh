#!/usr/bin/env bash
# port_equivalence.sh — proves a file ported from analysis/ into code/ differs
# from its original ONLY in path handling.
#
#   usage: code/checks/port_equivalence.sh <new-file> <old-file>
#
# Two independent checks:
#   1. Body diff. Strip every line that is allowed to change (the header, local
#      ROOT/HELPERS definitions, includes, and the RAW/PROC/OUT constants) plus
#      blank lines, from BOTH files, then diff. The remainder must be identical.
#   2. Include set. Because check 1 strips include lines wholesale, it cannot see
#      a DROPPED include. So compare the set of included basenames separately.
#
# TEMPORARY: delete this together with analysis/.
set -uo pipefail

new="$1"; old="$2"

strip_paths() {
  grep -vE '^[[:space:]]*$' "$1" \
  | grep -vE '^[[:space:]]*using Pkg[[:space:]]*$' \
  | grep -vE '^[[:space:]]*Pkg\.activate\(' \
  | grep -vE '^[[:space:]]*const[[:space:]]+(ROOT|HELPERS|DATA|PLOTS)[[:space:]]*=' \
  | grep -vE '^[[:space:]]*isfile\(joinpath\(ROOT,' \
  | grep -vE '^[[:space:]]*error\("ROOT = ' \
  | grep -vE '^[[:space:]]*include\(joinpath\(' \
  | grep -vE '^[[:space:]]*const[[:space:]]+(RAW|PROC|OUT|OUTDIR)[[:space:]]*='
}

# Included basenames, sorted — order of includes may legitimately change.
#
# Select lines that call include(joinpath(...)), then pull out any quoted
# "*.jl" token from that line. Matching is deliberately NOT anchored to
# immediately follow "joinpath(": the analysis/ originals nest a dirname(...)
# call before the path segments (e.g. include(joinpath(dirname(@__DIR__),
# "helpers", "types.jl"))), and a naive "no ) between joinpath( and the
# string" pattern can't cross that inner ")" — it would silently find zero
# includes on the old side for the majority of real files.
include_set() {
  grep -E 'include\(joinpath\(' "$1" \
  | grep -oE '"[A-Za-z0-9_]+\.jl"' | tr -d '"' | sort
}

rc=0

if ! diff <(strip_paths "$old") <(strip_paths "$new") > /tmp/port_body_diff.$$; then
  echo "BODY DIFFERS: $new"
  echo "  (< is analysis/, > is code/)"
  sed 's/^/  /' /tmp/port_body_diff.$$
  rc=1
fi
rm -f /tmp/port_body_diff.$$

# paths.jl is included by the new file and by nothing in the old — expected.
if ! diff <(include_set "$old") \
          <(include_set "$new" | grep -v '^paths\.jl$') > /tmp/port_inc_diff.$$; then
  echo "INCLUDE SET DIFFERS: $new"
  echo "  (< is analysis/, > is code/)"
  sed 's/^/  /' /tmp/port_inc_diff.$$
  rc=1
fi
rm -f /tmp/port_inc_diff.$$

[ $rc -eq 0 ] && echo "equivalent: $new"
exit $rc
