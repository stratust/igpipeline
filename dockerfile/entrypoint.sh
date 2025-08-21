#!/usr/bin/env bash
set -euo pipefail

: "${PERLBREW_ROOT:=/opt/perlbrew}"
export PATH="$PERLBREW_ROOT/bin:$PATH"

# Activate the perl without sourcing any rc file
eval "$("$PERLBREW_ROOT/bin/perlbrew" env perl-5.38.3)"

# Run user command or an interactive shell
if [[ $# -gt 0 ]]; then
  exec "$@"
else
  exec /bin/bash
fi
