#!/bin/bash
# Source the perlbrew environment
source /opt/perlbrew/etc/bashrc

# Activate desired Perl version without spawning subshell
perlbrew use perl-5.38.3

# Run whatever the user provides (e.g., /bin/bash)
exec "/bin/bash"
