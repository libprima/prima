#!/bin/bash
# Usage: bash intelvars_unix.sh

source /opt/intel/oneapi/setvars.sh
ifort --version
where ifort
echo PATH="$PATH":"$(dirname "$(command -v ifort)")" >> "$GITHUB_ENV"
