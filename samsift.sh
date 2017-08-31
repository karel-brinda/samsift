#! /usr/bin/env bash

set -o pipefail

readonly PROGDIR=$(dirname $0)
readonly ARGS="$@"

${PROGDIR}/samsift/samsift.py $ARGS
