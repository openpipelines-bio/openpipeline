#!/bin/bash

set -eo pipefail

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

/usr/local/bin/cirro launch --port $VIASH_PAR_PORT --host 0.0.0.0 $VIASH_PAR_INPUT
