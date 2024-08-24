#!/usr/bin/env bash
# Note: This test doesn't work yet

set -ex

./run_cirrocumulus --port 5005
wget -O page.html echo http://localhost:5005
[[ ! -f page.html ]] && echo "Page not found! Webserver may not be running!" && exit 1

echo ">>> Test finished successfully"
