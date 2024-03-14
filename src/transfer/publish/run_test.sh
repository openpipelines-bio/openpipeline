#!/usr/bin/env bash
set -ex

touch test_file.txt

echo ">>> Testing if publish in current local dir works"
./publish \
  --input test_file.txt \
  --output another_file.txt

[[ ! -f another_file.txt ]] && echo "It seems no output file is generated" && exit 1

echo ">>> Testing if publish in new local dir works"
./publish \
  --input test_file.txt \
  --output adir/yadir/another_file.txt

[[ ! -d adir/yadir ]] && echo "It seems no output directory is generated" && exit 1
[[ ! -f adir/yadir/another_file.txt ]] && echo "It seems no output file is generated" && exit 1

echo ">>> Test finished successfully"
