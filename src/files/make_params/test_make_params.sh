#!/bin/bash

## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_name-XXXXXXXX")
# function clean_up {
#     rm -rf "$tmpdir"
# }
# trap clean_up EXIT

echo "> Running $meta_executable."
$meta_executable \
  --base_dir "./src" \
  --pattern "*.vsh.yaml" \
  --n_dirname_drop 1 \
  --n_basename_id 1 \
  --output "$tmpdir/output.yaml" \
  --path_name "path" \
  --group_name "param_list" \
  --id_name "id"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/output.yaml" ]] && echo "Output file could not be found!" && exit 1

if ! grep -qw param_list "$tmpdir/output.yaml"; then
    echo "Yaml key 'param_list" not found && exit 1
fi


# this component always be present because this test is executed.
if ! grep -qw "id: src_files_make_params" "$tmpdir/output.yaml"; then
    echo "Yaml key 'id: src_files_make_params" not found && exit 1
fi
echo "> Test succeeded!"