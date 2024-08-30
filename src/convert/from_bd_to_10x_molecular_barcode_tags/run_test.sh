#!/bin/bash



echo "> Testing missing input file can be detected."
"$meta_executable" \
  --input "oops.BAM" \
  --output "foo.bam" \
  -t 8 \
  --bam > /dev/null

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo "Test failed: non-existent input file was not detected." && exit 1
fi

echo "> Testing if BAM output can be created"
"$meta_executable" \
  --input "$meta_resources_dir/output_raw/Combined_sample_Bioproduct.bam" \
  --output "foo.bam" \
  -t 8 \
  --bam

[ ! -f "foo.bam" ] && { echo Output file could not be found; exit 1; }
readarray -t output_tags < <( samtools view --no-header "foo.bam" | grep -oP "(?<=UB:Z:).*[\s]" )
readarray -t input_tags < <( samtools view --no-header "$meta_resources_dir/output_raw/Combined_sample_Bioproduct.bam" | grep -oP "(?<=MA:Z:).*[\s]" )
[ "${output_tags[*]}" == "${input_tags[*]}" ] || { echo "Input tags differ from output tags!"; exit 1; }

echo "> Testing if SAM output can be created"
"$meta_executable" \
  --input "$meta_resources_dir/output_raw/Combined_sample_Bioproduct.bam" \
  --output "foo.sam" \
  -t 8
[ ! -f "foo.sam" ] && { echo Output file could not be found; exit 1; }

readarray -t output_tags_sam < <( samtools view --no-header "foo.sam" | grep -oP "(?<=UB:Z:).*[\s]" )
[ "${output_tags_sam[*]}" == "${input_tags[*]}" ] || { echo "Input tags differ from output tags!"; exit 1; }


echo "> Test succeeded!"