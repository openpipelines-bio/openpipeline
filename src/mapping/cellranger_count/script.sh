
## VIASH START
par_id="sample_id"
par_input="fastq_dir"
par_transcriptome="reference_dir"
par_libraries=""
par_feature_ref=""
par_log="log.txt"
par_output="output"
## VIASH END

if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

extra_params=( )

# process libraries parameter
if [ ! -z "$par_libraries" ]; then 
  extra_params+=( "--libraries=$par_libraries" )
fi

# process feature ref parameter
if [ ! -z "$par_feature_ref" ]; then 
  extra_params+=( "--feature-ref=$par_feature_ref" )
fi

cellranger count \
  --id="$par_id" \
  --fastqs="$par_input" \
  --transcriptome="$par_transcriptome" \
  "${extra_params[@]}" \
  --disable-ui \
  2>&1 | tee -a "$par_log"

if [ -d "$par_id/outs/" ]; then
  mv "$par_id"/outs/* "$par_output/"
fi