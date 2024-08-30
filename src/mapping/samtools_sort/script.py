import tempfile
import subprocess
from pathlib import Path
import yaml

## VIASH START
par = {
    'input': ['resources_test/cellranger_tiny_fastq/bam/possorted_genome_bam.bam'],
    'output_bam': 'test_output.bam',
    'output_bai': 'test_output.bam.bai'
}
meta = {
    'cpus': 2,
    'temp_dir': '/tmp',
    'config': 'src/mapping/htseq/config.vsh.yaml'
}
## VIASH END

def generate_args(par, config):
    # fetch arguments from config
    arguments = [
        arg
        for group in config["argument_groups"]
        for arg in group["arguments"]
    ]

    cmd_args = []

    for arg in arguments:
        arg_val = par.get(arg["name"].removeprefix("--"))
        # The info key is always present (changed in viash 0.7.4) 
        # in the parsed config (None if not specified in source config)
        info = arg["info"] or {}
        orig_arg = info.get("orig_arg")
        if arg_val and orig_arg:
            if not arg.get("multiple", False):
                arg_val = [arg_val]

            if arg["type"] in ["boolean_true", "boolean_false"]:
                # if argument is a boolean_true or boolean_false, simply add the flag
                arg_val = [orig_arg]
            elif orig_arg.startswith("-"):
                # if the orig arg flag is not a positional,
                # add the flag in front of each element and flatten
                arg_val = [str(x) for val in arg_val for x in [orig_arg, val]]

            cmd_args.extend(arg_val)

    return cmd_args

# read config arguments
config = yaml.safe_load(Path(meta["config"]).read_text())

print(">> Constructing command", flush=True)
cmd_args = [ "samtools", "sort" ] + generate_args(par, config)

# manually process cpus parameter
if 'cpus' in meta and meta['cpus']:
    cmd_args.extend(["--threads", str(meta["cpus"])])
# add memory
if 'memory_mb' in meta and meta['memory_mb']:
    import math
    mem_per_thread = math.ceil(meta['memory_mb'] * .8 / meta['cpus'])
    cmd_args.extend(["-m", f"{mem_per_thread}M"])

with tempfile.TemporaryDirectory(prefix="samtools-", dir=meta["temp_dir"]) as temp_dir:
    # add tempdir
    cmd_args.extend(["-T", str(temp_dir + "/")])

    # run command
    print(">> Running samtools sort with command:", flush=True)
    print("+ " + ' '.join([str(x) for x in cmd_args]), flush=True)
    subprocess.run(cmd_args, check=True)

if par.get("output_bai"):
    print(">> Running samtools index with command:", flush=True)
    cmd_index_args = ["samtools", "index", "-b", par["output_bam"], par["output_bai"]]
    print("+ " + ' '.join([str(x) for x in cmd_index_args]), flush=True)
    subprocess.run(cmd_index_args, check=True)