import subprocess

## VIASH START
par = {
    "input": [
        "resources_test/10x_5k_anticmv/fastqc/",
        "resources_test/10x_5k_anticmv/fastqc/",
    ],
    "output": "output",
}
## VIASH END

# Run MultiQC
subprocess.run(["multiqc", "-o", par["output"]] + par["input"])
