import sys
import pytest
import gzip
import os
import re

## VIASH START
meta = {
    "name": "cellrnger_mkgtf",
    "resources_dir": "resources_test/",
    "executable": "target/docker/reference/cellranger_mkgtf/cellranger_mkgtf",
    "config": "src/reference/cellranger_mkgtf/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def subset_input_gtf(random_path):
    subset_input_path = random_path(extension="gtf.gz")
    with gzip.open(
        f"{meta['resources_dir']}/reference_gencodev41_chr1/reference.gtf.gz", "rt"
    ) as f_in:
        with gzip.open(subset_input_path, "wt") as f_out:
            for line in f_in:
                fields = line.split("\t")
                if len(fields) >= 4 and int(fields[3]) < 50001:
                    f_out.write(line)
    return subset_input_path


@pytest.mark.parametrize(
    "attributes", [["miRNA"], ["transcribed_unprocessed_pseudogene", "miRNA"]]
)
def test_gene_type_column(run_component, subset_input_gtf, random_path, attributes):
    output_gtf = random_path(extension="gtf.gz")
    args = ["--input_gtf", subset_input_gtf, "--output_gtf", output_gtf, "--attribute"]
    args.append(";".join([f"gene_type:{attribute}" for attribute in attributes]))

    print(args, flush=True)
    run_component(args)

    assert os.path.isfile(output_gtf), "Output GTF could not be found."

    with gzip.open(output_gtf, "rt") as f:
        unique_gene_types = {
            match
            for line in f
            for match in re.findall(r'gene_type "([^"]*)"', line.split("\t")[8])
        }
    assert (
        set(attributes) == unique_gene_types
    ), "Output GTF does not contain exactly the expected gene types."


def test_different_columns(run_component, subset_input_gtf, random_path):
    output_gtf = random_path(extension="gtf.gz")
    args = [
        "--input_gtf",
        subset_input_gtf,
        "--output_gtf",
        output_gtf,
        "--attribute",
        "gene_type:transcribed_unprocessed_pseudogene;transcript_id:ENST00000456328.2",
    ]

    run_component(args)
    assert os.path.isfile(output_gtf), "Output GTF could not be found."

    with gzip.open(output_gtf, "rt") as f:
        wrong_attributes_count = sum(
            1
            for line in f
            if dict(re.findall(r'(\S+) "([^"]*)"', line.split("\t")[8])).get(
                "gene_type"
            )
            != "transcribed_unprocessed_pseudogene"
            and dict(re.findall(r'(\S+) "([^"]*)"', line.split("\t")[8])).get(
                "transcript_id"
            )
            != "ENST00000456328.2"
        )
    assert (
        wrong_attributes_count == 0
    ), "Output GTF contains unexpected attribute values."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
