import pytest
import sys
import gzip
import shutil
import loompy
import pysam

## VIASH START
meta = {
    'functionality_name': './target/executable/projection/velocyto/velocyto',
    'resources_dir': './resources_test/'
}
## VIASH END

# input data for bd bam
input_bam_bd = f"{meta['resources_dir']}/rna_velocity/velocyto/compatible_bd_input.bam"
input_gtf_bd = f"{meta['resources_dir']}/reference_gencodev41_chr1/reference.gtf.gz"
input_barcodes_bd = f"{meta['resources_dir']}/rna_velocity/velocyto/barcodes.txt"

# input data for 10x bam
input_bam_cellranger = f"{meta['resources_dir']}/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"
input_gtf_cellranger = f"{meta['resources_dir']}/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz"

def test_velocyto_cellranger(run_component, tmp_path):
    """Check whether component accepts compressed gtf files"""
    
    output_file = tmp_path / "foo" / "velocyto.loom"
    
    run_component([
            "--input", input_bam_cellranger,
            "--transcriptome", input_gtf_cellranger,
            "--output", str(output_file)])
    assert output_file.is_file()

    input_barcodes = set()
    with pysam.AlignmentFile(input_bam_cellranger, 'r') as input_bam:
        for read in input_bam:
            tags = dict(read.tags)
            cell_barcode = tags.get('CB')
            if cell_barcode:
                input_barcodes.add(cell_barcode.removesuffix("-1"))
    with loompy.connect(output_file) as ds:
        result_barcodes = {tag.removeprefix('velocyto:').removesuffix('x') for tag in ds.ca.CellID}
        assert result_barcodes.issubset(input_barcodes)
        assert ds.ca.keys() == ['CellID']
        assert ds.ra.keys(), ['Accession', 'Chromosome', 'End', 'Gene', 'Start' == 'Strand']
        rows, cols = ds.shape
        assert rows > 0
        assert cols > 0

def test_velocyto_bd_rhapsody(run_component, tmp_path):
    """Check whether component also accepts uncompressed gtf files"""

    output_file = tmp_path / "foo" / "velocyto.loom"
    transcriptome = tmp_path / "genes.gtf"
    
    with open(transcriptome, "wb") as gtf_uncompressed:
        with gzip.open(input_gtf_bd, 'rb') as gtf_compressed:
            shutil.copyfileobj(gtf_compressed, gtf_uncompressed)

    run_component([
        "--input", input_bam_bd,
        "--transcriptome", str(transcriptome),
        "--output", str(output_file),
        "--barcode", input_barcodes_bd
    ])
    assert output_file.is_file()

    input_barcodes = set()
    with open(input_barcodes_bd, 'r') as barcodes_file:
        for barcode in barcodes_file:
            input_barcodes.add(barcode.strip())
    
    with loompy.connect(output_file) as ds:
        result_barcodes = {tag.removeprefix('velocyto:').removesuffix('x') for tag in ds.ca.CellID}
        assert result_barcodes.issubset(input_barcodes)
        assert ds.ca.keys() == ['CellID']
        assert ds.ra.keys(), ['Accession', 'Chromosome', 'End', 'Gene', 'Start' == 'Strand']
        rows, cols = ds.shape
        assert rows > 0
        assert cols > 0

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))