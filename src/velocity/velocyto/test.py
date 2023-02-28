import unittest
import subprocess
from pathlib import Path
import gzip
import shutil
import tempfile
import loompy
import pysam


## VIASH START
meta = {
    'functionality_name': './target/native/projection/velocyto/velocyto',
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

class TestVelocyto(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_velocyto_cellranger(self):
        # check whether component accepts compressed gtf files
        self._run_and_check_output([
                "--input", input_bam_cellranger,
                "--transcriptome", input_gtf_cellranger,
                "--output", "./foo/velocyto.loom"])
        self.assertTrue(Path("./foo/velocyto.loom").is_file())
        input_barcodes = set()
        with pysam.AlignmentFile(input_bam_cellranger, 'r') as input_bam:
            for read in input_bam:
                tags = dict(read.tags)
                cell_barcode = tags.get('CB')
                if cell_barcode:
                    input_barcodes.add(cell_barcode.removesuffix("-1"))
        with loompy.connect("./foo/velocyto.loom") as ds:
            result_barcodes = {tag.removeprefix('velocyto:').removesuffix('x') for tag in ds.ca.CellID}
            self.assertTrue(result_barcodes.issubset(input_barcodes))
            self.assertListEqual(ds.ca.keys(), ['CellID'])
            self.assertListEqual(ds.ra.keys(), ['Accession', 'Chromosome', 'End', 'Gene', 'Start', 'Strand'])
            rows, cols = ds.shape
            self.assertTrue(rows > 0)
            self.assertTrue(cols > 0)


    def test_velocyto_bd_rhapsody(self):
        # check whether component also accepts uncompressed gtf files
        with tempfile.NamedTemporaryFile(suffix=".gtf", mode='wb') as genes_uncompressed:
            with gzip.open(input_gtf_bd, 'rb') as genes_compressed:
                shutil.copyfileobj(genes_compressed, genes_uncompressed)
                self._run_and_check_output([
                        "--input", input_bam_bd,
                        "--transcriptome", genes_uncompressed.name,
                        "--output", "./foo/velocyto.loom",
                        "--barcode", input_barcodes_bd]
                        )
        self.assertTrue(Path("./foo/velocyto.loom").is_file())
        input_barcodes = set()
        with open(input_barcodes_bd, 'r') as barcodes_file:
            for barcode in barcodes_file:
                input_barcodes.add(barcode.strip())
        with loompy.connect("./foo/velocyto.loom") as ds:
            result_barcodes = {tag.removeprefix('velocyto:').removesuffix('x') for tag in ds.ca.CellID}
            self.assertTrue(result_barcodes.issubset(input_barcodes))
            self.assertListEqual(ds.ca.keys(), ['CellID'])
            self.assertListEqual(ds.ra.keys(), ['Accession', 'Chromosome', 'End', 'Gene', 'Start', 'Strand'])
            rows, cols = ds.shape
            self.assertTrue(rows > 0)
            self.assertTrue(cols > 0)


if __name__ == "__main__":
    unittest.main()