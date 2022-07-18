# Demux Pipeline

Demultiplication using one of the tools supported

## Input/output



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to sequencer output folder <details><summary>Help</summary><small>The path to the output folder of the sequencing run.  </small></details>| `string` |  |  |  |
| `sample_sheet` | The sample sheet for the samples to be demultiplexed. <details><summary>Help</summary><small>The sample sheet to be used for the demultiplexing. This file contains the Lane, Sample and Index column. To select the lane use the corresponding lane value or the `*` for using all lanes. For Sample the id for the output fastqs has to be provided. For the Index the used 10x reference barcodes have to be supplied. For more information see [10x mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq).

```
Lane,Sample,Index
1,test_sample,SI-P03-C9
```

*Simple sample sheet*

```
[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,s1,test_sample,,,SI-TT-D9,SI-TT-D9,SI-TT-D9,SI-TT-D9,p1,
```
*Dual indexing sample sheet.*</small></details>| `string` |  |  |  |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `demultiplexer` |  | `string` | None |  |  |
| `publishDir` | Path to the output folder which will contain the cell ranger output. <details><summary>Help</summary><small>The S3 output folder to write the data to which will contain the result data.</small></details>| `string` |  |  |  |

