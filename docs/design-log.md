# Openpipeline - design


## Global overview

### Second general architecture (20211013)

![Global overview](figures/pipelines-target-p3.png)
_Overview single cell processing pipelines in openpipeline. Every rectangle is a pipeline by itself of which multiple versions can exist. Data aggregation is performed in the circles. The major parts of the pipeline consist of Ingestion, Per sample processing, Integration, Transformation and Reporting._ 

### First general architecture (20210929)

![Overview single cell processing pipelines - version 20210919](figures/pipelines-target-p1.png)
_Overview single cell processing pipelines._

Remarks:
1. What with multi-modal data integration? 
2. RNA-velocity is starting at the wrong location since it is only a different mapping step.
3. Muon objects to be used for the multi-modal data.


