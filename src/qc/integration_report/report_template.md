---
title: "Integration QC report"
format: pdf
tbl-colwidths: [75,25]
---


| metric  | value  |
|--------|--------|
| adjusted rand index (ari)  | {ari_score}   |
| normalized mutual invormation (nmi)   | {nmi_score}   |
| average silhouette width (asw, label) | {asw_label}   |
| average silhouette width (asw, batch) | {asw_batch}   |
| principal component regression (pcr) | {pcr_score}   |
| graph connectivity | {graph_conn_score}   |
| average bio | {avg_bio}   |
: {{.striped}}

![]({umap_label}){{ffig-align="center" width=90%}}

![]({umap_batch}){{fig-align="center" width=90%}}