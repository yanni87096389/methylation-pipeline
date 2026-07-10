# methylation-pipeline
This repository contains the analysis pipeline for identifying and prioritising candidate colorectal cancer-associated DNA methylation markers from public methylation datasets.

The overall aim of this project is to discover robust tumour-specific hypermethylated differentially methylated regions (DMRs), validate them across independent datasets, filter out markers with high background methylation in non-target normal tissues and whole blood, and prioritise a small panel of DMRs for downstream wet-lab validation.

---

## Project overview

The pipeline follows a stepwise marker discovery and prioritisation strategy:

1. Integrate three GEO EPIC/850K methylation datasets for DMR discovery.
2. Identify hypermethylated DMRs using DMRcate.
3. Validate candidate DMRs in five independent GEO 450K datasets.
4. Retain DMRs consistently validated in at least three independent datasets.
5. Filter retained DMRs against other normal tissues and whole blood methylation datasets.
6. Use machine learning models to assess the classification performance of the final DMR set.
7. Prioritise a small panel of DMRs for wet-lab validation.

---

## Workflow


Three GEO EPIC/850K discovery datasets
        ↓
CpG-level differential methylation analysis
        ↓
Meta-analysis across discovery datasets
        ↓
DMRcate DMR discovery
        ↓
409 hyper-DMRs
        ↓
Independent validation in five GEO 450K datasets
        ↓
234 DMRs validated in at least three datasets
        ↓
Background filtering in other normal tissues
        ↓
111 candidate DMRs
        ↓
Whole blood background filtering
        ↓
110 final candidate DMRs
        ↓
Machine learning training and testing
        ↓
Candidate DMR panel selection for wet-lab validation
