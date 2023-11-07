# Kasikova et al 2023

This repository host the code used in our paper titled "Tertiary lymphoid structures and B cells determine clinically relevant T cell phenotypes in ovarian cancer." currently in revision.

## Abstract

Intratumoral tertiary lymphoid structures (TLSs) have been associated with improved outcome in a variety of cancer patient cohorts, reflecting the key role of TLSs in the initiation of tumor-targeting immunity. Here, we demonstrate that high-grade serous ovarian carcinoma (HGSOC) contains distinct immune aggregates with varying degrees of organization and maturation. Specifically, mature TLSs (mTLS) as forming only in 16% of HGSOCs with relatively elevated tumor mutational burden (TMB) were associated with an increased intratumoral density of CD8+ effector T (TEFF) cells  and TIM3+PD1+, hence poorly immune checkpoint inhibitor (ICI)-sensitive, CD8+ T cells. Converesely, CD8+ T cells from immunologically hot tumors like non-small cell lung carcinoma (NSCLC) were enriched in ICI-responsive TCF1+ PD1+ T cells. Spatial B-cell profiling identified patterns of in situ maturation and differentiation that were associated with mTLSs. Moreover, B-cell depletion promoted signs of a dysfunctional CD8+ T cell compartment among tumor-infiltrating lymphocytes from freshly isolated HGSOC and NSCLC biopsies. Taken together, our data demonstrate that – at odds with NSCLC – HGSOC is associated with a low density of follicular helper T cells and thus develops a limited number of mTLS that might be insufficient to preserve a ICI-sensitive TCF1+PD1+ CD8+ T cell phenotype. These findings point to key quantitative and qualitative differences between mTLSs in ICI-responsive vs ICI-irresponsive neoplasms that may guide the development of alternative immunotherapies for patients with HGSOC.

Experimental design
![Kasikova_expDesign](https://github.com/pepap/sotio-NatComm-2023/assets/7227977/832feb74-d803-4fe6-9d19-5aa4c9fdd396)

### Pre-processing of the Spatial transcriptomic data

Script detailing the pre-processing of Visium data is “Spatial_transcriptomics_Seurat.R” and take as input the space ranger output (10X genomics).

### R data analysis

The R scripts used to generate the figures can be found in the individual folders.

Other scripts are ordered and labeled by the analysis.

Data used in our analysis are also provided :
+ annotations and target regions
  - agi.ensembl.gr.rda                                                          ...  target regions from AGILENT sequencing experiment, ENSEMBL format
  - agi.mm10.gr.rda                                                             ...  target regions from AGILENT sequencing experiment, UCSC format
  - gencode.vM24.chr_patch_hapl_scaff.annotation.corrSeqLevels.trx2gene.PC.rda  ...  GENCODE gene annotation
  - targetReg.Twist_Mouse_Exome_Target_Rev1_7APR20.mm10.gr.rda                  ...  annotated target regions from ILLUMINA sequencing experiment
  - targetReg.agi.mm10.gr.rda                                                   ...  annotated target regions from AGILENT sequencing experiment

+ data for the figure 5E : plot--figure5E-data.tar.gz
  - Extract data :
```
tar xzfv plot--figure5E-data.tar.gz
```
