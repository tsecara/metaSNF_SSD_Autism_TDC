# Beyond Diagnosis: Transdiagnostic Neuro-Behavioural Profiles of Social Cognition via Multimodal Data Integration Across Autism and Schizophrenia Spectrum Disorders 

## Project Overview

This repository contains code and data used in the analyses for the manuscript titled **"Beyond Diagnosis: Transdiagnostic Neuro-Behavioural Profiles of Social Cognition via Multimodal Data Integration Across Autism and Schizophrenia Spectrum Disorders."** It contains scripts for processing, visualizing, and analyzing multimodal data to identify biologically informed subgroups across individuals with SSD, autism, and typically developing controls (TDC) using Similiarity Network Fusion (SNF) and metaSNF.   

In order:
- **Data Prep** – Preparation of data for main analysis; neuroCombat to remove scanner effects. 
- **SNF/metaSNF Analysis** – Primary analysis involving main multimodal data integration steps to produce clusters.
- **Cluster Comparison** – Statistical comparison of input features and out-of-model features (for cluster validation) acorss clusters (and diagnostic reference groups).
- **Stability Testing** – Bootstrap resampling and visualization of clustering stability and agreement. 

This repository uses and builds off of code from previous studies:

1.  https://github.com/juliagallucci/ABCD_SNF/tree/main
2.  https://github.com/gracejacobs/SNF-NDD-Clustering
