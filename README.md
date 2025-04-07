# Gene-Expression and Regulatory Analysis of Acute Myeloid Leukemia(AML) Using R and Microarray Data
This project analyzes microarray gene expression data in Acute Myeloid Leukemia (AML) using R. Differential expression analysis identifies genes altered in AML compared to healthy samples. Functional enrichment highlights key pathways and transcription factors involved in AML-related immune and inflammatory responses.
### Introduction

**Acute Myeloid Leukemia (AML)** is a type of cancer that affects the blood and bone marrow—the spongy tissue inside bones where blood cells are made. It is called "acute" because it progresses rapidly and requires immediate treatment. In AML, the bone marrow produces abnormal white blood cells, called myeloblasts, which multiply uncontrollably and interfere with the production of normal blood cells.

In this project, we aim to compare cancerous cells with healthy ones to identify genetic alterations that are unique to Acute Myeloid Leukemia (AML). By uncovering these specific changes, we hope to gain insights that could improve early diagnosis, guide treatment decisions, and potentially reveal new therapeutic targets.

### **1.Environment Setup and Package Installation**

The script begins by preparing the R environment for gene expression analysis using data from the NCBI GEO database. Necessary packages such as `GEOquery`, `limma`, `pheatmap`, and others are optionally installed to provide tools for data retrieval, manipulation, and visualization. The working directory is set to a specific folder where data files will be accessed or stored. Several important libraries are then loaded to enable functions for downloading GEO datasets, handling expression data, and creating plots. The dataset to be analyzed is identified by its GEO accession number (`GSE9476`) and its associated microarray platform (`GPL96`). Using the `getGEO()` function, the expression dataset is downloaded in matrix format with gene annotations included, and saved to a specified local directory. A conditional check is performed to select the correct dataset in case multiple platforms are present, ensuring that the appropriate data corresponding to the specified platform is extracted for further analysis.

### **2.Sample Annotation, Data Normalization Assessment, and Visualization**

Then, group labels were assigned to each sample based on their biological condition (e.g., CD34, BM, AML, PB), and the gene expression matrix was extracted from the GEO dataset. The dimensions and value ranges of the matrix were checked to determine whether the data was already normalized. A boxplot of the expression data was generated for quality control to visually assess normalization status. Although log2 transformation and quantile normalization were included as optional steps (commented out), they were not applied because the data appeared already normalized. Finally, a correlation heatmap was created using `pheatmap()` to visualize the relationships between samples based on their expression profiles, with group labels included to help interpret sample clustering.
