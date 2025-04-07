# Gene-Expression Analysis of AML Using R and Microarray Data
This project analyses microarray gene expression data in Acute Myeloid Leukemia (AML) using R. Differential expression analysis identifies genes altered in AML compared to healthy samples. Functional enrichment highlights key pathways and transcription factors involved in AML-related immune and inflammatory responses.
### Introduction

**Acute Myeloid Leukemia (AML)** is a type of cancer that affects the blood and bone marrow—the spongy tissue inside bones where blood cells are made. It is called "acute" because it progresses rapidly and requires immediate treatment. In AML, the bone marrow produces abnormal white blood cells, called myeloblasts, which multiply uncontrollably and interfere with the production of normal blood cells.The first step of the project involved searching the NCBI GEO database to find a relevant dataset for comparing AML cells with healthy cells. The selected dataset, identified by its GEO accession number (GSE9476), provided gene expression data from microarray analysis, which was crucial for further analysis of genetic differences between AML and healthy samples.


---

### **1.Environment Setup and Package Installation**

The script begins by preparing the R environment for gene expression analysis using data from the NCBI GEO database. Necessary packages such as `GEOquery`, `limma`, `pheatmap`, and others are optionally installed to provide tools for data retrieval, manipulation, and visualization. The working directory is set to a specific folder where data files will be accessed or stored. Several important libraries are then loaded to enable functions for downloading GEO datasets, handling expression data, and creating plots. The dataset to be analyzed is identified by its GEO accession number (`GSE9476`) and its associated microarray platform (`GPL96`). Using the `getGEO()` function, the expression dataset is downloaded in matrix format with gene annotations included, and saved to a specified local directory. A conditional check is performed to select the correct dataset in case multiple platforms are present, ensuring that the appropriate data corresponding to the specified platform is extracted for further analysis.

---

### **2.Sample Annotation, Data Normalization Assessment, and Visualization**

Then, group labels were assigned to each sample based on their biological condition (e.g., CD34, BM, AML, PB), and the gene expression matrix was extracted from the GEO dataset.

PB (Peripheral Blood):
This sample comes from the peripheral blood, which is the blood circulating throughout the body. It is often used as a comparison to other samples because it's typically more readily available.

BM (Bone Marrow):
This sample is derived from the bone marrow, the tissue inside bones where blood cells are produced. Bone marrow samples provide insight into normal blood cell development and function.

CD34:
CD34 refers to a specific marker on stem cells, often used to isolate hematopoietic stem cells (HSCs), which are involved in blood cell production. This sample could be focused on hematopoietic stem cells or progenitor cells.

AML (Acute Myeloid Leukemia):
These samples are taken from patients with Acute Myeloid Leukemia, a type of cancer that affects the blood and bone marrow, leading to the production of abnormal white blood cells (myeloblasts). 

The dimensions and value ranges of the matrix were checked to determine whether the data was already normalized. A boxplot of the expression data was generated for quality control to visually assess normalization status. Although log2 transformation and quantile normalization were included as optional steps (commented out), they were not applied because the data appeared already normalized. Finally, a correlation heatmap was created using `pheatmap()` to visualize the relationships between samples based on their expression profiles, with group labels included to help interpret sample clustering.

---

### 3.Analysis of heatmap:

![Image](https://github.com/user-attachments/assets/b5c8a758-b225-40d4-a8a3-181091ebc17f)

From the correlation heatmap , we can understand how similar the gene expression profiles are between different samples based on their group labels (CD34, BM, AML, PB). Looking at the **hierarchical clustering**, we can see that the **Peripheral Blood (PB)** samples cluster closely together and are clearly separated from the other groups, suggesting they are biologically distinct and potentially not ideal for comparison with the other sample types. In contrast, **CD34** and **AML** samples are positioned closer together on the heatmap, indicating they share more similar gene expression patterns — making them more appropriate for comparative analysis. The **Bone Marrow (BM)** samples appear more distant from both AML and CD34. Interestingly, although AML samples form a broad cluster, the **correlation between AML samples themselves is not very tight**, likely reflecting the **high heterogeneity of cancer cells**, as individual AML cases can vary significantly in their molecular profiles. Overall, this analysis supports focusing on CD34 and AML samples for analysis.

---

### 4. Principal Component Analysis (PCA)

In the next step, Principal Component Analysis (PCA) was performed on the gene expression matrix to reduce dimensionality and identify patterns or clusters , where genes are represented as rows.

![Image](https://github.com/user-attachments/assets/cd148917-7e0e-4b60-9bb9-fd8198b9e37d)

In this PCA plot, a single point for each gene is plotted, and the analysis is on the variation of gene expression across all the samples. PC1 (x-axis) represents the highest contributor to variance in the data set, and PC2 (y-axis) represents the second-highest, which is orthogonal to PC1. A close horizontal spread shows that PC1 is populated with genes having overall high or low expression levels, regardless of whether they are varying between conditions. This can occur because always-expressed genes—e.g., housekeeping genes or genes of minimal biological interest—have still high numerical values, which dominate PCA if data isn't centered and scaled.

As a result, PC1 may reflect absolute expression magnitudes rather than meaningful variation across conditions, while PC2 may capture smaller, more subtle patterns. Without proper preprocessing (e.g., centering, scaling, or filtering low-variance genes), PCA can focus more on expression intensity than on biologically interesting differences in expression patterns between samples.

---

### 5.Centering each gene

```r
ex.scaled = t(scale(t(ex), scale = FALSE))
```

This line:

- Subtracts the **mean expression** of each gene across all samples.
- So each gene now has an average expression of **zero** across samples.
- This removes the **baseline level** of each gene, and keeps only the **variation pattern**.


### 💡 By centering each gene:

- We **remove the bias** of “high expression = more important.”
- PCA will now focus on **how genes behave differently across samples**, not just how big their numbers are.
- It’s especially useful when comparing samples or conditions, because now PCA reveals **relative expression patterns**, which are more biologically meaningful.


After removing the mean expression per gene, the total variance is more evenly distributed, and PC1 no longer soaks up the majority.  It means PCA is now more sensitive to **relative differences in expression across samples**, rather than just "how highly expressed" a gene was.

then we run another PCA on **mean-centered gene expression**.

![Image](https://github.com/user-attachments/assets/8ee20bfa-468e-48f5-9f3e-1a7fb2bb2280)

- The result shows that variance is **more evenly distributed**, not just captured in PC1.
- This improves clarity — PCA now focuses on **patterns of variation**, not raw expression levels.
  
---

### 6.Plotting PC of samples

Next, PCA was applied to the samples to uncover global expression patterns, highlight potential clustering or outlier samples, and simplify the data structure for easier analysis and visualization.

![Image](https://github.com/user-attachments/assets/9eca4a4e-6c12-4ccf-a730-21194295149c)

This PCA plot shows how samples cluster based on gene expression profiles, using only the first two principal components (PC1 and PC2).

### Key Observations

1. **PB Samples (Purple): Distinct and Tightly Clustered**
    
    Peripheral blood (PB) samples form a compact cluster that is clearly separated from all other groups. This indicates that their gene expression patterns are **consistently distinct**, suggesting PB samples are **biologically different** and may not be suitable for direct comparison with the other sample types.
    
2. **BM Samples (Green): Consistent Expression**
    
    Bone marrow (BM) samples group closely together near the top center of the plot. This tight clustering reflects **uniform gene expression** within this group. While there is slight overlap with AML, BM samples are largely distinct, indicating stable expression patterns typical of normal bone marrow.
    
3. **CD34 Samples (Cyan): Separated but Relevant**
    
    CD34 samples are positioned on the far right, clearly separated from PB and BM. Their proximity to AML samples, however, suggests some **biological similarity or transitional features**, making them a potentially **informative group for comparison** in understanding AML development.
    
4. **AML Samples (Red): Scattered and Heterogeneous**
    
    AML samples are more widely dispersed across the plot, mostly between BM and CD34. This **scattered distribution** reflects the known **heterogeneity of AML**, as gene expression in cancer can vary greatly between individuals. The lack of tight clustering highlights the complexity and diversity within AML cases.
    


This plot confirms that:

- There are **distinct gene expression signatures** between sample groups.
- **CD34 and AML** are **relatively close**, making them ideal for differential expression comparison.
- **PB samples** are very different → they may not be biologically comparable to AML or CD34 (maybe worth excluding in some analyses).
- **AML’s diversity** is biologically expected due to its variable nature across patients or subtypes.

---

  ### **7.Identifying Differentially Expressed Genes Between AML and CD34 Samples**

A differential expression analysis was performed using the `limma` package to compare gene expression between **AML** and **CD34** samples.

 First, the group labels (`gr`) are converted into a factor and added to the `gset` object. A design matrix is then created without an intercept to model each group separately. A linear model is fitted to the expression data using `lmFit()`, estimating average gene expression per group. A contrast is defined to compare AML against CD34 (`AML - CD34`), and the model is refitted accordingly. 

Next, empirical Bayes moderation is applied with `eBayes()` to improve statistical reliability. The `topTable()` function extracts the list of genes from the analysis results, showing them in order of how statistically significant their expression differences are , with adjusted p-values (FDR) and log fold changes. Finally, a simplified table containing only gene symbols, IDs, adjusted p-values, and log fold changes is saved to a file. This process identifies genes that are significantly up- or down-regulated in AML compared to CD34.

---

### **8.Finding Upregulated and Downregulated Genes in AML**

Then the results of the differential expression analysis were filtered to identify genes that are significantly upregulated or downregulated in **AML** compared to **CD34** samples. First, it selects genes with a **log fold change greater than 1** and **adjusted p-value below 0.05**, indicating that these genes are **significantly overexpressed in AML**; the number of such genes is checked using `dim()`. To ensure uniqueness, gene symbols are extracted and duplicates are removed using `unique()`. Since some entries may include multiple gene names separated by "///", the `strsplit2()` function is used to split them into individual genes. The final list of upregulated genes is saved to a text file. 

The same steps are repeated for **downregulated genes**, defined by a **log fold change less than -1** and an adjusted p-value below 0.05, indicating genes that are **less expressed in AML (i.e., higher in CD34)**. This results in two clean gene lists: one for genes upregulated in AML and another for those downregulated, both ready for downstream analysis such as pathway enrichment or biomarker exploration.

---

### 9.Gene antology pathway analysis

After generating the gene lists (e.g., `aml.up.genes`), we can perform functional enrichment analysis using the **Enrichr** website. To do this, we simply copy the gene symbols from the list (e.g., upregulated genes in AML) and paste them into the input box on the Enrichr platform. Enrichr compares our gene list against a wide range of biological databases, including **TRANSFAC** and **JASPAR**, which contain information about transcription factors (TFs) and their known gene targets. When we select these databases, Enrichr shows us a list of transcription factors that are likely to regulate the genes we submitted.

 By hovering the mouse over a transcription factor in the results, we can see exactly which genes in our list are predicted to be its targets. For example, if the transcription factor **NFE2** appears near the top of the results, it suggests that NFE2 may play a key regulatory role in controlling the AML-associated genes we identified. This provides biological insight into potential upstream regulators driving gene expression changes in AML.

By clicking on a specific pathway in the Enrichr results, we can explore the **biological pathways associated with the upregulated genes**, such as those involved in the **immune system**. This helps us understand the broader functional context in which these genes operate.

Additionally,reviewing **Gene Ontology (GO) terms** reveals that many of the upregulated genes are linked to processes like the **inflammatory response**, providing further evidence that **inflammation plays a significant role in AML** pathogenesis. These insights offer a clearer view of the molecular mechanisms driving the disease.

---

### **Conclusion of Gene Expression and Functional Enrichment Analysis in AML**

The analysis revealed distinct gene expression patterns across sample types, with PCA showing that **AML samples are heterogeneous**, clustering between **CD34** and **BM** samples, while **PB** samples were clearly distinct. Differential expression analysis between AML and CD34 identified a set of significantly up- and down-regulated genes in AML. 

Functional enrichment of the upregulated genes using Enrichr highlighted key **biological pathways**, notably those related to the **immune system**, and **transcription factors** such as **NFE2** that may regulate these genes. Gene Ontology analysis further indicated a strong association with the **inflammatory response**, suggesting that immune and inflammatory processes are actively involved in the molecular landscape of AML. These results provide valuable insights into the regulatory mechanisms and biological functions underlying gene expression changes in AML.

## ✉ Contact
If you have any questions or feedback, feel free to reach out!
- **LinkedIn:** www.linkedin.com/in/elaheh-p-9918432a6

 **If you find this project useful, please ⭐ star the repository!** ⭐
