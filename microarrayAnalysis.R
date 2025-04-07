
#### Load useful packages
install.packages(c("GEOquery" , 'limma', 'pheatmap','ggplot2' , 'gplots' , 'reshape2' , 'plyr'))

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

#### Set Working Directory
setwd('D:/bioinf/dr.sharifi/Dr sharifi/Bioinf pishrafte 1/Sharifworkshop/ ')

#### specify dataset info
series = 'GSE9476'    # accession number of the data 
platform = 'GPL96'    #  microarray platform used (type of chip)

#### Download the Data
gset = getGEO( series , GSEMatrix= TRUE , AnnotGPL= TRUE , destdir= 'Data/')
length(gset)

- Downloads the gene expression data from GEO.
- GSEMatrix=TRUE` gets it in matrix format (easier to work with).
- AnnotGPL=TRUE` adds gene annotation (gene symbols, etc.).
- destdir='Data/'` saves the data locally in a folder called `Data/`.

#### Select the Correct Platform (if multiple)
if (length(gset) > 1) idx = grep(platform, attr(gset, "names")) else idx = 1
gset = gset[[idx]]

#### Assign group labels to samples
gr = c('CD34' , rep('BM' , 10) , rep('CD34', 7) , rep('AML', 26) , rep('PB', 10), rep('CD34', 10))


#### Extract the expression matrix
ex = exprs(gset)  
dim(ex)
max(ex)       #to see if it's normalized or not. if not we use log2
min(ex)

#### log2 scale, if required
#ex = log2( ex+ 1)
# exprs(gset) = ex

##### Boxplot for quality control
pdf("results/boxplot.pdf" , width= 64)
boxplot(ex)
dev.off()
#we saw that it is normalized


####Normalize , if required
# ex = normalizeQuantiles(ex)
# exprs(gset) = ex


#### Create a correlation heatmap
pdf('results/CorHeatmap.pdf' , width = 15 , height= 15 )
pheatmap(cor(ex), labels_row = gr , labels_col = gr , color= bluered(256) , border_color = NA)  # it calculates correlation between every 2 samples in the ex matrix
dev.off()

####  Initial PCA on raw expression data
pc = prcomp(ex)
pdf("results/PCA.pdf")
plot(pc)
plot(pc$x[ , 1:2])    #only pc1 and pc2
dev.off()
# shows pca1 is more important

names(pc)
dim(pc$x) # showing the genes
dim(pc$rotation) #showing the samples
dim(pc)

####Center each gene (row) manually
#we substract each row from the mean expression of it
ex.scaled = t(scale(t(ex) , scale= FALSE))  # scale makes the mean of each column equal to zero
 # since the genes are on row we should use t(ex) so they will go the columns and
# then we do transpose againe to bring them back to rows after scaling


####PCA after centering the genes
pc = prcomp(ex.scaled)
pdf('results/PC_scaled.pdf')
plot(pc)
plot(pc$x[ , 1:2])
dev.off()

#### then we plot PC of samples
pcr = data.frame(pc$rotation[ , 1:2] , Group = gr)
pdf('results/PCA_samples.pdf')
ggplot( pcr , aes( x= PC1 , y= PC2 , color= Group )) + geom_point(size = 3) + theme_bw()
dev.off()

#### Differential Expression Analysis
gr = factor(gr)
gset$description = gr
design = model.matrix(~description + 0, gset)
colnames(design) = levels(gr)
head(design)

fit = lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cont.matrix = makeContrasts(AML - CD34, levels=design)
fit2 = contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
head(tT)
colnames(tT)  # shows all the column names that are obtained from annotation

tT = subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
head(tT)
write.table(tT, 'results/ AML_CD34.txt', row.names=F, sep="\t", quote=F)

aml.up = subset(tT , logFC > 1 & adj.P.Val < 0.05)
# To find the genes that have a high expression in cancerous cells
dim(aml.up)
aml.up.genes = unique(aml.up$Gene.symbol) # to delete the repetitive ones
aml.up.genes =as.character( strsplit2(aml.up.genes ,'///'))
write.table(aml.up.genes , file= 'results/AML_CD34_UP.txt' , quote= F , row.names = F , col.names = F)


aml.down = subset(tT , logFC < -1 & adj.P.Val < 0.05)
# To find the genes that have a low expression in cancerous cells
dim(aml.down)
aml.down.genes = unique(aml.down$Gene.symbol) # to delete the repetitive ones
aml.down.genes = as.character(strsplit2( aml.down.genes ,'///'))
write.table(aml.down.genes , file= 'results/AML_CD34_DOWN.txt' , quote= F , row.names = F , col.names = F)

#### Gene antology pathway analysis

