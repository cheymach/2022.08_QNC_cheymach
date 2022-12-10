#import my data ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #organism specific
library(tximport) # package for getting Kallisto results into R
targets <- read_tsv("studydesign_11_14_2.txt")# read in your study design
path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
#filter and normalize data. Create violin plots showing the distribution 
library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = NPC:WTC11_2, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +# scale_fill_manual("blue","blue","yellow","yellow") +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = " ",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized") +
  scale_fill_manual(values=c("#F88379","#F88379","#0D98BA","#0D98BA")) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = NPC:WTC11_2, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = " ",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized") +
  scale_fill_manual(values=c("#F88379","#F88379","#0D98BA","#0D98BA")) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = NPC:WTC11_2, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = " ",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized") +
  scale_fill_manual(values=c("#F88379","#F88379","#0D98BA","#0D98BA")) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

#create PCA plot ----
library(DT)
library(gt)
library(plotly)

group <- targets$group
group <- factor(group)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=1) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

#create heat map ----
#NOTE: unlike my presentation, this version has a functioning heat map
#differential gene expression for heat map
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(plotly) 
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(differentiation = NPC - WTC11,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

#import libraries for generating the heatmap
library(limma) #we only use limma in this script for the 'avearrays' function
library(RColorBrewer) #need colors to make heatmaps
library(gplots)
library(heatmaply)
#decide color scheme for heat map
myheatcolors2 <- brewer.pal(name="RdYlBu", n=11)
#create a differential genes matrix to use for the heat map
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)

#hierarchical clustering 
#cluster the genes (rows) using the pearson method
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
#now cluster your samples (columns) using Spearman
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation

#define 2 clusters with k 
module.assign <- cutree(clustRows, k=2)

#assign a color to each module 
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
#
dev.new(width = 100, height = 50, unit = "cm")
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors2), scale='row', labRow=NA,
          density.info="none", trace="none",
          key = TRUE,
          cexRow=1, cexCol=1, margins=c(8,20)) 

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)




#volcano plot
library(EnhancedVolcano)
geneNames = rownames(myTopHits)
EnhancedVolcano(myTopHits.df, 
                lab = rownames(myTopHits),
                selectLab = c("RP11-777B9.5","GREB1L","LINC00371","MTRNR2L6"),
                x='logFC',
                y='adj.P.Val',
                legendLabels = c("neither FC nor p-value", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC)),
                title = "Differential gene expression (iPSCs vs NPCs)",
                pCutoff = .01,
                FCcutoff = 1.2)
#gene ontology analysis
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
#library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
#library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
#library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
#library(enrichplot) # great for making the standard GSEA enrichment plots
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = T, capped = T)

mygostplot <- gostplot(gost.res, interactive = F, capped = T)

publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

#try again with only list of genes with at least 1.2 fold change absolute value
#and with a p-value cutoff
#used excel to generate a table for this



myTopFoldChange <- myTopHits[1:4833,1:6]

myTopFCP <- read.csv(file = "myTopFoldChange_pvaluecutoff.csv")

myFCP <- myTopFoldChange[1:2171,1:6]
myTopFoldChange


gostFC.res <- gost(rownames(myTopFoldChange), organism = "hsapiens", correction_method = "fdr")
gostplot(gostFC.res, interactive = T, capped = T)


#GOs of interest: GO:0000902, GO:0015630, GO:0048468,HPA:0100000

myFCgostplot <- gostplot(gostFC.res, interactive = F, capped = T)

publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0000902", "GO:0015630", "GO:0048468", "HPA:0100000"),
  filename = NULL,
  width = NA,
  height = NA)

#create CSV of TopFoldChange for Gene Ontology analysis online, using geneontology.org
write.csv(myTopFoldChange, file = "myTopFoldChange.csv")
#write reference file for gene ontology analysis
write.csv(Tx, file = 'ReferenceGeneOntology.csv')
