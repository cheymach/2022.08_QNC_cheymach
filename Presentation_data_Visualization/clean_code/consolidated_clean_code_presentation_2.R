#import my data ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #organism specific
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
