### RNA seq analysis from N2 and skpo worm strains, collaboration with Danielle

# In this script we will analyse the RNA seq with the following conditions:
# 	- N2 with OP50
# 	- N2 with OGRF
# 	- skpo with OP50 
# 	- ep2 with OGRF

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# maybe use ulimit -s 16384 before start R console

# analysis directory:  "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/RNA_seq"


library(tximport)
library(DESeq2)


library(tidyverse)


here::set_here()
library(here)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(openxlsx)
library(viridis)


# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name


# prepare a list with file names
files = file.path(dir,"quants", samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist

# create an object with all the reference genes and transcripts
# THIS DOES NOT SEEM TO BE WORKING
# txdb = TxDb.Celegans.UCSC.ce11.refGene::TxDb.Celegans.UCSC.ce11.refGene
# k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
# tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# let's make our database from ensembldb 
ah = AnnotationHub::AnnotationHub(localHub = FALSE, proxy='127.0.0.1:10801')
ahDb = AnnotationHub::query(ah, pattern = c("Caenorhabditis elegans", "EnsDb", 98))
ahEdb = ahDb[[1]]
# generate the database 

# generate the database 
tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>% 
  tbl_df() %>%
  dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
  unnest

# join transcription info with gene ids and entrezids
info.join = tx2gene.complete %>% 
  tbl_df() %>%
  dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
  left_join(info)

write.csv(info, here('summary','gene_ids_mapping.csv'))

# subset to have tx_id in first column, and gene_id in second
tx2gene = tx2gene.complete[,c(1,7)]

# import quantification data 
txi = tximport(files, type = "salmon", tx2gene = tx2gene)



### starting analysis with DESeq2
# create DESeq data type to be analysed
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, design = ~  Sample)

# # prefilter, but that might not be necessary
# keep = rowSums(counts(ddsTxi)) >= 10
# ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "OP50")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)



### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% tbl_df()

gene_counts = gene_counts %>% 
  gather(Name, counts, OP50_1: MG1655_PDXJK_4) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(tbl_df(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 


### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
norm_gene_counts = counts(dds, normalized = TRUE)
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
g_counts = counts(dds, normalized = FALSE)

write.table(g_counts, here('summary','gene_counts.txt'), quote = FALSE, sep = '\t') 
write.table(norm_gene_counts, here('summary','norm_gene_counts.txt'), quote = FALSE, sep = '\t') 




# get results and tidy it
res = results(dds) 



# results with different shape of contrasts, tidy
# contrasts are written in this script that they show regulation in pathogenic species rather than in OP50 or N2
# Danielle suggested that
res.OP50= results(dds,   contrast = c("Sample", "OP50" , "OP50_DM"))  
res.OP50 = lfcShrink(dds, contrast = c("Sample", "OP50" , "OP50_DM"), res = res.OP50, type = 'ashr')

res.MG = results(dds,  contrast = c("Sample",  "MG1655", "MG1655_DM")) 
res.MG = lfcShrink(dds, contrast = c("Sample",  "MG1655", "MG1655_DM"), res = res.MG, type = 'ashr')

res.OP50.MG = results(dds,  contrast = c("Sample",  "OP50", "MG1655"))   
res.OP50.MG = lfcShrink(dds, contrast = c("Sample", "OP50", "MG1655"), res = res.OP50.MG, type = 'ashr')

res.OP50.MG.dm = results(dds, contrast = c("Sample",   "OP50_DM", "MG1655_DM")) 
res.OP50.MG.dm = lfcShrink(dds, contrast = c("Sample", "OP50_DM", "MG1655_DM"), res = res.OP50.MG.dm, type = 'ashr')



# tidying the results
res.OP50.tidy = as_tibble(res.OP50, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50',
  Contrast_description = 'Comparison of OP50 vs OP50_DM') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.MG.tidy = as_tibble(res.MG, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'MG',
  Contrast_description = 'Comparison of MG vs MG_DM') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OP50.MG.tidy = as_tibble(res.OP50.MG, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50_MG',
  Contrast_description = 'Comparison of OP50 vs MG') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OP50.MG.dm.tidy = as_tibble(res.OP50.MG.dm, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50_MG_DM',
  Contrast_description = 'Comparison of OP50_DM vs MG_DM') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

results.complete = res.N2.tidy %>% rbind(res.skpo.tidy, res.OGRF.tidy, res.OP50.tidy)

# write results in excel files
list_of_datasets = list('Comparison of OP50 vs OP50_DM' = res.OP50.tidy, 
                        'Comparison of MG vs MG_DM' = res.MG.tidy, 
                        'Comparison of OP50 vs MG' = res.OP50.MG.tidy,
                        'Comparison of OP50_DM vs MG_DM' = res.OP50.MG.dm.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'), colNames = T, rowNames = F) 



# transofrm data

vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)


# plot differences between different transformation data methods
df = bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("mean", "sd")  

ggplot(df, aes(x = mean, y = sd)) + 
  geom_hex(bins = 100) +
  # coord_fixed() + 
  facet_grid( . ~ transformation) +
  theme_light()


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'transformation_comparison.pdf'),
             height = 4.5, width = 12, useDingbats = FALSE)


###
# Sample distances

sampleDists = dist(t(assay(rld)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
names = colnames(sampleDists) %>%
  str_split('_', simplify = T) %>%
  data.frame %>% tbl_df() %>%
  unite(sample, X1, X2, sep = " - ") %>%
  dplyr::select(sample) %>%
  t %>% as.vector

colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
        col = colors)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Euclidean_distances_samples.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


##########
# PCA data

# plotPCA(rld, intgroup = c("Bacteria", "Worm"))

pcaData = plotPCA(rld, intgroup = c("Bacteria"), returnData = TRUE)
pcaData


# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}

# get info for the ellipses
ell = pcaData %>% group_by(Bacteria) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Bacteria)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, linetype = Bacteria), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, 
                               linetype = Bacteria, fill = Bacteria), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme_classic()


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)



## ploting some data
###
# quick check on a couple genes of interest
# WBGene00045401 = eol-1
# WBGene00018997 = F57B9.3
gns = c('WBGene00045401', 'WBGene00018997')

# acdh-1 and acdh-2
# WBGene00016943 = acdh-1
# WBGene00015894 = acdh-2
gns = c('WBGene00016943', 'WBGene00015894')

# WBGene00009221 = acs-2
# WBGene00010759 = cysl-2 

gns = c('WBGene00009221', 'WBGene00010759')


# WBGene00020343 = acs-3
# WBGene00003623 = nhr-25
gns = c('WBGene00016943','WBGene00009221', 'WBGene00020343', 'WBGene00003623')



# for Filipe, for QQ review
# hlh-30
gns = c('WBGene00020930')

# fmo-2
gns = c('WBGene00001477')

# ftn-1
gns = c('WBGene00001500')

# hif-1
gns = c('WBGene00001851')

# sqst-1
gns = c('WBGene00011737')

# fat-7
gns = c('WBGene00001399')

# fat-1
gns = c('WBGene00001393')

# cbp-1 (PE300 ortholog)
gns = c('WBGene00000366')

# nol-6 (nucleolar protein 6)
gns = c('WBGene00000608')

gene_counts %>%
  dplyr::filter(gene_id %in% gns) %>%
  ggplot(aes(y = counts, x = Sample)) +
  geom_boxplot(aes(fill = Bacteria)) +
  # scale_y_continuous(trans='log2') +
  # facet_wrap(~symbol) +
  facet_wrap(~gene_name, scales = 'free_y') +
  labs(x = 'Sample',
       y = 'Normalised counts') +
  theme_classic()

dev.copy2pdf(device = cairo_pdf,
             file =  here('summary', 'gene_counts.pdf'),
             height = 7, width = 6, useDingbats = FALSE)







#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################











#####################
### Batch effect ####
#####################


samples.batch = samples %>% mutate(Batch = as.factor(Batch))

dds.batch = DESeqDataSetFromMatrix(countData = g_counts,
                              colData = samples.batch,
                              design = ~Batch + Bacteria)
dds.batch = DESeq(dds.batch)

resultsNames(dds.batch)

res = results(dds.batch, name="Bacteria_MG1655_PDXJK_vs_MG1655")
summary(res)







##########
# PCA data

# plotPCA(rld, intgroup = c("Bacteria", "Worm"))



rld.batch = rlog(dds.batch, blind = FALSE)

assay(rld.batch) <- limma::removeBatchEffect(assay(rld.batch), rld.batch$Batch)

pcaData = plotPCA(rld.batch, intgroup = c("Bacteria"), returnData = TRUE)
pcaData


# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}

# get info for the ellipses
ell = pcaData %>% group_by(Bacteria) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Bacteria)) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, linetype = Bacteria), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, 
                               linetype = Bacteria, fill = Bacteria), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme_classic()


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld_batch.pdf'),
             height = 8, width = 9, useDingbats = FALSE)




###
# Sample distances

sampleDists = dist(t(assay(rld.batch)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
names = colnames(sampleDists) %>%
  str_split('_', simplify = T) %>%
  data.frame %>% tbl_df() %>%
  unite(sample, X1, X2, sep = " - ") %>%
  dplyr::select(sample) %>%
  t %>% as.vector

colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
        col = colors)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Euclidean_distances_samples_batch.pdf'),
             height = 8, width = 9, useDingbats = FALSE)




### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds.batch, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% tbl_df()

gene_counts = gene_counts %>% 
  gather(Name, counts, OP50_1: MG1655_PDXJK_4) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(tbl_df(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 


### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
norm_gene_counts = counts(dds, normalized = TRUE)
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
g_counts = counts(dds, normalized = FALSE)

write.table(g_counts, here('summary','gene_counts.txt'), quote = FALSE, sep = '\t') 
write.table(norm_gene_counts, here('summary','norm_gene_counts.txt'), quote = FALSE, sep = '\t') 




# get results and tidy it
res = results(dds.batch) 



# results with different shape of contrasts, tidy
# contrasts are written in this script that they show regulation in pathogenic species rather than in OP50 or N2
# Danielle suggested that
res.OP50 = results(dds,   contrast = c("Sample", "OP50_DM" , "OP50"))  
res.OP50 = lfcShrink(dds, contrast = c("Sample", "OP50_DM" , "OP50"), res = res.OP50, type = 'ashr')

res.MG = results(dds,  contrast = c("Sample",  "MG1655_DM", "MG1655")) 
res.MG = lfcShrink(dds, contrast = c("Sample",  "MG1655_DM", "MG1655"), res = res.MG, type = 'ashr')

res.OP50.MG = results(dds,  contrast = c("Sample",  "MG1655", "OP50"))   
res.OP50.MG = lfcShrink(dds, contrast = c("Sample", "MG1655", "OP50"), res = res.OP50.MG, type = 'ashr')

res.OP50.MG.dm = results(dds, contrast = c("Sample",   "MG1655_DM", "OP50_DM")) 
res.OP50.MG.dm = lfcShrink(dds, contrast = c("Sample",   "MG1655_DM", "OP50_DM"), res = res.OP50.MG.dm, type = 'ashr')




# tidying the results
res.OP50.tidy = as_tibble(res.OP50, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50',
  Contrast_description = 'Comparison of OP50_DM vs OP50') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.MG.tidy = as_tibble(res.MG, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'MG',
  Contrast_description = 'Comparison of MG_DM vs MG') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OP50.MG.tidy = as_tibble(res.OP50.MG, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50_MG',
  Contrast_description = 'Comparison of MG vs OP50') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OP50.MG.dm.tidy = as_tibble(res.OP50.MG.dm, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50_MG_DM',
  Contrast_description = 'Comparison of MG_DM vs OP50_DM') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

results.complete = res.OP50.tidy %>% rbind(res.MG.tidy, res.OP50.MG.tidy, res.OP50.MG.dm.tidy)

# write results in excel files
list_of_datasets = list('Comparison of OP50_DM vs OP50' = res.OP50.tidy, 
                        'Comparison of MG_DM vs MG' = res.MG.tidy, 
                        'Comparison of MG vs OP50' = res.OP50.MG.tidy,
                        'Comparison of MG_DM vs OP50_DM' = res.OP50.MG.dm.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats_batch.xlsx'), colNames = T, rowNames = F) 


results.complete %>%
  filter(gene_name == 'fat-7')





## ploting some data
###
# quick check on a couple genes of interest
# WBGene00045401 = eol-1
# WBGene00018997 = F57B9.3
gns = c('WBGene00045401', 'WBGene00018997')

# acdh-1 and acdh-2
# WBGene00016943 = acdh-1
# WBGene00015894 = acdh-2
gns = c('WBGene00016943', 'WBGene00015894')

# WBGene00009221 = acs-2
# WBGene00010759 = cysl-2 

gns = c('WBGene00009221', 'WBGene00010759')


# WBGene00020343 = acs-3
# WBGene00003623 = nhr-25
gns = c('WBGene00016943','WBGene00009221', 'WBGene00020343', 'WBGene00003623')



# for Filipe, for QQ review
# hlh-30
gns = c('WBGene00020930')

# fmo-2
gns = c('WBGene00001477')

# ftn-1
gns = c('WBGene00001500')

# hif-1
gns = c('WBGene00001851')

# sqst-1
gns = c('WBGene00011737')

# fat-7
gns = c('WBGene00001399')

# fat-1
gns = c('WBGene00001393')

# cbp-1 (PE300 ortholog)
gns = c('WBGene00000366')

# nol-6 (nucleolar protein 6)
gns = c('fat-7')

gene_counts %>%
  dplyr::filter(gene_name %in% gns) %>%
  ggplot(aes(y = counts, x = Sample)) +
  geom_boxplot(aes(fill = Bacteria)) +
  # scale_y_continuous(trans='log2') +
  # facet_wrap(~symbol) +
  facet_wrap(~gene_name, scales = 'free_y') +
  labs(x = 'Sample',
       y = 'Normalised counts') +
  theme_classic()




###########################
### Enrichment analysis ###
###########################


### this is a test to see the differentially enriched genes in https://amp.pharm.mssm.edu/WormEnrichr/
### check this website as well: http://www.wormcat.com/index


# positive threshold
pos = 0
# negative threshold
neg = 0

OP50.gns.up = res.OP50.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange >= pos) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(OP50.gns.up)), 'OP50_UP_genes.txt', quote = FALSE, col.names = F, row.names = F)
OP50.gns.down = res.OP50.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange <= neg) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(OP50.gns.down)), 'OP50_DOWN_genes.txt', quote = FALSE, col.names = F, row.names = F)

# if log2FC is > 0, this means that the gene is DOWN-regulated in GCB
res.MG.up = res.MG.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange >= pos) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(res.MG.up)), 'MG_UP_genes.txt', quote = FALSE, col.names = F, row.names = F)
res.MG.down = res.MG.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange <= neg) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(res.MG.down)), 'MG_DOWN_genes.txt', quote = FALSE, col.names = F, row.names = F)


OGRF.gns.up = res.OGRF.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange >= pos) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(OGRF.gns.up)), 'OG1RF_UP_genes.txt', quote = FALSE, col.names = F, row.names = F)
OGRF.gns.down = res.OGRF.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange <= neg) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(OGRF.gns.down)), 'OG1RF_DOWN_genes.txt', quote = FALSE, col.names = F, row.names = F)


skpo.gns.up = res.skpo.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange >= pos) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(skpo.gns.up)), 'skpo_UP_genes.txt', quote = FALSE, col.names = F, row.names = F)
skpo.gns.down = res.skpo.tidy %>% dplyr::filter(padj <= 0.05, log2FoldChange <= neg) %>% dplyr::select(gene_id) %>% t %>% as.vector
write.table(c('Wormbase.ID', unique(skpo.gns.down)), 'skpo_DOWN_genes.txt', quote = FALSE, col.names = F, row.names = F)

res.OP50.tidy %>% 
  filter(gene_id == 'WBGene00000033')







