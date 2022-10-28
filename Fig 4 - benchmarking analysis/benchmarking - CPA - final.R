devtools::install("../smidaSMICellTyping")
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(MLEcell)
library(clue)
library(pheatmap)
library(SingleCellExperiment)
library(scater)
library(mclust)
library(scran)
library(igraph)
library(SC3)
library(Seurat)
library(CHETAH)
rm(list = ls())
pal = tableau_color_pal()(10) #color palette
setwd("~/hades/pdanaher/Insitutype manuscript/Fig 4 - benchmarking in CPA")

#### load CPA data: ---------------
load("../data/Cell Pellet Array annotation and raw counts.RData")
load("../data/CPA16_RNAseq.RData")
counts = t(as.matrix(raw))
rm(raw)
badprobes = read.csv("../data/genes with efficiency 8-fold below average in old CPA panel.csv")[,2]
counts = counts[, !is.element(colnames(counts), badprobes)]
boxplot(log(annot$raw_totalCount)~annot$cell_line,las=2)

abline(h = 150, col = 2)
# remove the cell lines with failed FOVs:
failed.cell.lines = names(which(by(annot$raw_totalCount, annot$cell_line, median) < 150))
remove = is.element(annot$cell_line, failed.cell.lines)
annot = annot[!remove, ]
counts = counts[!remove, ]

#### unsupervised clustering -------------------------
# without cohorts 
unsup.nocohorts.13.10start10k = insitutype(counts = counts,
                                           neg = annot$negmean,
                                           n_clusts = 13,
                                           n_starts = 10,
                                           n_phase1 = 10000)

# if cohorts
ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort = fastCohorting(mat = ifdata)

unsup.ifcohorts.13.10start10k = insitutype(counts = counts,
                                           neg = annot$negmean,
                                           n_clusts = 13,
                                           n_starts = 10,
                                           n_phase1 = 10000,
                                           cohort = ifcohort)

## Preprocess for other methods
sce = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                           colData = annot)
sce = logNormCounts(sce) #normalize
sce = scater::runPCA(sce)
pcs = reducedDim(sce)

## Mclust
mclust_eee_13 = Mclust(pcs, G = 13, modelNames = "EEE")
mclust::adjustedRandIndex(mclust_eee_13$classification, annot$cell_line)

## Leiden
g = buildSNNGraph(sce, use.dimred="PCA", k = 20) #increase k for fewer clusters
leiden.clusts = sapply(c(0.01, 0.03, seq(0.05,0.2,0.05)), function(x){cluster_leiden(g, resolution = x)$membership})
table(leiden.clust$membership)
leiden.ari = mclust::adjustedRandIndex(leiden.clust$membership, annot$cell_line)
leiden.aris = apply(leiden.clusts, 2, function(x){adjustedRandIndex(x, annot$cell_line)})
mclust::adjustedRandIndex(leiden.clusts[,2], annot$cell_line)

sce$leiden = as.factor(leiden.clust$membership)

# k means
km.13 =kmeans(pcs, centers = 13)
mclust::adjustedRandIndex(km.13$cluster, annot$cell_line)

## SC3
rowData(sce)$feature_symbol = rownames(sce)
sce <- sc3(sce, ks = 13)
sce <- sc3_run_svm(sce, ks = 13)
mclust::adjustedRandIndex(sce$sc3_13_clusters, annot$cell_line)

tibble(Truth = annot$cell_line,
       "Insitutype (no cohorts)" = unsup.nocohorts.13.10start10k$clust,
       "Insitutype (IF cohorts)" = unsup.ifcohorts.13.10start10k$clust,
       "Mclust" = mclust_eee_13$classification,
       "Leiden" = leiden.clusts[,2],
       "K-means" = km.13$cluster,
       "SC3" = sce$sc3_13_clusters, 
       transcripts = rowSums(counts)) ->
  data_all

ari_ci = function(data){
  1-(data %>%
       summarise(across(`Insitutype (no cohorts)`:SC3, ~ attr(ari(table(Truth, .x), digits = 7)$ari, "ci"))) %>%
       pivot_longer(everything(), values_to = "CI", names_to = "Method"))$CI
}
getCIs = ari_ci(data_all)
data_all %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "All cells")->
  data_all_all_ari
getCIs = ari_ci(data_all %>% filter(transcripts <=50))  
data_all %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "0-50")->
  data_all_0_50_ari
getCIs = ari_ci(data_all %>% filter(transcripts <=100 & transcripts >50))  
data_all %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "51-100")->
  data_all_51_100_ari
getCIs = ari_ci(data_all %>% filter(transcripts <=200 & transcripts >101))  
data_all %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "101-200")->
  data_all_101_200_ari
getCIs = ari_ci(data_all %>% filter(transcripts >200))  
data_all %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = ">200")->
  data_all_201_ari
rbind(data_all_all_ari, data_all_0_50_ari,
      data_all_51_100_ari, data_all_101_200_ari,
      data_all_201_ari) %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `1-ARI`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(y = "1 - ARI", title = "All genes", x = "Transcripts",
       fill = "Unsupervised\nmethods") +
  scale_fill_manual(values = pal) ->
  plot_all

## Subsample genes
# Half genes
subsample_half = sample.int(ncol(counts), ncol(counts)/2)

# without cohorts 
unsup.nocohorts.13.10k10starts_half = insitutype(counts = counts[,subsample_half],
                                                 neg = annot$negmean,
                                                 n_clusts = 13,
                                                 n_phase1 = 10000,
                                                 n_starts = 10)
mclust::adjustedRandIndex(unsup.nocohorts.13.10k_half$clust, annot$cell_line)

# if cohorts
ifdata = as.matrix(annot[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort = fastCohorting(mat = ifdata)

unsup.ifcohorts.13.10k10starts_half = insitutype(counts = counts[,subsample_half],
                                                 neg = annot$negmean,
                                                 n_clusts = 13,
                                                 n_phase1 = 10000,
                                                 n_starts = 10,
                                                 cohort = ifcohort)

mclust::adjustedRandIndex(unsup.ifcohorts.13.10k_half$clust, annot$cell_line)

## Preprocess
sce_half = SingleCellExperiment(assays = list(counts = t(counts[,subsample_half])), #make SCE
                                colData = annot)
sce_half = logNormCounts(sce_half) #normalize
sce_half = scater::runPCA(sce_half)
pcs_half = reducedDim(sce_half)

## Mclust
mclust_eee_13_half = Mclust(pcs_half, G = 13, modelNames = "EEE")
mclust::adjustedRandIndex(mclust_eee_13_half$classification, annot$cell_line)


## Leiden
g_half = buildSNNGraph(sce_half, use.dimred="PCA", k = 20) #increase k for fewer clusters
leiden.clusts_half = sapply(c(0.01, 0.03, seq(0.05,0.2,0.05)), function(x){cluster_leiden(g_half, resolution = x)$membership})
leiden.aris_half = apply(leiden.clusts_half, 2, function(x){adjustedRandIndex(x, annot$cell_line)})
mclust::adjustedRandIndex(leiden.clusts_half[,3], annot$cell_line)


# k means
km.13_half =kmeans(pcs_half, centers = 13)
mclust::adjustedRandIndex(km.13_half$cluster, annot$cell_line)

## SC3
rowData(sce_half)$feature_symbol = rownames(sce_half)
sce_half <- sc3(sce_half, ks = 13)
sce_half <- sc3_run_svm(sce_half, ks = 13)

mclust::adjustedRandIndex(sce_half$sc3_13_clusters, annot$cell_line)

data.frame(Method = c("Insitutype (no cohorts)",
                      "Insitutype (IF cohorts)",
                      "Mclust",
                      "Leiden",
                      "K-means",
                      "SC3"),
           ARI = c(adjustedRandIndex(annot$cell_line, unsup.nocohorts.13.10k_half$clust),
                   adjustedRandIndex(annot$cell_line, unsup.ifcohorts.13.10k_half$clust),
                   NA,
                   adjustedRandIndex(annot$cell_line, leiden.clusts_half[,3]),
                   adjustedRandIndex(annot$cell_line, km.13_half$cluster),
                   adjustedRandIndex(annot$cell_line, sce_half$sc3_13_clusters)
           )) %>%
  mutate(ARI.transform = -log10(1-ARI)) ->
  ARI.summary_half


tibble(Truth = annot$cell_line,
       "Insitutype (no cohorts)" = unsup.nocohorts.13.10k10starts_half$clust,
       "Insitutype (IF cohorts)" = unsup.ifcohorts.13.10k10starts_half$clust,
       "Mclust" = NA,
       "Leiden" = leiden.clusts_half[,3],
       "K-means" = km.13_half$cluster,
       "SC3" = sce_half$sc3_13_clusters, 
       transcripts = rowSums(counts[,subsample_half])) ->
  data_half
getCIs = ari_ci(data_half %>% select(-Mclust))  
data_half %>%
  select(-Mclust) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:5],
         ci_high = getCIs[6:10],
         Cells = "All cells")->
  data_half_all_ari
getCIs = ari_ci(data_half %>% select(-Mclust) %>% filter(transcripts <=50))  
data_half %>%
  select(-Mclust) %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:5],
         ci_high = getCIs[6:10],
         Cells = "0-50")->
  data_half_0_50_ari
getCIs = ari_ci(data_half  %>% select(-Mclust) %>% filter(transcripts <=100 & transcripts >50))  
data_half %>%
  select(-Mclust) %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:5],
         ci_high = getCIs[6:10],
         Cells = "51-100")->
  data_half_51_100_ari
getCIs = ari_ci(data_half %>% select(-Mclust) %>% filter(transcripts <=200 & transcripts >101))  
data_half %>%
  select(-Mclust) %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:5],
         ci_high = getCIs[6:10],
         Cells = "101-200")->
  data_half_101_200_ari
getCIs = ari_ci(data_half %>% select(-Mclust) %>% filter(transcripts >200))  
data_half %>%
  select(-Mclust) %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:5],
         ci_high = getCIs[6:10],
         Cells = ">200")->
  data_half_201_ari
rbind(data_half_all_ari, data_half_0_50_ari,
      data_half_51_100_ari, data_half_101_200_ari,
      data_half_201_ari) %>%
  filter(Method != "Mclust") %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `1-ARI`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(y = "1 - ARI", title = "Half genes", x = "Transcripts") +
  guides(fill = "none") + 
  scale_fill_manual(values = pal[c(1:4,6)]) ->
  plot_half


ggplot(ARI.summary_half, aes(x = Method, y = ARI.transform, fill = Method))+
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "-log(1-ARI)") +
  #coord_cartesian(ylim = c(0.8, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(ARI.summary_half, aes(x = Method, y = ARI, fill = Method))+
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Adjusted Rand Index (ARI)") +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# 96 HVGs
## Preprocess
sce_hvg = SingleCellExperiment(assays = list(counts = t(counts)), #make SCE
                               colData = annot)
sce_hvg = logNormCounts(sce_hvg) #normalize
top = getTopHVGs(sce_hvg, n = 96)
sce_hvg = sce_hvg[top, ]
sce_hvg = scater::runPCA(sce_hvg, subset_row = top)

counts_hvg = counts[,top]
hvg_filter = rowSums(counts_hvg) >= 5
annot_hvg = annot[hvg_filter,]
counts_hvg = counts_hvg[hvg_filter,]
sce_hvg = sce_hvg[,hvg_filter]
pcs_hvg = reducedDim(sce_hvg)

# without cohorts 
unsup.nocohorts.13.10k10starts_hvg = insitutype(counts = counts_hvg,
                                                neg = annot_hvg$negmean,
                                                n_clusts = 13,
                                                n_phase1 = 10000,
                                                n_starts = 10)
mclust::adjustedRandIndex(unsup.nocohorts.13.10k_hvg$clust, annot_hvg$cell_line)

# if cohorts
ifdata_hvg = as.matrix(annot_hvg[, paste0(c("Blue", "Red", "Yellow", "Green"), "Int")])
ifcohort_hvg = fastCohorting(mat = ifdata_hvg)

unsup.ifcohorts.13.10k10starts_hvg = insitutype(counts = counts_hvg,
                                                neg = annot_hvg$negmean,
                                                n_clusts = 13,
                                                n_phase1 = 10000,
                                                n_starts = 10,
                                                cohort = ifcohort_hvg)

mclust::adjustedRandIndex(unsup.ifcohorts.13.10k_hvg$clust, annot_hvg$cell_line)

## Mclust
mclust_eee_13_hvg = Mclust(pcs_hvg, G = 13, modelNames = "EEE")
mclust::adjustedRandIndex(mclust_eee_13_hvg$classification, annot_hvg$cell_line)

## Leiden
g_hvg = buildSNNGraph(sce_hvg, use.dimred="PCA", k = 20) #increase k for fewer clusters
leiden.clusts_hvg = sapply(c(0.01, 0.03, seq(0.05,0.2,0.05)), function(x){cluster_leiden(g_hvg, resolution = x)$membership})
leiden.aris_hvg = apply(leiden.clusts_hvg, 2, function(x){adjustedRandIndex(x, annot_hvg$cell_line)})
mclust::adjustedRandIndex(leiden.clusts_hvg[,3], annot_hvg$cell_line)


# k means
km.13_hvg =kmeans(pcs_hvg, centers = 13)
mclust::adjustedRandIndex(km.13_hvg$cluster, annot_hvg$cell_line)

## SC3
rowData(sce_hvg)$feature_symbol = rownames(sce_hvg)
sce_hvg <- sc3(sce_hvg, ks = 13)
sce_hvg <- sc3_run_svm(sce_hvg, ks = 13)

mclust::adjustedRandIndex(sce_hvg$sc3_13_clusters, annot_hvg$cell_line)

data.frame(Method = c("Insitutype (no cohorts)",
                      "Insitutype (IF cohorts)",
                      "Mclust",
                      "Leiden",
                      "K-means",
                      "SC3"),
           ARI = c(adjustedRandIndex(annot_hvg$cell_line, unsup.nocohorts.13.10k_hvg$clust),
                   adjustedRandIndex(annot_hvg$cell_line, unsup.ifcohorts.13.10k_hvg$clust),
                   adjustedRandIndex(annot_hvg$cell_line, mclust_eee_13_hvg$classification),
                   adjustedRandIndex(annot_hvg$cell_line, leiden.clusts_hvg[,3]),
                   adjustedRandIndex(annot_hvg$cell_line, km.13_hvg$cluster),
                   adjustedRandIndex(annot_hvg$cell_line, sce_hvg$sc3_13_clusters)
           )) %>%
  mutate(ARI.transform = -log10(1-ARI)) ->
  ARI.summary_hvg

tibble(Truth = annot_hvg$cell_line,
       "Insitutype (no cohorts)" = unsup.nocohorts.13.10k10starts_hvg$clust,
       "Insitutype (IF cohorts)" = unsup.ifcohorts.13.10k10starts_hvg$clust,
       "Mclust" = mclust_eee_13_hvg$classification,
       "Leiden" = leiden.clusts_hvg[,3],
       "K-means" = km.13_hvg$cluster,
       "SC3" = sce_hvg$sc3_13_clusters, 
       transcripts = rowSums(counts_hvg)) ->
  data_hvg

getCIs = ari_ci(data_hvg)  
data_hvg %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "All cells")->
  data_hvg_all_ari
getCIs = ari_ci(data_hvg %>% filter(transcripts <=50))  
data_hvg %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "0-50")->
  data_hvg_0_50_ari
getCIs = ari_ci(data_hvg %>% filter(transcripts <=100 & transcripts >50))  
data_hvg %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "51-100")->
  data_hvg_51_100_ari
getCIs = ari_ci(data_hvg %>% filter(transcripts <=200 & transcripts >101))  
data_hvg %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = "101-200")->
  data_hvg_101_200_ari
getCIs = ari_ci(data_hvg %>% filter(transcripts >200))  
data_hvg %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:SC3, ~ adjustedRandIndex(Truth, .x))) %>%
  pivot_longer(everything(), values_to = "ARI", names_to = "Method") %>%
  mutate(`1-ARI` = 1-ARI,
         Method = relevel(factor(Method), ref = "Insitutype (no cohorts)"),
         ci_low = getCIs[1:6],
         ci_high = getCIs[7:12],
         Cells = ">200")->
  data_hvg_201_ari
rbind(data_hvg_all_ari, data_hvg_0_50_ari,
      data_hvg_51_100_ari, data_hvg_101_200_ari,
      data_hvg_201_ari) %>%
  # filter(Method != "Mclust") %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `1-ARI`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(y = "1 - ARI", title = "10% highly variable genes", x = "Transcripts") +
  guides(fill = "none") + 
  scale_fill_manual(values = pal[c(1:6)]) ->
  plot_hvg
plot_hvg
p_unsup = (plot_all / plot_half / plot_hvg) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')

ggplot(ARI.summary_hvg, aes(x = Method, y = ARI.transform, fill = Method))+
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "-log(1-ARI)") +
  #coord_cartesian(ylim = c(0.8, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(ARI.summary_hvg, aes(x = Method, y = ARI, fill = Method))+
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Adjusted Rand Index (ARI)") +
  coord_cartesian(ylim = c(0.8, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#### supervised cell typing -----------------------
alpha_ref = read_csv("../data/A1182_S1_CellMatrix.csv")
alpha_ref = alpha_ref[,]
alpha_ref %>%
  select(-`...1`) %>%
  remove_rownames() %>%
  column_to_rownames("CellId") %>%
  select(!contains("FalseCode"),
         !contains("NegPrb")) ->
  alpha_ref_filtered
fov_ref = sapply(strsplit(rownames(alpha_ref_filtered), "_"),"[[", 2)
fov_key = data.frame(FOV = paste0("FOV", stringr::str_pad(1:16, 3, pad = "0")),
                     cell_line = c("IGROV1", "MOLT4", "HOP92", "HS578T",
                                   "COLO205", "CCRF-CEM", "RPMI-8226", "DU145", 
                                   "PC3", "M14", "SUDHL4", "HCT116", 
                                   "EKVX", "HL60", "SKMEL2", "MDA-MB-468"))
annot_ref = left_join(data.frame(FOV = fov_ref),
                      fov_key) 
cell_filter = annot_ref$cell_line %in% annot$cell_line
annot_ref = annot_ref[cell_filter,]
alpha_ref %>%
  select(contains("NegPrb")) %>%
  rowMeans() ->
  neg_ref
neg_ref = neg_ref[cell_filter]
gene_intersect = intersect(colnames(alpha_ref_filtered), colnames(counts))
counts_ref = as.matrix(alpha_ref_filtered[cell_filter,gene_intersect])
# get reference:
meanprofiles = MLEcell:::Estep(counts = counts_ref, clust = annot_ref$cell_line, neg = neg_ref)


# without cohorts 
sup.nocohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = meanprofiles,
                             cohort = NULL)$clust

# with IF cohorts
sup.ifcohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = meanprofiles, 
                             cohort = ifcohort)$clust
table(sup.nocohorts == annot$cell_line)
table(sup.ifcohorts == annot$cell_line)
table(sup.ifcohorts, annot$cell_line)

#SingleR
library(SingleR)

singleR = SingleR(test = sce,
                  ref = meanprofiles,
                  labels = colnames(meanprofiles))
table(singleR$pruned.labels == annot$cell_line)

singleR_half = SingleR(test = sce_half,
                       ref = meanprofiles,
                       labels = colnames(meanprofiles))
singleR_hvg = SingleR(test = sce_hvg,
                      ref = meanprofiles,
                      labels = colnames(meanprofiles))

#Seurat
library(Seurat)
rownames(annot) = rownames(counts)
seurat.test = CreateSeuratObject(t(counts), meta.data = annot)
rownames(annot_ref) = rownames(counts_ref)
seurat.ref = CreateSeuratObject(t(counts_ref), meta.data = annot_ref)
seurat.test = NormalizeData(seurat.test, verbose = FALSE)
seurat.ref = NormalizeData(seurat.ref, verbose = FALSE)
anchors <- FindTransferAnchors(reference = seurat.ref, query = seurat.test,
                               features = gene_intersect)
seurat.predictions <- TransferData(anchorset = anchors, refdata = seurat.ref$cell_line)
table(seurat.predictions$predicted.id == annot$cell_line)

#half genes 
gene_intersect_half = intersect(colnames(alpha_ref_filtered), colnames(counts)[subsample_half])

seurat.test.half = CreateSeuratObject(t(counts[,subsample_half]), meta.data = annot)
seurat.test.half = NormalizeData(seurat.test.half, verbose = FALSE)
anchors_half <- FindTransferAnchors(reference = seurat.ref, query = seurat.test.half,
                                    features = gene_intersect_half)
seurat.predictions.half <- TransferData(anchorset = anchors_half, refdata = seurat.ref$cell_line)

#hvg
gene_intersect_hvg = intersect(colnames(alpha_ref_filtered), colnames(counts_hvg))
seurat.test.hvg = CreateSeuratObject(t(counts_hvg), meta.data = annot_hvg)
seurat.test.hvg = NormalizeData(seurat.test.hvg, verbose = FALSE)
anchors_hvg <- FindTransferAnchors(reference = seurat.ref, query = seurat.test.hvg,
                                   features = gene_intersect_hvg)
seurat.predictions.hvg <- TransferData(anchorset = anchors_hvg, refdata = seurat.ref$cell_line)


#CHETAH
## For the reference we define a "counts" assay and "celltypes" metadata
sce_ref = SingleCellExperiment(assays = list(counts = t(counts_ref)),
                               colData = annot_ref)
sce_ref = sce_ref[,colSums(counts(sce_ref))>0]
sce_ref = logNormCounts(sce_ref)
## For the input we define a "counts" assay and "TSNE" reduced dimensions
sce = runTSNE(sce)
sce <- CHETAHclassifier(input = sce,
                        ref_cells = sce_ref,
                        ref_ct = "cell_line")
table(sce$celltype_CHETAH == annot$cell_line)

sce_half = runTSNE(sce_half)
sce_half <- CHETAHclassifier(input = sce_half,
                             ref_cells = sce_ref,
                             ref_ct = "cell_line")

sce_ref_hvg = sce_ref[intersect(rownames(sce_ref), top),]
sce_ref_hvg = sce_ref[,colSums(counts(sce_ref_hvg))>0]
sce_ref_hvg = logNormCounts(sce_ref_hvg)
sce_hvg = runTSNE(sce_hvg)
sce_hvg <- CHETAHclassifier(input = sce_hvg,
                            ref_cells = sce_ref_hvg,
                            ref_ct = "cell_line")

## subsample
sup.nocohorts.half = insitutypeML(counts = counts[,subsample_half],
                                  neg = annot$negmean,
                                  reference_profiles = meanprofiles,
                                  cohort = NULL)$clust

# with IF cohorts
sup.ifcohorts.half = insitutypeML(counts = counts[,subsample_half],
                                  neg = annot$negmean,
                                  reference_profiles = meanprofiles, 
                                  cohort = ifcohort)$clust
# hvg
sup.nocohorts.hvg = insitutypeML(counts = counts_hvg,
                                 neg = annot_hvg$negmean,
                                 reference_profiles = meanprofiles,
                                 cohort = NULL)$clust

# with IF cohorts
sup.ifcohorts.hvg = insitutypeML(counts = counts_hvg,
                                 neg = annot_hvg$negmean,
                                 reference_profiles = meanprofiles, 
                                 cohort = ifcohort_hvg)$clust
table(sup.nocohorts == annot$cell_line)
table(sup.ifcohorts == annot$cell_line)
table(sup.nocohorts.half == annot$cell_line)
table(sup.ifcohorts.half == annot$cell_line)
table(sup.nocohorts.hvg == annot$cell_line)
table(sup.ifcohorts.hvg == annot$cell_line)

#Plot
data.frame(Method = c("Insitutype (no cohorts)",
                      "Insitutype (IF cohorts)",
                      "SingleR",
                      "Seurat",
                      "CHETAH"),
           Accuracy = c(sum(sup.nocohorts == annot$cell_line)/nrow(annot)*100,
                        sum(sup.ifcohorts == annot$cell_line)/nrow(annot)*100,
                        sum(singleR$pruned.labels == annot$cell_line, na.rm = T)/nrow(annot)*100,
                        sum(seurat.predictions$predicted.id == annot$cell_line)/nrow(annot)*100,
                        sum(sce$celltype_CHETAH == annot$cell_line)/nrow(annot)*100)
) ->
  accuracy_summary

tibble(Truth = annot$cell_line,
       "Insitutype (no cohorts)" = sup.nocohorts,
       "Insitutype (IF cohorts)" = sup.ifcohorts,
       "SingleR" = singleR$pruned.labels,
       "Seurat" = seurat.predictions$predicted.id,
       "CHETAH" = sce$celltype_CHETAH,
       transcripts = rowSums(counts)) ->
  sup_all


sup_all %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_all %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "All cells")->
  sup_all_all_ari
sup_all %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_all %>% filter(transcripts<=50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "0-50")->
  sup_all_0_50_ari
sup_all %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_all %>% filter(transcripts <=100 & transcripts >50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "51-100")->
  sup_all_51_100_ari
sup_all %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_all %>% filter(transcripts <=200 & transcripts >101) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "101-200")->
  sup_all_101_200_ari
sup_all %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_all %>% filter(transcripts >200) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = ">200")->
  sup_all_201_ari
rbind(sup_all_all_ari, sup_all_0_50_ari,
      sup_all_51_100_ari, sup_all_101_200_ari,
      sup_all_201_ari) %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `Error (%)`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(title = "All genes", x = "Transcripts", fill = "Supervised\nmethods") +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(values = pal[c(1:2, 7:9)]) ->
  plot_all_sup
plot_all_sup
#plot half
tibble(Truth = annot$cell_line,
       "Insitutype (no cohorts)" = sup.nocohorts.half,
       "Insitutype (IF cohorts)" = sup.ifcohorts.half,
       "SingleR" = singleR_half$pruned.labels,
       "Seurat" = seurat.predictions.half$predicted.id,
       "CHETAH" = sce_half$celltype_CHETAH,
       transcripts = rowSums(counts)) ->
  sup_half

sup_half %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_half %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "All cells")->
  sup_half_all_ari
sup_half %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_half %>% filter(transcripts <=50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "0-50")->
  sup_half_0_50_ari
sup_half %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_half %>% filter(transcripts <=100 & transcripts >50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "51-100")->
  sup_half_51_100_ari
sup_half %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_half %>% filter(transcripts <=200 & transcripts >101) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "101-200")->
  sup_half_101_200_ari
sup_half %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_half %>% filter(transcripts >200) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = ">200")->
  sup_half_201_ari
rbind(sup_half_all_ari, sup_half_0_50_ari,
      sup_half_51_100_ari, sup_half_101_200_ari,
      sup_half_201_ari) %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `Error (%)`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(title = "Half genes", x = "Transcripts") +
  guides(fill = "none") + 
  scale_fill_manual(values = pal[c(1:2, 7:9)]) ->
  plot_half_sup
plot_half_sup

#plot hvg
tibble(Truth = annot_hvg$cell_line,
       "Insitutype (no cohorts)" = sup.nocohorts.hvg,
       "Insitutype (IF cohorts)" = sup.ifcohorts.hvg,
       "SingleR" = singleR_hvg$pruned.labels,
       "Seurat" = seurat.predictions.hvg$predicted.id,
       "CHETAH" = NA,
       transcripts = rowSums(counts_hvg)) ->
  sup_hvg

sup_hvg %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_hvg %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "All cells")->
  sup_hvg_all_ari
sup_hvg %>%
  filter(transcripts <=50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_hvg %>% filter(transcripts <=50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "0-50")->
  sup_hvg_0_50_ari
sup_hvg %>%
  filter(transcripts <=100 & transcripts >50) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_hvg %>% filter(transcripts <=100 & transcripts >50) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "51-100")->
  sup_hvg_51_100_ari
sup_hvg %>%
  filter(transcripts <=200 & transcripts >101) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_hvg %>% filter(transcripts <=200 & transcripts >101) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = "101-200")->
  sup_hvg_101_200_ari
sup_hvg %>%
  filter(transcripts >200) %>%
  summarise(across(`Insitutype (no cohorts)`:CHETAH, ~ sum(Truth == .x, na.rm = T)/n()*100)) %>%
  pivot_longer(everything(), values_to = "Accuracy", names_to = "Method") %>%
  mutate(`Error (%)` = 100-Accuracy,
         Method = factor(Method, levels = c("Insitutype (no cohorts)",
                                            "Insitutype (IF cohorts)",
                                            "CHETAH",
                                            "Seurat",
                                            "SingleR")),
         n = (sup_hvg %>% filter(transcripts >200) %>% summarise(n = n()))$n,
         wrong = `Error (%)` * n/100,
         ci_low = binom::binom.confint(wrong, n, methods = "wilson")$lower*100,
         ci_high = binom::binom.confint(wrong, n, methods = "wilson")$upper*100,
         Cells = ">200")->
  sup_hvg_201_ari
rbind(sup_hvg_all_ari, sup_hvg_0_50_ari,
      sup_hvg_51_100_ari, sup_hvg_101_200_ari,
      sup_hvg_201_ari) %>%
  filter(Method != "CHETAH") %>%
  mutate(Cells = factor(Cells, levels = c("All cells", "0-50", "51-100", "101-200", ">200"))) %>%
  ggplot(aes(x = Cells, y = `Error (%)`, fill = Method))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, position = position_dodge(0.9))+
  theme_bw() +
  labs(title = "10% highly variable genes", x = "Transcripts") +
  guides(fill = "none") + 
  scale_fill_manual(values = pal[c(1:2, 8:9)]) ->
  plot_hvg_sup
plot_hvg_sup
p_sup = (plot_all_sup / plot_half_sup / plot_hvg_sup) + 
  plot_layout(guides = "collect")  & 
  theme(legend.position = 'bottom')
p_final = (p_unsup | p_sup) + plot_annotation(tag_levels = "a")
ggsave(filename = "benchmark_all_10-10.png", p_final & theme(text = element_text(size = 10)), 
       width = 9, height = 10, device = png)
