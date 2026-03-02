#Bartek Puzio 1146247 Assignment 2
# Bulk RNA-seq Analysis: Velum Formation in S. cerevisiae
# This script analyzes differential gene expression across three timepoints
# (Early, Thin, Mature) during velum biofilm formation in S. cerevisiae.
# Salmon quantifications are imported via tximport and analyzed with DESeq2.
# Pairwise contrasts and LRT identify genes responding to biofilm progression.
# Functional annotation is performed using GO enrichment, KEGG pathway
# analysis, and GSEA to characterize the biological processes driving
# velum maturation.

library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(styler)
#style_file("genome2.R")

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

#### DATA LOADING ####

quant_files <- list.files("salmon_output", pattern = "quant.sf", 
                          recursive = TRUE, full.names = TRUE)

meta <- data.frame(
  sample     = basename(dirname(quant_files)),
  quant_path = quant_files
)

# Assign timepoint factor with ordered levels
meta$timepoint <- factor(c("Mature", "Mature", "Mature", "Thin", "Thin", "Thin", "Early", "Early", "Early"),
                         levels = c("Early", "Thin", "Mature"))
write.csv(meta, "sample_info.csv", row.names = FALSE)
meta <- meta[nrow(meta):1, ]

# Build transcript-to-gene map from GTF
txdb    <- makeTxDbFromGFF("GCF_000146045.2_R64_genomic.gtf", format = "gtf")
tx2gene <- AnnotationDbi::select(txdb,
                                 keys    = keys(txdb, keytype = "TXNAME"),
                                 columns = c("TXNAME", "GENEID"),
                                 keytype = "TXNAME")
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)

# Import Salmon quantifications
quant_paths <- setNames(meta$quant_path, meta$sample)

txi <- tximport(
  quant_paths,
  type           = "salmon",
  tx2gene        = tx2gene,
  ignoreAfterBar = TRUE
)

#### DESEQ2 OBJECT & GENE MAPPING ####

dds <- DESeqDataSetFromTximport(
  txi,
  colData = meta,
  design  = ~ timepoint
)

# Filter lowly expressed genes: require >= 10 counts in at least 3 samples
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)

# Map ORF IDs to gene names and Entrez IDs
all_gene_ids <- rownames(results(dds))

gene_map <- bitr(all_gene_ids,
                 fromType = "ORF",
                 toType   = c("GENENAME", "ENTREZID"),
                 OrgDb    = org.Sc.sgd.db)

#### PAIRWISE DEG RESULTS ####

ref        <- levels(meta$timepoint)[1]
timepoints <- levels(meta$timepoint)[-1]

# Run shrinkage-corrected pairwise contrasts against Early
deg_list <- lapply(timepoints, function(tp) {
  res <- results(dds, contrast = c("timepoint", tp, ref), alpha = 0.05)
  res <- lfcShrink(dds, contrast = c("timepoint", tp, ref), res = res, type = "normal")
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj) %>%
    mutate(comparison = paste0(tp, "_vs_", ref))
})

# Attach gene names to each contrast
deg_list <- lapply(deg_list, function(df) {
  left_join(df, gene_map, by = c("gene_id" = "ORF")) %>%
    mutate(label = ifelse(is.na(GENENAME), gene_id, GENENAME))
})

names(deg_list) <- paste0(timepoints, "_vs_", ref)

# Save DEG tables
for (nm in names(deg_list)) {
  write.csv(deg_list[[nm]], file = paste0("results/DEG_", nm, ".csv"), row.names = FALSE)
}

#### EXPLORATORY PLOTS ####
# PCA
vsd      <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "timepoint", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = timepoint, label = name)) +
  geom_point(size = 4) +
  geom_text_repel() +
  labs(
    title = "PCA — Velum Formation Over Time",
    x     = paste0("PC1: ", pct_var[1], "% variance"),
    y     = paste0("PC2: ", pct_var[2], "% variance")
  ) +
  theme_classic()

ggsave("plots/PCA_velum_timepoints.pdf", width = 7, height = 5)

# MA plot: overall expression vs fold change for Mature vs Early
pdf("plots/MA_Mature_vs_Early.pdf", width = 7, height = 5)
plotMA(results(dds, contrast = c("timepoint", "Mature", "Early")),
       ylim = c(-5, 5),
       main = "MA Plot: Mature vs Early")
dev.off()

pdf("plots/MA_Thin_vs_Early.pdf", width = 7, height = 5)
plotMA(results(dds, contrast = c("timepoint", "Thin", "Early")),
       ylim = c(-5, 5),
       main = "MA Plot: Thin vs Early")
dev.off()

#### VOLCANO PLOTS ####

# Thin vs Early
df1 <- deg_list[["Thin_vs_Early"]] %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

ggplot(df1, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  ggrepel::geom_text_repel(
    data = filter(df1, sig != "NS") %>% slice_min(padj, n = 15),
    aes(label = label), size = 3, max.overlaps = 5
  ) +
  scale_color_manual(values = c(Up = "#e63946", Down = "#457b9d", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Thin vs Early", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_classic()

ggsave("plots/volcano_Thin_vs_Early.pdf", width = 7, height = 5)

# Mature vs Early
df2 <- deg_list[["Mature_vs_Early"]] %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

ggplot(df2, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  ggrepel::geom_text_repel(
    data = filter(df2, sig != "NS") %>% slice_min(padj, n = 15),
    aes(label = label), size = 3, max.overlaps = 5
  ) +
  scale_color_manual(values = c(Up = "#e63946", Down = "#457b9d", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Mature vs Early", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_classic()

ggsave("plots/volcano_Mature_vs_Early.pdf", width = 7, height = 5)

# Thin vs Mature (computed separately as it's not in deg_list)
df3 <- results(dds, contrast = c("timepoint", "Thin", "Mature"), alpha = 0.05)
df3 <- lfcShrink(dds, contrast = c("timepoint", "Thin", "Mature"), res = df3, type = "normal")
df3 <- as.data.frame(df3) %>%
  tibble::rownames_to_column("gene_id") %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  ))

df3 <- left_join(df3, gene_map, by = c("gene_id" = "ORF")) %>%
  mutate(label = ifelse(is.na(GENENAME), gene_id, GENENAME))

ggplot(df3, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  ggrepel::geom_text_repel(
    data = filter(df3, sig != "NS") %>% slice_min(padj, n = 15),
    aes(label = label), size = 3, max.overlaps = 5
  ) +
  scale_color_manual(values = c(Up = "#e63946", Down = "#457b9d", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Thin vs Mature", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_classic()

ggsave("plots/volcano_Thin_vs_Mature.pdf", width = 7, height = 5)

write.csv(df3, "results/DEG_Thin_vs_Mature.csv", row.names = FALSE)


#### HEATMAP (PAIRWISE DEGs) ####

# Take top 30 DEGs per contrast by adjusted p-value
top_genes <- bind_rows(deg_list) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  group_by(comparison) %>%
  slice_min(padj, n = 30) %>%
  pull(gene_id) %>%
  unique()

mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))

annotation_col <- data.frame(Timepoint = meta$timepoint, row.names = meta$sample)

row_labels2 <- data.frame(gene_id = rownames(mat)) %>%
  left_join(gene_map, by = c("gene_id" = "ORF")) %>%
  mutate(label = ifelse(is.na(GENENAME), gene_id, GENENAME)) %>%
  pull(label)

pheatmap(
  mat,
  annotation_col = annotation_col,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  labels_row     = row_labels2,
  fontsize_row   = 7,
  cellwidth      = 40,
  cellheight     = 7,
  color          = colorRampPalette(c("#457b9d", "white", "#e63946"))(100)
)

#### LRT (LIKELIHOOD RATIO TEST) ####

# LRT tests which genes vary significantly across all timepoints simultaneously,
# complementing the pairwise contrasts above
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)

lrt_results <- results(dds_lrt) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

write.csv(lrt_results, "results/LRT_timepoint.csv", row.names = FALSE)

sig_lrt <- filter(lrt_results, padj < 0.05)
cat("Significant genes across timepoints:", nrow(sig_lrt), "\n")

# Heatmap of top 50 LRT genes
top_lrt <- sig_lrt %>% slice_min(padj, n = 50) %>% pull(gene_id)

mat_lrt <- assay(vsd)[top_lrt, ]
mat_lrt <- t(scale(t(mat_lrt)))

row_labels <- data.frame(gene_id = rownames(mat_lrt)) %>%
  left_join(gene_map, by = c("gene_id" = "ORF")) %>%
  mutate(label = ifelse(is.na(GENENAME), gene_id, GENENAME)) %>%
  pull(label)

pheatmap(
  mat_lrt,
  annotation_col = annotation_col,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  labels_row     = row_labels,
  fontsize_row   = 7,
  color          = colorRampPalette(c("#457b9d", "white", "#e63946"))(100)
)

#### FUNCTIONAL ENRICHMENT ####

# Gene lists derived from Mature vs Early contrast
res_df <- as.data.frame(deg_list[["Mature_vs_Early"]]) %>%
  left_join(gene_map, by = c("gene_id" = "ORF"))

up_genes   <- res_df %>% filter(padj < 0.05 & log2FoldChange >  1) %>% pull(gene_id) %>% unique()
down_genes <- res_df %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(gene_id) %>% unique()
sig_genes  <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_id) %>% unique()
all_genes  <- res_df$gene_id %>% unique()
cat("Up:", length(up_genes), "| Down:", length(down_genes), "| Sig total:", length(sig_genes), "\n")

# GO ORA — Biological Process, Molecular Function, Cellular Component
ego_bp <- enrichGO(gene          = sig_genes,
                   universe      = all_genes,
                   OrgDb         = org.Sc.sgd.db,
                   keyType       = "ORF",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_mf <- enrichGO(gene          = sig_genes,
                   universe      = all_genes,
                   OrgDb         = org.Sc.sgd.db,
                   keyType       = "ORF",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_cc <- enrichGO(gene          = sig_genes,
                   universe      = all_genes,
                   OrgDb         = org.Sc.sgd.db,
                   keyType       = "ORF",
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

# Prepare gene lists for Thin vs Early
res_thin_early <- as.data.frame(deg_list[["Thin_vs_Early"]]) %>%
  left_join(gene_map, by = c("gene_id" = "ORF"))

sig_genes_te  <- res_thin_early %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_id) %>% unique()
all_genes_te  <- res_thin_early$gene_id %>% unique()

ego_bp_te <- enrichGO(gene          = sig_genes_te,
                      universe      = all_genes_te,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = "ORF",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# Prepare gene lists for Thin vs Mature
res_thin_mature <- df3 %>%
  left_join(gene_map, by = c("gene_id" = "ORF"))

sig_genes_tm  <- res_thin_mature %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_id) %>% unique()
all_genes_tm  <- res_thin_mature$gene_id %>% unique()

ego_bp_tm <- enrichGO(gene          = sig_genes_tm,
                      universe      = all_genes_tm,
                      OrgDb         = org.Sc.sgd.db,
                      keyType       = "ORF",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# Dotplots
dotplot(ego_bp_te, showCategory = 15,
        title = "GO Biological Process Enrichment (Thin vs Early)") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_GO_BP_Thin_vs_Early.pdf", width = 8, height = 6)

dotplot(ego_bp_tm, showCategory = 15,
        title = "GO Biological Process Enrichment (Thin vs Mature)") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_GO_BP_Thin_vs_Mature.pdf", width = 8, height = 6)

# Export
if (nrow(as.data.frame(ego_bp_te)) > 0) write.csv(as.data.frame(ego_bp_te), "results/GO_BP_Thin_vs_Early.csv")
if (nrow(as.data.frame(ego_bp_tm)) > 0) write.csv(as.data.frame(ego_bp_tm), "results/GO_BP_Thin_vs_Mature.csv")
# Compare GO enrichment between up and downregulated gene sets
compare_go <- compareCluster(
  geneCluster   = list(Upregulated = up_genes, Downregulated = down_genes),
  fun           = "enrichGO",
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# KEGG pathway enrichment
sig_entrez <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% pull(gene_id) %>% unique()
all_entrez <- res_df %>% pull(gene_id) %>% unique()

kegg_enrich <- enrichKEGG(gene         = sig_entrez,
                          universe     = all_entrez,
                          organism     = "sce",
                          keyType      = "kegg",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
# KEGG enrichment — Thin vs Early
sig_entrez_te <- res_thin_early %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(gene_id) %>% unique()

all_entrez_te <- res_thin_early %>%
  pull(gene_id) %>% unique()

kegg_te <- enrichKEGG(gene         = sig_entrez_te,
                      universe     = all_entrez_te,
                      organism     = "sce",
                      keyType      = "kegg",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

dotplot(kegg_te, showCategory = 15,
        title = "KEGG Pathway Enrichment (Thin vs Early)") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_KEGG_Thin_vs_Early.pdf", width = 8, height = 6)

if (nrow(as.data.frame(kegg_te)) > 0) write.csv(as.data.frame(kegg_te), "results/KEGG_Thin_vs_Early.csv")

# GSEA ranked by log2FC
gene_ranks <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  { setNames(.$log2FoldChange, .$gene_id) }

gsea_go <- gseGO(geneList     = gene_ranks,
                 OrgDb        = org.Sc.sgd.db,
                 keyType      = "ORF",
                 ont          = "BP",
                 minGSSize    = 15,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05)

#### ENRICHMENT PLOTS #### 

dotplot(ego_bp, showCategory = 15,
        title = "GO Biological Process Enrichment (Mature vs Early)") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_GO_BP.pdf", width = 8, height = 6)

dotplot(compare_go, showCategory = 10,
        title = "GO BP: Upregulated vs Downregulated in Mature Biofilm") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_compare_GO.pdf", width = 8, height = 6)

dotplot(kegg_enrich, showCategory = 15,
        title = "KEGG Pathway Enrichment (Mature vs Early)") +
  theme_minimal(base_size = 12)
ggsave("plots/dotplot_KEGG.pdf", width = 8, height = 6)

#### GENE OF INTEREST EXPRESSION ####

# Look up ORF IDs by gene name
gene_map %>% filter(GENENAME %in% c("FLO11", "TDH1", "OLE1"))

genes_of_interest <- c("YIR019C", "YJL052W", "YGL055W")  # FLO11, TDH1, OLE1

# Extract normalized counts per sample for each gene
plot_df <- purrr::map_dfr(genes_of_interest, function(g) {
  d <- plotCounts(dds, gene = g, intgroup = "timepoint", returnData = TRUE)
  d$gene_id <- g
  d
})

ggplot(plot_df, aes(x = timepoint, y = count, color = timepoint)) +
  geom_point(position = position_jitter(width = 0.12), size = 2) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black") +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
  facet_wrap(~ gene_id, scales = "free_y") +
  scale_color_manual(values = c(Early = "#457b9d", Thin = "grey60", Mature = "#e63946")) +
  theme_classic() +
  xlab("Timepoint") +
  ylab("Normalized expression")

ggsave("plots/GOI_expression_timepoints.pdf", width = 8, height = 4)

#### EXPORT RESULTS ####

if (nrow(as.data.frame(ego_bp)) > 0)      write.csv(as.data.frame(ego_bp),      "results/GO_BP_enrichment.csv")
if (nrow(as.data.frame(ego_mf)) > 0)      write.csv(as.data.frame(ego_mf),      "results/GO_MF_enrichment.csv")
if (nrow(as.data.frame(ego_cc)) > 0)      write.csv(as.data.frame(ego_cc),      "results/GO_CC_enrichment.csv")
if (nrow(as.data.frame(kegg_enrich)) > 0) write.csv(as.data.frame(kegg_enrich), "results/KEGG_enrichment.csv")

if (nrow(as.data.frame(gsea_go)) > 0)     write.csv(as.data.frame(gsea_go),     "results/GSEA_GO_BP.csv")





