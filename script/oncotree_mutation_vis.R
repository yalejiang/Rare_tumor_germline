#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: oncotree_mutation_vis.R
# Description: Analysis and Visualization Pipeline for Rare Tumor Germline Mutation
# Version: 1.0
# Dependencies: See README.md for R package requirements
# Input Files: Figure1.rds, figure2.rds, figure4.rdata, gene_info.xlsx
# Output: Four core figures for rare tumor germline mutation analysis (tiff, 300dpi)
# Note: All input files should be placed in the project root directory
# ==============================================================================

# Clean environment & set seed for reproducibility
rm(list = ls())
set.seed(12345) # Fix random seed for plot consistency

# ==============================================================================
# 1. Load required R packages
# ==============================================================================
required_packages <- c(
  "tidyverse", "readxl", "ggprism", "viridis", "scico", "paletteer",
  "ggsci", "RColorBrewer", "patchwork", "scales", "purrr", "writexl",
  "maftools", "ComplexHeatmap", "circlize", "grid", "org.Hs.eg.db",
  "clusterProfiler", "conflicted", "ggsankey" # ggsankey for dot_sankey function
)

# Install packages if not installed
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    library(pkg, character.only = TRUE)
  }
}))

# Resolve function conflicts (dplyr priority)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")

# ==============================================================================
# Figure 1: Cancer Type Distribution (Pie + Bubble + Bar Plot)
# Core: Sample distribution by cancer type/subtype & mutation rate
# Input: Figure1.rds (contains clin_info, cancer_level, total_counts)
# Output: figure1_cancer_distribution.tiff
# ==============================================================================
load(file = "figur1.rdata") # Load clin_info, cancer_level, total_counts

# Prepare pie chart data
pie_data <- clin_info %>%
  count(Cancer_Type) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(prop)) %>%
  mutate(
    xend = case_when(
      prop > 0.35 ~ 1.6,
      prop > 0.25 ~ 1.9,
      prop > 0.05 ~ 2.2,
      prop > 0.04 ~ 2.4,
      prop > 0.03 ~ 2.6,
      prop > 0.02 ~ 2.8,
      TRUE ~ 3
    )
  ) %>%
  mutate(Cancer_Type = factor(Cancer_Type, levels = Cancer_Type)) %>%
  arrange(desc(Cancer_Type)) %>%
  mutate(ypos = cumsum(prop) - 0.5 * prop)

# Define cancer type color palette (consistent across all figures)
cancer_color <- c(
  "Nervous system tumor"    = "#1F78B4",
  "Soft tissue tumor"       = "#33A02C",
  "Gastrointestinal tumor"  = "#E31A1C",
  "Gynecologic tumor"       = "#FB9A99",
  "Neuroendocrine tumor"    = "#FF7F00",
  "Head and neck tumor"     = "#6A3D9A",
  "Mesothelioma"            = "#B15928",
  "Others"                  = "#CAB2D6",
  "Bone tumor"              = "#A6CEE3"
)

# Plot 1: Pie chart for cancer type proportion
p1 <- ggplot(pie_data, aes(x = "", y = prop, fill = Cancer_Type)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  coord_polar("y", start = 0, clip = "off") +
  geom_segment(aes(x = 1.5, y = ypos, xend = xend, yend = ypos)) +
  geom_label(aes(y = ypos, x = xend, label = paste(Cancer_Type, scales::percent(prop))),
             size = 4, color = "white", show.legend = F) +
  annotate(geom = "text", x = -0.5, y = 0,
           label = paste0("n = ", sum(pie_data$n)),
           size = 5, color = "black", fontface = "bold") +
  scale_fill_manual(values = cancer_color) +
  theme_void() +
  theme(legend.position = "none", panel.background = element_blank())

# Prepare bubble & bar plot data
df_bubble <- clin_info %>% count(Cancer_Type, Tumor_Subtype)
mut_counts <- clin_info %>% count(Cancer_Type) %>% rename(mut_count = n)
mut_summary <- full_join(total_counts, mut_counts, by = "Cancer_Type") %>%
  replace_na(list(total = 0, mut_count = 0)) %>%
  mutate(mutation_rate = ifelse(total > 0, mut_count / total * 100, 0)) %>%
  mutate(Cancer_Type = factor(Cancer_Type, levels = rev(cancer_level)))
df_bubble$Cancer_Type <- factor(df_bubble$Cancer_Type, levels = rev(cancer_level))

# Plot 2: Bar plot for mutation rate
p_bar <- ggplot(mut_summary, aes(x = mutation_rate, y = Cancer_Type, fill = Cancer_Type)) +
  geom_col() +
  scale_fill_paletteer_d(palette = "ggsci::default_igv", guide = "none") +
  geom_text(aes(label = sprintf("%.2f%%", mutation_rate), x = mutation_rate),
            hjust = -0.1, size = 3.5, color = "black") +
  theme_prism(border = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Mutation rate (%)", y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 11))

# Plot 3: Bubble plot for tumor subtype distribution
p_bubble <- ggplot(data = df_bubble, aes(x = Tumor_Subtype, y = Cancer_Type)) +
  geom_point(aes(size = n, color = Tumor_Subtype), show.legend = TRUE) +
  scale_color_paletteer_d(palette = "ggsci::default_igv", guide = "none") +
  scale_size_continuous(range = c(3.6, 8)) +
  theme_prism(border = TRUE) +
  labs(legend.title = "Sample number") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank(),
        legend.position = "top",
        plot.margin = margin(5, 5, 5, 0)) +
  guides(size = guide_legend(title = "Sample number"))

# Combine Figure 1 & save
final_plot1 <- p_bubble + p_bar + plot_layout(widths = c(9, 1.5))
ggsave("figure1_cancer_distribution.tiff", final_plot1,
       width = 18, height = 10, dpi = 300, units = "in", compression = "lzw")

# ==============================================================================
# Figure 2: Mutation Landscape (Heatmap + Bar + Tile + OR Plot)
# Core: Gene mutation type/rate, cancer-gene mutation distribution, OR analysis
# Input: figure2.rds (contains mat, clin_info, mut_data)
# Output: figure2_mutation_landscape.tiff, figure2_cancer_gene_or.tiff
# ==============================================================================
load("figure2.rdata") # Load mat, clin_info, mut_data

# Prepare gene mutation frequency data
gene_freq <- mut_data %>%
  select(-Tumor_Sample_Barcode) %>%
  summarise(across(everything(), ~ mean(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Mutation_Frequency") %>%
  arrange(desc(Mutation_Frequency))
genes_ordered <- gene_freq$Gene[gene_freq$Gene %in% rownames(mat)]
mat <- mat[genes_ordered, , drop = FALSE]

# Define mutation type color palette
color_mapping <- c(
  "Frameshift"    = "#4DAF4A",
  "Missense"      = "#a6bddb",
  "Nonsense"      = "#fdc086",
  "Other"         = "#984EA3",
  "Promoter"      = "#377EB8",
  "Splice_Site"   = "#F781BF",
  "Silent"        = "#FF7F00",
  "Span"          = "#00CED1",
  "Multi_Hit"     = "#E41A1C"
)

# Prepare mutation rate/type annotation for heatmap
mutation_pct <- round(gene_freq$Mutation_Frequency * 100, 2)
mutation_type_list <- lapply(rownames(mat), function(gene) {
  table(factor(mat[gene, ], levels = names(color_mapping)))
})
mutation_type_mat <- do.call(rbind, mutation_type_list)
rownames(mutation_type_mat) <- rownames(mat)

# Right annotation: Mutation rate text + mutation type barplot
text_anno <- row_anno_text(paste0(mutation_pct, "%"), just = "left", location = 0, gp = gpar(fontsize = 10))
bar_anno <- row_anno_barplot(mutation_type_mat,
                             gp = gpar(fill = color_mapping, col = NA),
                             border = FALSE, bar_width = 0.8, axis = T,
                             width = unit(2, "cm"),
                             axis_param = list(side = "top", gp = gpar(fontsize = 10)))
right_anno <-  rowAnnotation(
  "Mutation Rate" = text_anno,
  "Mutation Type" = bar_anno,
  annotation_name_side = "bottom",
  show_annotation_name = c(FALSE, FALSE),
  width = unit(3.2, "cm")
)

# Top annotation: Sample-level mutation type distribution
mutation_type_col_list <- lapply(colnames(mat), function(sample) {
  table(factor(mat[, sample], levels = names(color_mapping)))
})
mutation_type_col_mat <- do.call(cbind, mutation_type_col_list)
colnames(mutation_type_col_mat) <- colnames(mat)
top_anno <- columnAnnotation(
  "Mutation Type Distribution" = anno_barplot(t(mutation_type_col_mat),
                                              gp = gpar(fill = color_mapping, col = NA),
                                              border = FALSE, bar_width = 0.8, axis = TRUE,
                                              axis_param = list(side = "left", gp = gpar(fontsize = 8))),
  annotation_name_gp = gpar(fontsize = 9),
  annotation_name_side = "left",
  show_annotation_name = FALSE,
  height = unit(1, "cm")
)

# Column annotation: Clinical information
anno_cols <- c("Cancer_Type","MSI_Status","TMB_Status","Metastasis", "MutationType")
group_info <- clin_info[clin_info$Patient_ID %in% colnames(mat), ]
group_info <- group_info[match(colnames(mat), group_info$Patient_ID), ]
anno_df <- group_info %>% select(any_of(anno_cols)) %>% as.data.frame()
rownames(anno_df) <- group_info$Patient_ID

# Clinical annotation color palette
anno_colors <- list(
  Cancer_Type = cancer_color,
  MSI_Status = c("MSS" = "#4DAF4A", "MSI-H" = "#E41A1C", "NA" = "#F0F0F0"),
  TMB_Status = c("TMB-L" = "#377EB8", "TMB-H" = "#984EA3", "NA" = "#F0F0F0"),
  Metastasis = c("No" = "#FF7F00", "Yes" = "#A65628", "NA" = "#F0F0F0"),
  MutationType = c("P" = "#F781BF", "LP" = "#66C2A5", "NA" = "#F0F0F0")
)
column_anno <- HeatmapAnnotation(df = anno_df, col = anno_colors, na_col = "grey90", show_annotation_name = T)

# Plot mutation heatmap (p1)
factor <- factor(anno_df$Cancer_Type[match(colnames(mat), rownames(anno_df))], levels = cancer_level)
mat <- as.matrix(mat)
p_heatmap <- Heatmap(
  mat, name = "Mutation Type", col = color_mapping, na_col = "white",
  row_title = "Gene", column_title = "Patients",
  show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 14),
  show_column_names = FALSE, cluster_rows = FALSE,
  right_annotation = right_anno, bottom_annotation = column_anno,
  column_split = factor, top_annotation = top_anno
)

# Save mutation heatmap
tiff("figure2_mutation_heatmap.tiff", width = 24, height = 16, units = "in", res = 300, compression = "lzw")
draw(p_heatmap)
dev.off()

# Prepare cancer-gene mutation distribution data
clin_info_temp <- clin_info %>%
  separate_rows(Germ_Mut_Genes, sep = ",") %>%
  filter(Germ_Mut_Genes != "") %>%
  distinct(Patient_ID, Cancer_Type, Germ_Mut_Genes, .keep_all = TRUE)

# Plot 2: Gene mutation count barplot (p2)
bar_data <- clin_info_temp %>% count(Germ_Mut_Genes, MutationType, name = "total_mut_count") %>% arrange(desc(total_mut_count))
total_label <- bar_data %>% group_by(Germ_Mut_Genes) %>% summarise(total = sum(total_mut_count), .groups = "drop") %>% arrange(desc(total))
gene_order <- unique(total_label$Germ_Mut_Genes)
mutation_colors <- c("P" = "#CC0010", "LP" = "#5293C1")

p2 <- ggplot(bar_data, aes(x = factor(Germ_Mut_Genes, levels = gene_order),
                           y = total_mut_count, fill = MutationType)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = mutation_colors) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05)), labels = number_format(accuracy = 1)) +
  geom_text(data = total_label, aes(x = Germ_Mut_Genes, y = total + 1, label = total),
            inherit.aes = F, color = "black", size = 3, fontface = "bold") +
  labs(x = NULL, y = "Patients", fill = "P/LP") +
  theme_minimal() +
  theme(legend.position = "top", panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.line.y = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"))

# Plot 3: Cancer-gene mutation proportion tile plot (p3)
heatmap_data <- clin_info_temp %>%
  group_by(Cancer_Type) %>%
  mutate(total_in_category = n_distinct(Patient_ID)) %>%
  count(Germ_Mut_Genes,Cancer_Type,total_in_category,name = "mut_count") %>%
  mutate(proportion = mut_count / n_distinct(clin_info_temp$Patient_ID),
         category_label = paste0(Cancer_Type, " (n=", total_in_category, ")" )) %>%
  mutate(category_label = factor(category_label, levels =  c(
    "Nervous system tumor (n=124)", "Soft tissue tumor (n=78)",
    "Gastrointestinal tumor (n=43)", "Gynecologic tumor (n=21)",
    "Neuroendocrine tumor (n=16)", "Head and neck tumor (n=13)",
    "Mesothelioma (n=12)", "Bone tumor (n=2)", "Others (n=4)"
  )))

p3 <- ggplot(heatmap_data, aes(y = category_label,
                               x = factor(Germ_Mut_Genes, levels = gene_order),
                               fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(proportion > 0, mut_count, "")), color = "black", size = 3) +
  scale_fill_gradient(low = "#F4B3C2", high = "#D53E4F", labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_discrete(limits = rev) +
  labs(x = "", y = "Cancer Type", fill = "Proportion") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"))

# Plot 4: Cancer-gene OR analysis tile plot (p4)
mutation_data <- clin_info_temp %>%
  select(Patient_ID,Cancer_Type,Germ_Mut_Genes,MutationType) %>%
  mutate(MutationType = case_when(MutationType %in% c("P", "LP") ~"1", TRUE ~ NA_character_),
         MutationType = as.numeric(MutationType))

all_combinations <- expand_grid(Patient_ID = unique(mutation_data$Patient_ID),
                                Germ_Mut_Genes = unique(mutation_data$Germ_Mut_Genes)) %>%
  left_join(clin_info %>% select(Patient_ID, Cancer_Type), by = "Patient_ID") %>%
  distinct()

mutation_full <- all_combinations %>%
  mutate(Has_Mutation = if_else(paste(Patient_ID, Germ_Mut_Genes) %in%
                                  paste(mutation_data$Patient_ID, mutation_data$Germ_Mut_Genes), 1, 0))

mutation_table <- mutation_full %>%
  group_by(Cancer_Type, Germ_Mut_Genes) %>%
  summarise(Mutated = sum(Has_Mutation), Not_Mutated = n() - sum(Has_Mutation), .groups = 'drop')

# Function: Calculate Odds Ratio (OR) and 95% CI
calculate_or <- function(data, target_cancer, gene) {
  target_data <- data %>% filter(Cancer_Type == target_cancer, Germ_Mut_Genes == gene)
  other_data <- data %>% filter(Cancer_Type != target_cancer, Germ_Mut_Genes == gene) %>%
    summarise(Mutated = sum(Mutated), Not_Mutated = sum(Not_Mutated))
  
  a <- target_data$Mutated; b <- target_data$Not_Mutated
  c <- other_data$Mutated; d <- other_data$Not_Mutated
  
  or <- (a * d) / (b * c)
  log_or_se <- sqrt(1/a + 1/b + 1/c + 1/d)
  log_or <- log(or)
  ci_lower <- exp(log_or - 1.96 * log_or_se)
  ci_upper <- exp(log_or + 1.96 * log_or_se)
  
  tibble(Cancer_Type = target_cancer, Germ_Mut_Genes = gene,
         OR = or, CI_lower = ci_lower, CI_upper = ci_upper,
         p_value = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value)
}

# Calculate OR for all cancer-gene pairs
cancer_types <- unique(mutation_full$Cancer_Type)
genes <- unique(mutation_full$Germ_Mut_Genes)
or_results <- map_df(cancer_types, function(cancer) {
  map_df(genes, function(gene) {
    calculate_or(mutation_table, cancer, gene)
  })
}) %>% mutate(p_adj = p.adjust(p_value, method = "fdr"))

# Prepare OR plot data
plot_data <- or_results %>%
  mutate(OR = ifelse(is.infinite(OR) | OR == 0, NA_real_, OR),
         CI_lower = ifelse(is.na(OR), NA_real_, CI_lower),
         CI_upper = ifelse(is.na(OR), NA_real_, CI_upper)) %>%
  mutate(OR_plot = case_when(OR > 10 ~ 10, OR < 0.1 ~ 0.1, TRUE ~ OR),
         significance = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ ""
         )) %>%
  mutate(Germ_Mut_Genes = factor(Germ_Mut_Genes, levels = gene_order),
         Cancer_Type = factor(Cancer_Type, levels = cancer_level))

# Plot OR heatmap
p4 <- ggplot(plot_data, aes(x = Germ_Mut_Genes, y = Cancer_Type)) +
  geom_tile(aes(fill = log2(OR_plot)), color = "white") +
  geom_text(aes(label = significance), color = "black", size = 4, na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
  scale_y_discrete(limits = rev) +
  labs(x = "Gene", y = "Cancer Type", fill = "log2(OR)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 12, face = "bold", colour = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 14, face = "bold", colour = "black"))

# Combine cancer-gene plots & save
combined_plot2 <- p2 / p3 / p4 + plot_layout(heights = c(1, 1, 1))
ggsave("figure2_cancer_gene_or.tiff", combined_plot2,
       width = 20, height = 18, dpi = 300, units = "in", compression = "lzw")

# ==============================================================================
# Figure 3: KEGG Enrichment Analysis (Dot + Sankey Plot)
# Core: KEGG enrichment of germline mutated genes with dot-sankey combined plot
# Input: gene_info.xlsx (contains germline mutated genes in "Gene" column)
# Output: figure3_kegg_enrichment.tiff
# ==============================================================================
# Load germline genes and convert ID (SYMBOL → ENTREZID)
germ_genes <- read_xlsx("gene_info.xlsx") %>%
  select(Gene) %>%
  drop_na() %>%
  unique() %>%
  pull(Gene)

s2e <- bitr(germ_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG enrichment analysis
kegg <- enrichKEGG(
  gene = s2e$ENTREZID,
  keyType = "ncbi-geneid",
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.1
)

# Function: Dot plot + Sankey plot for enrichment result
dot_sankey <- function(data,
                       top = 10,
                       sankey.text.size = 3,
                       sankey.x.axis.text.position = 1,
                       sankey.x.axis.text.size = 10,
                       dot.color = "RdYlBu",
                       dot.position = 8,
                       dot_label = "Description",
                       dot.text.size = 10,
                       dot.scale = 0.9,
                       dot.x = 0.45,
                       dot.y = 0.17,
                       dot.width = 0.5,
                       dot.height = 0.63,
                       ...) {
  # Data format conversion
  if (isS4(data)) {
    data_df <- data@result %>% drop_na()
  } else if (is.data.frame(data)) {
    data_df <- data %>% drop_na()
  } else {
    stop("Error: Data format must be S4 object (clusterProfiler) or data frame!")
  }
  
  # Filter top pathways and calculate rich factor
  data_top <- data_df %>%
    slice(1:top) %>%
    mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
           richFactor = round(richFactor, 2)) %>%
    arrange(desc(richFactor)) %>%
    select(ID, Description, p.adjust, Count, richFactor, geneID)
  
  # Dot plot for enrichment result
  p1 <- ggplot(data_top, aes(x = richFactor, y = reorder(.data[[dot_label]], richFactor))) +
    geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    scale_colour_distiller(palette = dot.color, direction = -1, name = "-log10(FDR)") +
    labs(x = "Rich Factor", y = "") +
    theme(axis.text.x = element_text(size = dot.text.size, colour = "black"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = dot.text.size + 2, colour = "black"),
          legend.title = element_text(size = dot.text.size - 2),
          legend.text = element_text(size = dot.text.size - 3),
          legend.position = "bottom", legend.box = "horizontal")
  
  # Prepare sankey plot data
  data_sankey <- data_top %>%
    select(geneID, all_of(dot_label)) %>%
    separate_rows(geneID, sep = "/") %>%
    distinct()
  
  df_sankey <- make_long(data_sankey, colnames(data_sankey))
  df_sankey$node <- factor(df_sankey$node, levels = c(unique(data_sankey[[dot_label]]),
                                                      unique(data_sankey$geneID)))
  
  # Sankey plot for gene-pathway connection
  node_cols <- ggsci::pal_npg()(length(unique(df_sankey$node)))
  p2 <- ggplot(df_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                              fill = node, label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "black", linewidth = 0.1, show.legend = FALSE) +
    geom_sankey_text(size = sankey.text.size, color = "black", hjust = 1,
                     position = position_nudge(x = -0.03)) +
    theme_sankey(base_family = "sans") +
    scale_fill_manual(values = node_cols) +
    scale_x_discrete(labels = c("Pathway", "Gene")) +
    theme(axis.text.x = element_text(size = sankey.x.axis.text.size, colour = "black",
                                     hjust = sankey.x.axis.text.position, vjust = 10),
          axis.title = element_blank(),
          plot.margin = unit(c(0, dot.position, 0, 0), "cm"))
  
  # Combine plots
  p4 <- ggdraw() +
    draw_plot(p2) +
    draw_plot(p1, scale = dot.scale, x = dot.x, y = dot.y,
              width = dot.width, height = dot.height)
  
  return(p4)
}

# Plot KEGG enrichment dot-sankey plot & save
figure3 <- dot_sankey(data = kegg, top = 8)
ggsave("figure3_kegg_enrichment.tiff", figure3,
       width = 22, height = 12, dpi = 300, units = "in", compression = "lzw")

# ==============================================================================
# Figure 4: Second Hit Mutation Analysis (Stacked Bar Plot)
# Core: Somatic second hit mutation distribution by gene/cancer type
# Input: figure4.rdata (contains somic_data, germ_data)
# Output: figure4_second_hit_mutation.tiff
# ==============================================================================
load("figure4.rdata") # Load somic_data, germ_data

# Plot 1: Somatic mutation type distribution by gene
combo_stats <- somic_data %>%
  count(Hugo_Symbol, Variant_Classification, name = "Count") %>%
  group_by(Hugo_Symbol) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  arrange(Hugo_Symbol, -Count)

p_secondhit_type <- ggplot(combo_stats, aes(x = Hugo_Symbol, y = Percentage, fill = Variant_Classification)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            position = position_stack(vjust = 0.5), size = 2.5, color = "black", fontface = "bold") +
  labs(title = "Distribution of Second Hit Mutation",
       x = "Gene", y = "Percentage (%)", fill = "Somatic Mutation Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(), legend.background = element_blank()) +
  scale_fill_brewer(palette = "Set1")

# Prepare cancer-specific second hit data
second_hit_data <- germ_data %>% filter(Second_Hit == "Yes")
gene_hit_by_cancer <- second_hit_data %>%
  distinct(Patient_ID, Germ_Mut_Genes, Cancer_Type) %>%
  group_by(Cancer_Type, Germ_Mut_Genes) %>%
  summarise(sample_count = n(), .groups = "drop") %>%
  arrange(Cancer_Type, desc(sample_count))

# Calculate total samples per cancer type
total_samples_by_cancer <- germ_data %>%
  distinct(Patient_ID, Cancer_Type) %>%
  group_by(Cancer_Type) %>%
  summarise(total_count = n())

# Merge data and calculate percentage
gene_hit_by_cancer <- gene_hit_by_cancer %>%
  left_join(total_samples_by_cancer, by = "Cancer_Type") %>%
  mutate(percentage = sample_count / total_count * 100)

# Total second hit percentage per cancer type
total_hits_by_cancer <- gene_hit_by_cancer %>%
  group_by(Cancer_Type) %>%
  summarise(total_hits = sum(sample_count)) %>%
  left_join(total_samples_by_cancer, by = "Cancer_Type") %>%
  mutate(percentage = total_hits / total_count * 100)

# Define color palette for second hit plot
my_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#6B4E91", "#A6CE39", "#FFB32F" , "#999999"
)

# Plot 2: Second hit mutation by cancer type (stacked bar)
p_secondhit_cancer <- ggplot(gene_hit_by_cancer, aes(x = Cancer_Type, y = sample_count, fill = Germ_Mut_Genes)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sample_count), position = position_stack(vjust = 0.5), size = 3) +
  geom_text(data = total_hits_by_cancer, aes(x = Cancer_Type, y = total_hits, fill = NULL,
                                             label = sprintf("%.1f%%", percentage)),
            vjust = -0.5, size = 4, fontface = "bold") +
  theme_minimal() +
  labs(title = "Second Hit across Cancer Types",
       x = "Cancer Type", y = "Sample Count with Second Hit", fill = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = my_colors)

# Combine second hit plots & save
combined_plot4 <- p_secondhit_type / p_secondhit_cancer + plot_layout(heights = c(1, 1))
ggsave("figure4_second_hit_mutation.tiff", combined_plot4,
       width = 20, height = 14, dpi = 300, units = "in", compression = "lzw")

# ==============================================================================
# Session Information (for reproducibility)
# ==============================================================================
sessionInfo()