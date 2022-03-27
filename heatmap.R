library(tidyverse)
data_all <- list.files(path = "~/Desktop/scd6_data/accepted/log2fc",     
                       pattern = "*.csv", full.names = TRUE) %>%lapply(read_csv) %>%bind_cols
# (%>%) in order to perforn two steps in a single line
data_all <- data_all[,c(3,10,17,24)] # log2foldchange columns
colnames(data_all) <- c("scd6[delta]ntd_none_scd6[delta]ntd_3%-urea",
                        "scd6_wild_none_tif3[delta]ntd_3%-urea", "wild_3%-urea_scd6[delta]ntd_3%-urea", "wild_3%-urea_scd6[delta]ntd_none")
data_all <- data_all %>% column_to_rownames(var="gene") #column_to_rownames to presetve the column names 
data_all[is.na(data_all)] <- 0          # computes the gene variance 
library(ComplexHeatmap)
library(circlize)
Heatmap(data_all, name = "Expression",clustering_distance_columns=function(x) as.dist(1-cor(t(x))), clustering_method_columns="ward.D2",
        clustering_distance_rows="euclidean",clustering_method_rows="ward.D2", column_split = 6, col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        column_title="Heatmap of scd6 after 3% urea treatment",column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        width = unit(10, "cm"), height = unit(17, "cm"),column_names_rot = 45,
        column_names_gp = gpar(fontsize = 8))
# euclidean distance method has been used here for generating heatmap 