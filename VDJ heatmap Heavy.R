# Lukas Forssander-Song
# VDJ Heatmap Script for BC_Teen and Golden_Ticket Heavy filtered
# 2025-11-05

# Load packages----------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)


# Load data as csv-------------------------------------------------------------

# BC_Teen data
BC_HAneg_CS <- read.csv("BC_teen_HAneg_class_switched.csv")
BC_HAneg_N <- read.csv("BC_teen_HAneg_naive.csv")
BC_HApos_CS <- read.csv("BC_teen_HApos_class_switched.csv")
BC_HApos_N <- read.csv("BC_teen_HApos_naive.csv")

# Golden_ticket data
Gold_HAneg_CS <- read.csv("Golden_Ticket_HAneg_class_switched.csv")
Gold_HAneg_N <- read.csv("Golden_Ticket_HAneg_naive.csv")
Gold_HApos_CS <- read.csv("Golden_Ticket_HApos_class_switched.csv")
Gold_HApos_N <- read.csv("Golden_Ticket_HApos_naive.csv")


# V heatmap--------------------------------------------------------------------

# Group dataframes by V.name, filter for Heavy chain, and calculate frequency 
# as a proportion
BC_HAneg_CS_V <- BC_HAneg_CS |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(BC_HAneg_CS_freq = Total_clones/sum(Total_clones)*100) 

BC_HAneg_N_V <- BC_HAneg_N |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(BC_HAneg_N_freq = (Total_clones/sum(Total_clones))*100) 

BC_HApos_CS_V <- BC_HApos_CS |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(BC_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

BC_HApos_N_V <- BC_HApos_N |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(BC_HApos_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_CS_V <- Gold_HAneg_CS |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(Gold_HAneg_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_N_V <- Gold_HAneg_N |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(Gold_HAneg_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_CS_V <- Gold_HApos_CS |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(Gold_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_N_V <- Gold_HApos_N |>
  select(Clones, V.name) |>
  group_by(V.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(V.name, "IGH")) |>
  mutate(Gold_HApos_N_freq = Total_clones/sum(Total_clones)*100)

# Create list of V dataframes to join
V_list <- list(BC_HAneg_CS_V, BC_HAneg_N_V, BC_HApos_CS_V, BC_HApos_N_V,
               Gold_HAneg_CS_V, Gold_HAneg_N_V, Gold_HApos_CS_V, Gold_HApos_N_V
)

# Create condensed matrix from V dataframes where similar V.names are merged 
# for pheatmap
V_matrix_condensed <- V_list |>
  purrr::reduce(full_join, by = "V.name") |>
  select(V.name, ends_with("_freq")) |>
  mutate(across(everything(), ~replace(., is.na(.), 0))) |>
  mutate(V.name = str_replace(V.name, "-.*", "")) |>
  group_by(V.name) |>
  summarise(across(everything(), sum)) |>
  column_to_rownames("V.name") |>
  rename(BC_teen_HAneg_class_switched = BC_HAneg_CS_freq,
         BC_teen_Haneg_naive = BC_HAneg_N_freq,
         BC_teen_HApos_class_switched = BC_HApos_CS_freq,
         BC_teen_HApos_naive = BC_HApos_N_freq,
         Golden_ticket_HAneg_class_switched = Gold_HAneg_CS_freq,
         Golden_ticket_HAneg_naive = Gold_HAneg_N_freq,
         Golden_ticket_HApos_class_switched = Gold_HApos_CS_freq,
         Golden_ticket_HApos_naive = Gold_HApos_N_freq) |>
  as.matrix()

# Plot V heatmap using pheatmap
pheatmap(V_matrix_condensed, # Condensed V genes with "-XX" combined
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         color = colorRampPalette(c("white", "firebrick1", "firebrick4"))(100),
         breaks = seq(0, 60, length.out = 101))

# Plot V heatmap using ComplexHeatmap

## Identify which columns belong to BC_teen or Golden_ticket
BC_cols_V <- grep("^BC_teen_", colnames(V_matrix_condensed))
GT_cols_V <- grep("^Golden_ticket_", colnames(V_matrix_condensed))

## Define the custom colour palettes (values range roughly 0-60)
col_fun_BC_V <- colorRamp2(c(0, 30, 60), c("white","firebrick2", "firebrick4"))
col_fun_GT_V <- colorRamp2(c(0, 30, 60), c("white", "skyblue2", "skyblue4"))

## Create heatmaps for each subset
ht_BC_V <- Heatmap(
  V_matrix_condensed[, BC_cols_V],
  name = "BC_teen",
  col = col_fun_BC_V,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "V gene",
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

ht_GT_V <- Heatmap(
  V_matrix_condensed[, GT_cols_V],
  name = "Golden_ticket",
  col = col_fun_GT_V,
  cluster_rows = FALSE, # keeps same row order at ht_BC
  cluster_columns = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

## Combine and draw heatmaps side-by-side with aligned rows
draw(ht_BC_V + ht_GT_V)

# D heatmap--------------------------------------------------------------------

# Group dataframes by D.name, filter for Heavy chain, and calculate frequency 
# as a proportion
BC_HAneg_CS_D <- BC_HAneg_CS |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(BC_HAneg_CS_freq = Total_clones/sum(Total_clones)*100) 

BC_HAneg_N_D <- BC_HAneg_N |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(BC_HAneg_N_freq = (Total_clones/sum(Total_clones))*100) 

BC_HApos_CS_D <- BC_HApos_CS |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(BC_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

BC_HApos_N_D <- BC_HApos_N |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(BC_HApos_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_CS_D <- Gold_HAneg_CS |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(Gold_HAneg_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_N_D <- Gold_HAneg_N |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(Gold_HAneg_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_CS_D <- Gold_HApos_CS |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(Gold_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_N_D <- Gold_HApos_N |>
  select(Clones, D.name) |>
  group_by(D.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(D.name, "IGH")) |>
  mutate(Gold_HApos_N_freq = Total_clones/sum(Total_clones)*100)

# Create list of D dataframes to join
D_list <- list(BC_HAneg_CS_D, BC_HAneg_N_D, BC_HApos_CS_D, BC_HApos_N_D,
               Gold_HAneg_CS_D, Gold_HAneg_N_D, Gold_HApos_CS_D, Gold_HApos_N_D
)

# Create matrix from D dataframes for pheatmap
D_matrix <- D_list |>
  purrr::reduce(full_join, by = "D.name") |>
  select(D.name, ends_with("_freq")) |>
  mutate(across(everything(), ~replace(., is.na(.), 0))) |>
  mutate(across(everything(), ~replace(., . == "", 0))) |>
  mutate(D.name = str_replace(D.name, "-.*", "")) |>
  mutate(D.name = str_replace(D.name, "/.*", "")) |>
  group_by(D.name) |>
  summarise(across(everything(), sum)) |>
  column_to_rownames("D.name") |>
  rename(BC_teen_HAneg_class_switched = BC_HAneg_CS_freq,
         BC_teen_Haneg_naive = BC_HAneg_N_freq,
         BC_teen_HApos_class_switched = BC_HApos_CS_freq,
         BC_teen_HApos_naive = BC_HApos_N_freq,
         Golden_ticket_HAneg_class_switched = Gold_HAneg_CS_freq,
         Golden_ticket_HAneg_naive = Gold_HAneg_N_freq,
         Golden_ticket_HApos_class_switched = Gold_HApos_CS_freq,
         Golden_ticket_HApos_naive = Gold_HApos_N_freq) |>
  as.matrix()

# Plot D pheatmap
pheatmap(D_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9,
         color = colorRampPalette(c("white", "firebrick1", "firebrick4"))(100),
         breaks = seq(0, 40, length.out = 101))

# Plot D heatmap using ComplexHeatmap

## Identify which columns belong to BC_teen or Golden_ticket
BC_cols_D <- grep("^BC_teen_", colnames(D_matrix))
GT_cols_D <- grep("^Golden_ticket_", colnames(D_matrix))

## Define the custom colour palettes (values range roughly 0-40)
col_fun_BC_D <- colorRamp2(c(0, 20, 40), c("white","firebrick2", "firebrick4"))
col_fun_GT_D <- colorRamp2(c(0, 20, 40), c("white", "skyblue2", "skyblue4"))

## Create heatmaps for each subset
ht_BC_D <- Heatmap(
  D_matrix[, BC_cols_D],
  name = "BC_teen",
  col = col_fun_BC_D,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "D gene",
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 12)
)

ht_GT_D <- Heatmap(
  D_matrix[, GT_cols_D],
  name = "Golden_ticket",
  col = col_fun_GT_D,
  cluster_rows = FALSE, # keeps same row order at ht_BC
  cluster_columns = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 12)
)

## Combine and draw heatmaps side-by-side with aligned rows
draw(ht_BC_D + ht_GT_D)


# J heatmap--------------------------------------------------------------------

# Group dataframes by J.name, filter for Heavy chain, and calculate frequency
# as a proportion
BC_HAneg_CS_J <- BC_HAneg_CS |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(BC_HAneg_CS_freq = Total_clones/sum(Total_clones)*100) 

BC_HAneg_N_J <- BC_HAneg_N |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(BC_HAneg_N_freq = (Total_clones/sum(Total_clones))*100) 

BC_HApos_CS_J <- BC_HApos_CS |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(BC_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

BC_HApos_N_J <- BC_HApos_N |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(BC_HApos_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_CS_J <- Gold_HAneg_CS |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(Gold_HAneg_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_N_J <- Gold_HAneg_N |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(Gold_HAneg_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_CS_J <- Gold_HApos_CS |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(Gold_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_N_J <- Gold_HApos_N |>
  select(Clones, J.name) |>
  group_by(J.name) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(J.name, "IGH")) |>
  mutate(Gold_HApos_N_freq = Total_clones/sum(Total_clones)*100)

# Create list of J dataframes to join
J_list <- list(BC_HAneg_CS_J, BC_HAneg_N_J, BC_HApos_CS_J, BC_HApos_N_J,
               Gold_HAneg_CS_J, Gold_HAneg_N_J, Gold_HApos_CS_J, Gold_HApos_N_J
)

# Create matrix from J dataframes for pheatmap
J_matrix <- J_list |>
  purrr::reduce(full_join, by = "J.name") |>
  select(J.name, ends_with("_freq")) |>
  mutate(across(everything(), ~replace(., is.na(.), 0))) |>
  column_to_rownames("J.name") |>
  rename(BC_teen_HAneg_class_switched = BC_HAneg_CS_freq,
         BC_teen_Haneg_naive = BC_HAneg_N_freq,
         BC_teen_HApos_class_switched = BC_HApos_CS_freq,
         BC_teen_HApos_naive = BC_HApos_N_freq,
         Golden_ticket_HAneg_class_switched = Gold_HAneg_CS_freq,
         Golden_ticket_HAneg_naive = Gold_HAneg_N_freq,
         Golden_ticket_HApos_class_switched = Gold_HApos_CS_freq,
         Golden_ticket_HApos_naive = Gold_HApos_N_freq) |>
  as.matrix()

# Plot J pheatmap
pheatmap(J_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9,
         color = colorRampPalette(c("white", "firebrick1", "firebrick4"))(100),
         breaks = seq(0, 65, length.out = 101))


# Plot J heatmap using ComplexHeatmap

## Identify which columns belong to BC_teen or Golden_ticket
BC_cols_J <- grep("^BC_teen_", colnames(J_matrix))
GT_cols_J <- grep("^Golden_ticket_", colnames(J_matrix))

## Define the custom colour palettes (values range roughly 0-60)
col_fun_BC_J <- colorRamp2(c(0, 30, 60), c("white","firebrick2", "firebrick4"))
col_fun_GT_J <- colorRamp2(c(0, 30, 60), c("white", "skyblue2", "skyblue4"))

## Create heatmaps for each subset
ht_BC_J <- Heatmap(
  J_matrix[, BC_cols_J],
  name = "BC_teen",
  col = col_fun_BC_J,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "J gene",
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

ht_GT_J <- Heatmap(
  J_matrix[, GT_cols_J],
  name = "Golden_ticket",
  col = col_fun_GT_J,
  cluster_rows = FALSE, # keeps same row order at ht_BC
  cluster_columns = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

## Combine and draw heatmaps side-by-side with aligned rows
draw(ht_BC_J + ht_GT_J)


# C Gene heatmap---------------------------------------------------------------

# Group dataframes by c_gene, filter for Heavy chain and calculate frequency 
# as a proportion
BC_HAneg_CS_c_gene <- BC_HAneg_CS |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(BC_HAneg_CS_freq = Total_clones/sum(Total_clones)*100) 

BC_HAneg_N_c_gene <- BC_HAneg_N |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(BC_HAneg_N_freq = (Total_clones/sum(Total_clones))*100) 

BC_HApos_CS_c_gene <- BC_HApos_CS |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(BC_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

BC_HApos_N_c_gene <- BC_HApos_N |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(BC_HApos_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_CS_c_gene <- Gold_HAneg_CS |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(Gold_HAneg_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HAneg_N_c_gene <- Gold_HAneg_N |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(Gold_HAneg_N_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_CS_c_gene <- Gold_HApos_CS |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(Gold_HApos_CS_freq = Total_clones/sum(Total_clones)*100)

Gold_HApos_N_c_gene <- Gold_HApos_N |>
  select(Clones, c_gene) |>
  group_by(c_gene) |>
  summarise(Total_clones = sum(Clones)) |>
  filter(startsWith(c_gene, "IGH")) |>
  mutate(Gold_HApos_N_freq = Total_clones/sum(Total_clones)*100)

# Create list of c_gene dataframes to join
c_gene_list <- list(BC_HAneg_CS_c_gene, BC_HAneg_N_c_gene, BC_HApos_CS_c_gene, 
                    BC_HApos_N_c_gene, Gold_HAneg_CS_c_gene, 
                    Gold_HAneg_N_c_gene, Gold_HApos_CS_c_gene, 
                    Gold_HApos_N_c_gene)

# Create matrix from c_gene dataframes for pheatmap
c_gene_matrix <- c_gene_list |>
  purrr::reduce(full_join, by = "c_gene") |>
  select(c_gene, ends_with("_freq")) |>
  mutate(across(everything(), ~replace(., is.na(.), 0))) |>
  column_to_rownames("c_gene") |>
  rename(BC_teen_HAneg_class_switched = BC_HAneg_CS_freq,
         BC_teen_Haneg_naive = BC_HAneg_N_freq,
         BC_teen_HApos_class_switched = BC_HApos_CS_freq,
         BC_teen_HApos_naive = BC_HApos_N_freq,
         Golden_ticket_HAneg_class_switched = Gold_HAneg_CS_freq,
         Golden_ticket_HAneg_naive = Gold_HAneg_N_freq,
         Golden_ticket_HApos_class_switched = Gold_HApos_CS_freq,
         Golden_ticket_HApos_naive = Gold_HApos_N_freq) |>
  as.matrix()

# Plot c_gene pheatmap
pheatmap(c_gene_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9,
         color = colorRampPalette(c("white", "firebrick1", "firebrick4"))(100),
         breaks = seq(0, 100, length.out = 101))

# Plot D heatmap using ComplexHeatmap

## Identify which columns belong to BC_teen or Golden_ticket
BC_cols_c_gene <- grep("^BC_teen_", colnames(c_gene_matrix))
GT_cols_c_gene <- grep("^Golden_ticket_", colnames(c_gene_matrix))

## Define the custom colour palettes (values range roughly 0-60)
col_fun_BC_c_gene <- colorRamp2(c(0, 50, 100), c("white","firebrick2", "firebrick4"))
col_fun_GT_c_gene <- colorRamp2(c(0, 50, 100), c("white", "skyblue2", "skyblue4"))

## Create heatmaps for each subset
ht_BC_c_gene <- Heatmap(
  c_gene_matrix[, BC_cols_c_gene],
  name = "BC_teen",
  col = col_fun_BC_c_gene,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "Constant region",
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

ht_GT_c_gene <- Heatmap(
  c_gene_matrix[, GT_cols_c_gene],
  name = "Golden_ticket",
  col = col_fun_GT_c_gene,
  cluster_rows = FALSE, # keeps same row order at ht_BC
  cluster_columns = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

## Combine and draw heatmaps side-by-side with aligned rows
draw(ht_BC_c_gene + ht_GT_c_gene)
