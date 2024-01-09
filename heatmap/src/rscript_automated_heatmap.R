################################################################################
### HEATMAP CODE ###############################################################
################################################################################
### Structure
### 1. input and prepare parameters
### 2. prepare library
### 3. find, prepare and transform csv file
### 4. matrix transformation
### 5. prepare row and col names
### 6. generate heatmaps
### End
################################################################################

### 1. input and prepare parameters
#params prepared by "heatmap.py"
csv_file <- "<<INPUT_CSV>>"
len_target <- as.numeric("<<INPUT_LEN_TARGET>>")
str_probe_vector <- "<<PROBE_VECTOR_FOR_R>>"
name_output_file_pdf <- "<<NAME_PDF_OUTPUT>>"

#other params
target_start <- as.numeric("<<INPUT_TARGET_START>>")
target_end <- as.numeric("<<INPUT_TARGET_END>>")
probe_name <- "<<PROBE_NAME>>"

#heatmap params
y_axis_label <- "<<INPUT_Y_AXIS_LAB>>"
y_axis_increment <- as.numeric("<<INPUT_Y_AXIS_INCREMENT>>")
heatmap_units <- "kcal/mol/nt"
colnames_heatmap <- "<<COLNAMES>>"

#heatmap titles
title_hm_each_probe_mean <- paste("Heatmap of", y_axis_label, "\nand", probe_name, "(mean of each probe)")
title_hm_matrix_mean <- paste("Heatmap of", y_axis_label, "\nand", probe_name, "(mean of all probes)")

#transforming probe_vector (string) into numeric vector
str_probe_vector <- strsplit(str_probe_vector, ", ")
numeric_probe_vector <- c()
for (i in 1:length(unlist(str_probe_vector))) {
    numeric_probe_vector <- append(
        numeric_probe_vector,
        as.numeric(unlist(str_probe_vector)[i])
    )
}

#transforming colname_heatmap into character vector
colnames_heatmap <- strsplit(colnames_heatmap, ",")
col_vector <- c()
for (i in 1:length(unlist(colnames_heatmap))) {
    col_vector <- append(col_vector, as.character(unlist(colnames_heatmap)[i]))
}


### 2. prepare library
suppressPackageStartupMessages(library(tidyverse)) # manipulate data
suppressPackageStartupMessages(library(ComplexHeatmap)) # make heatmaps
suppressPackageStartupMessages(library(viridis)) # colours
suppressPackageStartupMessages(library(RColorBrewer)) # colours


### 3. find, prepare and transform csv file
#finding file
read_csv_file_probes <- read.csv(csv_file, header = TRUE)
#todo assert this exists
print("Found CSV File")

#correcting for probe length
for (col in 1:ncol(read_csv_file_probes)) {
    read_csv_file_probes[, col] <- read_csv_file_probes[, col] /
    numeric_probe_vector[col]
}
print("Corrected for length of probe")

#calculate mean by probe
means_by_probe <- c() # create empty vector "means_by_probe"
means_whole_matrix_sum_array <- c()
for (col in 1:ncol(read_csv_file_probes)) {
    sum_col <- sum(read_csv_file_probes[, col], na.rm = TRUE)
    col_mean <- (sum_col) / ((1 * len_target) - (1 * numeric_probe_vector[col]))
    means_by_probe <- append(means_by_probe, col_mean)
    means_whole_matrix_sum_array <- append(means_whole_matrix_sum_array, sum_col)
}
print("Created means array for each probe used:")

for (i in 1:length(means_by_probe)) {
    mean_rounded <- round(means_by_probe[i], digits = 2)
    print(
        paste("    -> mean of probe", i, "is:", mean_rounded)
    )
}

#calculate mean for whole matrix
num_mean_all <- sum(means_whole_matrix_sum_array)
denom_mean_all <- (len_target * ncol(read_csv_file_probes)) - (sum(numeric_probe_vector))

mean_whole_matrix <- num_mean_all / denom_mean_all

print(paste("And -> mean of entire matrix is:", round(mean_whole_matrix, digits = 2)))

### 4. transform matrix
matrix1_probe <- as.matrix(read_csv_file_probes)
for (col in 1:ncol(matrix1_probe)) {
    col_mean <- means_by_probe[col]
    matrix1_probe[matrix1_probe[, col] > col_mean, col] <- col_mean
}

matrix2_all <- as.matrix(read_csv_file_probes)
for (col in 1:ncol(matrix2_all)) {
    matrix_mean <- mean_whole_matrix
    matrix2_all[matrix2_all[, col] > matrix_mean, col] <- matrix_mean
}
print(paste("Created 2 matrices (mean of each probe, mean of entire matrix), with below the respective means transformed to the mean"))


### 5. prepare row and col names
#mean annotation (for matrix1_probe)
legend_each_mean <- c()
for (col in 1:ncol(matrix1_probe)) {
    mean_col <- round(means_by_probe[col], digits = 2)
    label <- paste0("(µ = ", mean_col, ")")
    legend_each_mean <- c(legend_each_mean, label)
}

gene_annotate_each_mean <- HeatmapAnnotation(
    Genes = anno_block(
        gp = gpar(fill = "white"),
        labels = c(legend_each_mean),
        labels_gp = gpar(col = "black", fontsize = 8)
    )
)

#remove rownames, set target_end and paste incremental row names
rownames(matrix1_probe)[1:dim(matrix1_probe)[1]] <- " "
rownames(matrix2_all)[1:dim(matrix2_all)[1]] <- " "

if (target_start < target_end) {
    target_next_increment <- target_start
    target_end <- target_start + dim(matrix1_probe)[1] - 1
    target_end <- target_start + dim(matrix2_all)[1] - 1

    for (n in seq(
        from = 1,
        to = dim(matrix1_probe)[1],
        by = y_axis_increment
        )
    ) {
        target_next_increment <- target_next_increment + y_axis_increment
        rownames(matrix1_probe)[n] <- paste0("==", target_next_increment, "=>")
        rownames(matrix2_all)[n] <- paste0("==", target_next_increment, "=>")
    }

} else {
    target_next_increment <- target_start
    target_end <- target_start - dim(matrix1_probe)[1] + 1
    target_end <- target_start - dim(matrix2_all)[1] + 1

    for (
        n in seq(
            from = y_axis_increment,
            to = target_start,
            by = y_axis_increment
        )
    ) {
        target_next_increment <- target_next_increment - y_axis_increment
        rownames(matrix1_probe)[n] <- paste0("==", target_next_increment, "=>")
        rownames(matrix2_all)[n] <- paste0("==", target_next_increment, "=>")
    }
}

#rownames for start and end of target
rownames(matrix1_probe)[1] <- paste0("Start: ", target_start)
rownames(matrix2_all)[1] <- paste0("Start: ", target_start)
rownames(matrix1_probe)[dim(matrix1_probe)[1]] <- paste0("End: ", target_end)
rownames(matrix2_all)[dim(matrix2_all)[1]] <- paste0("End: ", target_end)

#prepare colourblind friendly palette
brown_green <- brewer.pal(7, "YlGnBu")


### 6. generate heatmap
pdf(name_output_file_pdf)

#heatmap1 using matrix with mean of each probe
Heatmap(
    matrix1_probe,
    cluster_columns = FALSE,
    column_title = title_hm_each_probe_mean,
    column_title_gp = gpar(fontsize = 16),
    column_split = rep(1:ncol(matrix1_probe)), # nolint: seq_linter.
    bottom_annotation = gene_annotate_each_mean,
    column_labels = col_vector,
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_title = y_axis_label,
    row_title_gp = gpar(fontsize = 14),
    show_row_names = TRUE,
    row_names_side = "left",
    cluster_rows = FALSE,
    name = heatmap_units,
    col = viridis(10, direction = 1),
)
print("Generated a heatmap displaying the mean of each probe")

#heatmap2 using matrix with mean of all genes
Heatmap(
    matrix2_all,
    cluster_columns = FALSE,
    column_title = title_hm_matrix_mean,
    column_title_gp = gpar(fontsize = 16),
    column_split = rep(1:ncol(matrix2_all)),
    column_labels = col_vector,
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_title = y_axis_label,
    row_title_gp = gpar(fontsize = 14),
    show_row_names = TRUE,
    row_names_side = "left",
    cluster_rows = FALSE,
    name = paste0(
        heatmap_units,
        "\n(µ= ",
        round(mean_whole_matrix, digits = 2),
        ")"
    ),
    col = viridis(10, direction = 1)
)
print("Generated a heatmap displaying the mean of entire matrix")

dev.off()
