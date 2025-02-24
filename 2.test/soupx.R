# soupx.R 250117
# /opt/conda/bin/R
library(optparse)
library(DropletUtils)
library(SoupX)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

files_txt_path <- args[1]
samples_txt_path <- args[2]
minCG <- as.numeric(args[3])
tfidfMin <- as.numeric(args[4])
highestrho <- as.numeric(args[5])
before_method <- args[6]

input_files_0 <- readLines(files_txt_path, warn = FALSE)
input_files_0 <- unlist(strsplit(input_files_0, ',', fixed = TRUE))
input_files_0 <- input_files_0[!is.na(input_files_0) & input_files_0 != '']
print(input_files_0)

folder_name_list <- readLines(samples_txt_path, warn = FALSE)
folder_name_list <- unlist(strsplit(folder_name_list, ',', fixed = TRUE))
folder_name_list <- folder_name_list[!is.na(folder_name_list) & folder_name_list != '']
print(folder_name_list)


if (before_method == 'scRNA-seq_v3') {
    raw_matrix_list <- sapply(seq_along(input_files_0), function(i) {
        file.path(input_files_0[i], folder_name_list[i], '02.count', 'raw_matrix')
    })
    filter_matrix_list <- sapply(seq_along(input_files_0), function(i) {
        file.path(input_files_0[i], folder_name_list[i], '02.count', 'filter_matrix')
    })
} else if (before_method == 'scRNA-seq_v3.1.5') {
    raw_matrix_list <- sapply(seq_along(input_files_0), function(i) {
        file.path(input_files_0[i], '02.cDNAAnno', 'RawMatrix')
    })
    filter_matrix_list <- sapply(seq_along(input_files_0), function(i) {
        file.path(input_files_0[i], '04.Matrix', 'FilterMatrix')
    })
} else {
    stop("Unsupported before_method: ", before_method)
}

print(raw_matrix_list)
print(filter_matrix_list)

run_soupx <- function(rawMatrix, filterMatrix, outdir, minCG, tfidfMin, highestrho) {
    options(future.globals.maxSize = 100000 * 1024^3)
    tod <- Read10X(rawMatrix, gene.column=1)
    toc <- Read10X(filterMatrix, gene.column=1)
    
    all <- CreateSeuratObject(toc)
    all <- subset(all, subset = nFeature_RNA > minCG)
    all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
    all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(all)
    all <- ScaleData(all, features = all.genes)
    all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = FALSE)
    all <- FindNeighbors(all, dims = 1:30)
    all <- FindClusters(all, resolution = 0.5)
    all <- RunUMAP(all, dims = 1:30)
    
    matx <- all@meta.data
    toc <- GetAssayData(object = all, layer = "counts", assay = "RNA")
    
    raw <- CreateSeuratObject(tod)
    tod <- GetAssayData(object = raw, layer = "counts", assay = "RNA")
    tod <- tod[rownames(all),]
    
    sc <- SoupChannel(tod, toc)
    sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    pdfrho <- paste0(outdir, "_rho.pdf")
    pdf(pdfrho, width=9)
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE)
    dev.off()
    
    rho_value <- unique(sc$metaData$rho)
    rho_adjust <- ifelse(rho_value < highestrho, 'good', 'bad')
    
    out <- adjustCounts(sc)
    DropletUtils::write10xCounts(outdir, out, version="3")
    
    return(list(rho_value = rho_value, rho_adjust = rho_adjust))
}

rho_value_list <- list()
rho_adjust_list <- list()

for (i in 1:length(folder_name_list)) {
    result <- run_soupx(raw_matrix_list[i], filter_matrix_list[i], folder_name_list[i], minCG=minCG, tfidfMin=tfidfMin, highestrho=highestrho)
    rho_value_list[[i]] <- result$rho_value
    rho_adjust_list[[i]] <- result$rho_adjust
}

rho_value_list
rho_adjust_list

file_conn <- file('soupx_rho.txt', open = "w")
cat("highestrho:", highestrho, "\n", file = file_conn)
for (i in seq_along(folder_name_list)) {
    cat("sample:", folder_name_list[i], "\n", file = file_conn)
    cat("rho:", rho_value_list[[i]][1], "\n", file = file_conn)  # 提取列表中的向量的第一个值
    cat("how:", rho_adjust_list[[i]][1], "\n", file = file_conn)  # 提取列表中的向量的第一个值
}
close(file_conn)