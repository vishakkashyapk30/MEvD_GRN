#!/usr/bin/env Rscript
# Split the Domcke et al. fetal-lung Seurat object (GEO GSM4508936,
# GSM4508936_lung_filtered.seurat.RDS.gz) into the per-cell-type layout the
# scMultiomeGRN code expects ("lung_sep_data"):
#   <out>/<Cell_type>/atac/{matrix.mtx,barcodes.tsv,peaks.tsv}    <- `peaks` assay
#   <out>/<Cell_type>/scrna/{matrix.mtx,barcodes.tsv,genes.tsv,var_features.tsv}  <- `RNA` assay
# NB the `RNA` assay of this object is a gene-body(+2 kb) ACCESSIBILITY count
# matrix (GEO sample description), which is what the paper used as "scRNA-seq".
# Cell-type directory names = metadata cell type with spaces -> "_" (as the
# official build_graph.py does).
#
# Needs only base R + Matrix. If SeuratObject is not installed, S4 slots are
# read with attr() (readRDS keeps them even without the class definition).
#
# Usage: Rscript scripts/15_scmgrn_split_seurat.R <rds(.gz)> <out_dir> [celltype_col] [min_cells]
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: <rds> <out_dir> [celltype_col] [min_cells]")
rds <- args[1]; out <- args[2]
ct_col <- if (length(args) >= 3 && args[3] != "auto") args[3] else NA
min_cells <- if (length(args) >= 4) as.integer(args[4]) else 1L
have_so <- requireNamespace("SeuratObject", quietly = TRUE)
if (have_so) suppressPackageStartupMessages(library(SeuratObject))
slot_of <- function(o, s) { v <- tryCatch(methods::slot(o, s), error = function(e) NULL)
  if (is.null(v)) v <- attr(o, s, exact = TRUE); v }

message("reading ", rds, " (SeuratObject available: ", have_so, ")")
con <- if (grepl("\\.gz$", rds)) gzfile(rds, "rb") else file(rds, "rb")
obj <- readRDS(con); close(con)
assays <- slot_of(obj, "assays")
message("assays: ", paste(names(assays), collapse = ", "))
meta <- slot_of(obj, "meta.data")
message("metadata columns: ", paste(colnames(meta), collapse = ", "))
if (is.na(ct_col)) {
  cand <- c("cell_type", "celltype", "Main_cluster_name", "cell_type_name", "CellType")
  ct_col <- cand[cand %in% colnames(meta)][1]
  if (is.na(ct_col)) stop("no cell-type column found; pass it as the 3rd argument")
}
message("cell-type column: ", ct_col)
peak_name <- intersect(c("peaks", "ATAC", "peak"), names(assays))[1]
rna_name <- intersect(c("RNA", "GeneActivity", "rna"), names(assays))[1]
P <- slot_of(assays[[peak_name]], "counts"); R <- slot_of(assays[[rna_name]], "counts")
vf <- slot_of(assays[[rna_name]], "var.features")
message(sprintf("peaks: %d x %d, RNA: %d x %d, RNA var.features: %d",
                nrow(P), ncol(P), nrow(R), ncol(R), length(vf)))
cts <- as.character(meta[[ct_col]])
tab <- sort(table(cts), decreasing = TRUE); print(tab)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(cell_type = names(tab), n = as.integer(tab)),
            file.path(out, "cell_type_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cells <- rownames(meta)
for (ct in names(tab)) {
  if (tab[[ct]] < min_cells) next
  d <- file.path(out, gsub(" ", "_", ct))
  keep <- cells[cts == ct]
  for (mod in c("atac", "scrna")) {
    M <- if (mod == "atac") P else R
    M <- M[, intersect(keep, colnames(M)), drop = FALSE]
    md <- file.path(d, mod); dir.create(md, recursive = TRUE, showWarnings = FALSE)
    writeMM(as(M, "CsparseMatrix"), file.path(md, "matrix.mtx"))
    writeLines(colnames(M), file.path(md, "barcodes.tsv"))
    writeLines(rownames(M), file.path(md, if (mod == "atac") "peaks.tsv" else "genes.tsv"))
    if (mod == "scrna" && length(vf) > 0) writeLines(vf, file.path(md, "var_features.tsv"))
  }
  message(sprintf("wrote %s (%d cells)", d, length(keep)))
}
