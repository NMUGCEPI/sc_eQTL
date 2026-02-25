##coloc
library(data.table)
library(tidyfst)
library(coloc)

qtl_coloc <- function(gene_id, coloc_file, sample_size) {
    eqtl_coloc <- list(
        snp = coloc_file$SNP,
        position = coloc_file$bp,
        beta = coloc_file$beta_eqtl,
        varbeta = coloc_file$varbeta_eqtl,
        type = "quant",
        N = as.numeric(sample_size),
        MAF = coloc_file$maf
    )
    gwas_coloc <- list(
        snp = coloc_file$SNP,
        position = coloc_file$bp,
        beta = coloc_file$beta_gwas,
        varbeta = coloc_file$varbeta_gwas,
        type = "cc"
    )
    result <- coloc.abf(dataset1 = eqtl_coloc, dataset2 = gwas_coloc)
    res <- t(as.data.frame(result$summary))
    rownames(res) <- gene_id
    res <- as.data.frame(res)
    return(res)
}

dat <- fread('coloc_input.xls')
gene_list <- dat %>% select_dt(gene_id) %>% distinct_dt() %>% pull_dt(gene_id)
sample_size <- 203 #sample size of eqtl data

res_all <- data.frame()
    for (j in 1:length(gene_list)) {
        sub <- dat %>%
            filter_fst(gene_id == gene_list[j])
        if (nrow(sub) > 30) {
            res <- qtl_coloc(gene_list[j], sub, sample_size)
            res_all <- rbind(res_all, res)
        }
    }

