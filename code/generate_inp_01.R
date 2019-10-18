#!/usr/bin/Rscript

##read directly from matrix data

phen_dat <- read.delim("~/Syn_collab/test.08302019.xls.study.IDs.tsv", sep = "\t",
                       header = T)
phen_dat1 <- phen_dat[,c(1,3)]
phen_dat1[[2]] <- factor(phen_dat1[[2]], levels = levels(phen_dat1$code))

##phenotype x pat_id
phen_mat <- table(phen_dat1)
#head(colnames(phen_mat)[which(phen_mat[1,] > 0)]) (works)
##dim(phen_mat) = 92279 x 22790
#phen_mat <- readRDS("~/Syn_collab/phen_mat.rds")
##phenotypes reported gt 50 (filter)
##remove phenotypes with no sample hits or have greater
##than 50 samples in them; some IBD pairs are not present in the phen_samp mat
phen_mat <- phen_mat[,-c(which(colSums(phen_mat) <= 50))]
##plink IBD output
IBD_dist <- read.delim("~/Syn_collab/test.genome", sep = "\t", header = T)
IBD_dist1 <- IBD_dist[,c(1,3,10)]
##separate related(focus) and unrelated cases
IBD_dist1_rel <- IBD_dist1[IBD_dist1$PI_HAT >= 0.125,]
IBD_dist1_unrel <- IBD_dist1[IBD_dist1$PI_HAT < 0.125,]

##remove samples phenotype matrix that donot have IBD_dist
IBD_dist_relsamp <- unique(c(IBD_dist1_rel$FID1, IBD_dist1_rel$FID2)) ##77690 patients in total
phen_mat <- phen_mat[rownames(phen_mat) %in% as.character(IBD_dist_relsamp),]

##Adjust IBD distance dataframe
IBD_dist1_rel1 <- IBD_dist1_rel[IBD_dist1_rel$FID1 %in% rownames(phen_mat) & IBD_dist1_rel$FID2 %in% rownames(phen_mat),]

saveRDS(phen_mat, file = "~/Syn_collab/phen_mat.rds", compress = T)
saveRDS(IBD_dist1_rel1, file = "~/Syn_collab/IBD_dist1_rel1.rds", compress = T)
