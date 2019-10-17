##running step3
##preprocessing
p5_s50 <- readRDS("~/Syn_collab/para_IBD_phen_test.rds")
p5_s50_mat <- t(p5_s50)
class(p5_s50_mat) <- "numeric"
source("~/Syn_collab/check_rm_dup_03.R")
t2_mat <- dup_adj(p5_s50_mat) ##final output
saveRDS(t2_mat, file = "~/Syn_collab/phen_selected.rds", compress = T)