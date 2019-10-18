##test output
##QC
p5_s50 <- readRDS("~/Documents/Syn_collab/Data/para_IBD_phen_test.rds")
p5_s50_mat <- t(p5_s50)
class(p5_s50_mat) <- "numeric"
apply(p5_s50_mat[,-c(1:3)], 2, function(x)table(duplicated(x)))

#Dup removed QC
dup_rm_sel <- readRDS("~/Syn_collab/phen_selected.rds")
apply(dup_rm_sel[,-c(1:3)], 2, function(x)table(duplicated(x)))
apply(dup_rm_sel[,-c(1:3)], 2, function(x)length(unique(x)))

##check with positive control
pos_cont <- read.delim("~/Syn_collab/selected.txt", header = T, sep = "")
apply(dup_rm_sel[,-c(1:3)], 2, function(x)table(unique(x) %in% pos_cont$FID))