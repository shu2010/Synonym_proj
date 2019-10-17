#!/usr/bin/Rscript

##this script captures prioritises samples for PheWAS case-control split


phen_mat <- readRDS("~/Syn_collab/phen_mat.rds")
IBD_dist1_rel1 <- readRDS("~/Syn_collab/IBD_dist1_rel1.rds")

#function to prioritise the sample with phenotypic info.

pick_pinf <- function(ibd_dat){
  p1 <- ibd_dat[1]
  p2 <- ibd_dat[2]
  pi_hat <- ibd_dat[3]
  phen <- colnames(phen_mat)[1:5]
  c1 <- phen_mat[rownames(phen_mat) %in% c(as.character(p1), as.character(p2)), colnames(phen_mat) %in% phen]
  if(is.null(dim(c1))){
    print("sample pair not found")
    print(p1)
    print(p2)
    sel_samp <- NULL
  }
  else{
    # print("sample pair found")
    ##presence of one phenotype
    
    pheno_c1 <- c1[,which(colSums(c1) == 1)]
    if(is.table(pheno_c1)){
      pheno_sel_samp <- apply(pheno_c1, 2, function(x)names(which(x > 0)))
    }
    else {
      pheno_sel_samp <-  names(pheno_c1[which(pheno_c1 > 0)])
      names(pheno_sel_samp) <- colnames(c1)[which(colSums(c1) == 1)]
    }
    
    ##presence of none or both phenotypes
    
    pheno_c2 <- c1[,which(colSums(c1) != 1)]
    if(is.table(pheno_c2)){
      pheno_sel_samp2 <- apply(pheno_c2, 2, function(x)names(x)[sample(2, 1)])
    }
    else if(!is.table(pheno_c2)){
      # print(p1)
      # print(p2)
      pheno_sel_samp2 <-  names(pheno_c2[which(pheno_c2 != 1)])
      names(pheno_sel_samp2) <- colnames(c1)[which(colSums(c1) != 1)]
      
    }
    #pheno_sel_samp2 <- colnames(c1[,which(colSums(c1) != 1)])
    ##add criteria (Age/sex) for selection later
    
    sel_samp <- c(pheno_sel_samp, pheno_sel_samp2)
    sel_samp <- sel_samp[match(phen, names(sel_samp))]
    sel_samp <- c(p1, p2, pi_hat, sel_samp)
    #  sel_samp <- ifelse(as.character(p1) %in% sel_samp, 1, 2) ##reduce memory footprint
    # sel_samp <- sel_samp[match(phen, sel_samp)]
  }
  return(sel_samp)
}

##Parallelised version

library(parallel)
cl <- makeCluster(25)
clusterExport(cl, "phen_mat")
system.time(para_pheno <- parApply(cl, IBD_dist1_rel1, 1, pick_pinf))
stopCluster(cl)
saveRDS(para_pheno, file = "~/Syn_collab/para_IBD_phen_test.rds", compress = T)
