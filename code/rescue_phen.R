##Steps include pre-processing, rescue and duplicate removal
##read directly from matrix data

##This section is run only once to build sample x phenotype matrix
phen_dat <- read.delim("~/Syn_collab/test.08302019.xls.study.IDs.tsv", sep = "\t",
                      header = T)

phen_dat1 <- phen_dat[,c(1,3)]
phen_dat1[[2]] <- factor(phen_dat1[[2]], levels = levels(phen_dat1$code))

##phenotype x pat_id
phen_mat <- table(phen_dat1)
saveRDS(phen_mat, file = "~/Syn_collab/phen_mat.rds", compress = T)

##The saved matrix is read directly for subsequent rescue step; donot repeat the above steps
##after the matrix is generated. 
phen_mat <- readRDS("~/Syn_collab/phen_mat.rds")
##plink IBD output
IBD_dist <- read.delim("~/Syn_collab/test.genome", sep = "\t", header = T)

#function to prioritise the sample with phenotypic info.

pick_pinf <- function(ibd_dat){
  p1 <- ibd_dat[1]
  p2 <- ibd_dat[2]
  pi_hat <- ibd_dat[3]
  phen <- colnames(phen_mat)
  c1 <- phen_mat[rownames(phen_mat) %in% c(as.character(p1), as.character(p2)), colnames(phen_mat) %in% phen]
  if(is.null(dim(c1))){
    sel_samp <- NULL
  }
  else{
    
    ##presence of one phenotype
    
    pheno_c1 <- c1[,which(colSums(c1) == 1)]
    pheno_sel_samp <- apply(pheno_c1, 2, function(x)names(which(x > 0)))
    
    ##presence of none or both phenotypes; pick sample randomly
    
    pheno_c2 <- c1[,which(colSums(c1) != 1)]
    pheno_sel_samp2 <- apply(pheno_c2, 2, function(x)names(x)[sample(2, 1)])
    
    sel_samp <- c(pheno_sel_samp, pheno_sel_samp2)
    sel_samp <- sel_samp[match(phen, names(sel_samp))]
  }
  return(sel_samp)
}


IBD_dist1 <- IBD_dist[,c(1,3,10)]

test_pheno <- apply(IBD_dist1, 1, pick_pinf)
class(test_pheno) <- "numeric"
saveRDS(test_pheno, file = "~/Documents/Syn_collab/IBD_phen_df.rds", compress = T)


##Duplicate removal(under test)
# ##check duplication of samples in randomly picked cases

dup_adj <- function(tab_phe){
  print(dim(tab_phe))
  if(dim(tab_phe)[2] == 4){
    print("dimension eq 4")
    cname <- colnames(tab_phe)[4]
    print(cname)
    t_chk <- tab_phe[,-c(1:3)]
    
    t_chk[duplicated(t_chk)] <- ifelse(as.character(t_chk[duplicated(t_chk)]) %in% as.character(tab_phe$FID1)
                                       , tab_phe$FID2, tab_phe$FID1)
    if(grep("TRUE", duplicated(t_chk)) > 0){
      print("more duplicates detected")
      t_chk[duplicated(t_chk)] <- ifelse(as.character(t_chk[duplicated(t_chk)]) %in% as.character(tab_phe$FID1)
                                         , tab_phe$FID2, tab_phe$FID1)
      round1 <- cbind.data.frame(tab_phe[,c(1:3)], t_chk)
      colnames(round1)[4] <- cname
      return(round1)
    }
    else{
      # class(round1) <- "numeric"
      print("all duplicates removed")
      round1 <- cbind.data.frame(tab_phe[,c(1:3)], t_chk)
      colnames(round1)[4] <- cname
      return(round1)
    }
  }
  
  else{  
    print("dimension gt 1")
    t_chk <- tab_phe[,-c(1:3)]
    ind_list <- sapply(t_chk, function(x)as.character(x[duplicated(x)]))
    ind_list <- ind_list[lapply(ind_list,length)>0]
    t_chk <- t_chk[,colnames(t_chk) %in% names(ind_list)]
    
    if(dim(t_chk)[2] == 2){
      print("dimension eq 2")
      ind_list1 <- apply(t_chk, 2, function(x)which(grepl(paste(as.character(x[duplicated(x)]),collapse="|")
                                                          , as.character(x))))
      ind_list1 <- unlist(apply(ind_list1, 2, list), recursive = FALSE)
      
      ind_list2 <- apply(t_chk, 2, function(x)as.character(x)[grep(paste(as.character(x[duplicated(x)]),collapse="|")
                                                                   , as.character(x))])
      ind_list2 <- unlist(apply(ind_list2, 2, list), recursive = FALSE)
    }
    else if(dim(t_chk)[2] > 2){
      print("dimension gt 2")
      ind_list1 <- apply(t_chk, 2, function(x)which(grepl(paste(as.character(x[duplicated(x)]),collapse="|")
                                                          , as.character(x))))
      
      ind_list2 <- apply(t_chk, 2, function(x)as.character(x)[grep(paste(as.character(x[duplicated(x)]),collapse="|")
                                                                   , as.character(x))])
    }
    
    l1 <- lapply(seq_along(ind_list), function(y, n, i){ phen_mat[rownames(phen_mat) %in%
                                                                    y[[i]], colnames(phen_mat)
                                                                  %in% n[[i]] ]}, y = ind_list,n = names(ind_list))
    
    l2 <- Map(cbind.data.frame, ind_list2, l1, ind_list1)
    
    l2 <- lapply( l2 , setNames , nm = c("dup_samp", "phen_bin", "row_IBD") )
    l2 <- lapply(l2, function(x)x[x[,2] == 0,])
    l3 <- lapply(l2, function(x)x[duplicated(x[,1]),])
    l4 <- lapply(l3, function(x)tab_phe[x[,3],c(1:2)])
    l5 <- Map(cbind.data.frame, l3, l4)
    l6 <- lapply(l5, function(x)ifelse(x[,1] %in% x[,4], x[,5], x[,4]))
    l7 <- Map(cbind.data.frame, l5, l6)
    l7 <- lapply( l7 , setNames , nm = c("dup_samp", "phen_bin", "row_IBD", "FID1", "FID2", "replace") )
    
    #l8 <- lapply(l7, function(x)t_chk[x[,3],] <- x[,6])
    
    test <- lapply(t_chk, function(x)as.character(x))
    
    round1 <- mapply(x = test, y = l7, FUN = function(x, y) replace(x, y[,3], y[,6]))
    class(round1) <- "numeric"
    round1 <- cbind.data.frame(tab_phe[,c(1:3)], round1)
    return(round1)
  }
}

##test
r1 <- dup_adj(tab_phe)
r2 <- dup_adj(r1)
