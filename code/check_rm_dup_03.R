##duplication check
# ##check duplication of samples in randomly picked cases
##Function for checking duplicates and reducing the numbers taking into account the phenotypes
library(rlist)
phen_mat <- readRDS("~/Syn_collab/phen_mat.rds")
dup_adj <- function(tab_phe){
 
  print(dim(tab_phe))
  tab_phe <- as.data.frame(tab_phe)
  
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
    t_chk <- t_chk[,colnames(t_chk) %in% colnames(phen_mat)]
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
      # ind_list1 <- apply(t_chk, 2, function(x)which(grepl(paste(as.character(x[duplicated(x)])),collapse="|")
      #                                                    , as.character(x)))
      ind_list1 <- apply(t_chk, 2, function(x)which(as.character(x) %in% as.character(x[duplicated(x)])))
      
      
      # ind_list2 <- apply(t_chk, 2, function(x)as.character(x)[grep(paste(as.character(x[duplicated(x)]),collapse="|")
      #                                                               , as.character(x))])
      ind_list2 <- apply(t_chk, 2, function(x)as.character(x)[which(as.character(x) %in% 
                                                                      as.character(x[duplicated(x)]))])
    }
    
    
    ##Add pheno name into the list
    # l1 <- mapply(append, names(ind_list2), ind_list2, SIMPLIFY = FALSE)
    #  l1 <- lapply(l1, function(x)phen_mat[rownames(phen_mat) %in% unique(x[-1]), colnames(phen_mat) %in% x[1]])
    
    #  l2 <- Map(cbind.data.frame, ind_list2, ind_list1)
    l2 <- Map(cbind.data.frame, ind_list2, ind_list1, names(ind_list2))
    
    l3 <- lapply(l2, function(x)phen_mat[rownames(phen_mat) %in% as.character(x[,1]), 
                                         colnames(phen_mat) %in% unique(as.character(x[,3]))])
    l3_names <- lapply(l3, function(x)names(x))
    l3 <- Map(cbind.data.frame, l3, l3_names)
    l4 <- mapply(x = l3, y = l2, FUN = function(x, y) x[match(y[,1], x[,2]),1])
    
    l5 <- Map(cbind.data.frame, l2, l4)
    l5 <- lapply( l5 , setNames , nm = c("dup_samp","row_IBD", "phen", "phen_bin") )
    l5 <- lapply(l5, function(x)x[!is.na(x[,4]),])
    
    #l5 <- lapply(l5, function(x)x[x[,2] == 0,])
    
    ##Duplicated cases will be targeted using indices from row_IBD attribute
    l6 <- lapply(l5, function(x)x[duplicated(x[,1]),])
    #  l6_unique <- lapply(l5, function(x)x[!duplicated(x[,1]),])
    
    ##Extract IBD pairwise information  
    l7 <- lapply(l6, function(x)tab_phe[x[,2],c(1:2)])
    l7 <- Map(cbind.data.frame, l6, l7)
    
    ##Add phenotype of the paired case that is not present in dup_samp column
    l7_phen_pair <- lapply(l7, function(x)ifelse(as.character(x[,1]) == as.character(x[,5]), x[,6], x[,5]))
    l7_pp <- Map(cbind.data.frame, l7, "paired_phen" = l7_phen_pair)
    l7_pp1 <- lapply(l7_pp, function(x)phen_mat[rownames(phen_mat) %in% x[,7], colnames(phen_mat) %in% unique(x[,3])])
    l7_pp1_names <- lapply(l7_pp1, function(x)names(x))
    l7_pp1 <- Map(cbind.data.frame, l7_pp1, l7_pp1_names)
    l7_pp2 <- mapply(x = l7_pp1, y = l7_pp, FUN = function(x, y) x[match(y[,7], x[,2]),1])
    
    l7_pp <- Map(cbind.data.frame, l7_pp, "pp_phen_bin" = l7_pp2)
  
    l8 <- lapply(l7_pp, function(x)ifelse(x[,4] == x[,8], x[,7], 
                                          ifelse(x[,4] > x[,8], x[,1], x[,7])))

    
    l9 <- Map(cbind.data.frame, l7_pp, "replace" = l8)
    
    l10 <- lapply(t_chk, function(x)as.character(x))
    
    round1 <- mapply(x = l10, y = l9, FUN = function(x, y) replace(x, y[,2], y[,9]))
    class(round1) <- "numeric"
    round1 <- cbind.data.frame(tab_phe[,c(1:3)], round1)
    return(round1)
  }
}
