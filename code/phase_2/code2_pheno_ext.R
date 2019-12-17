##New parallelised version
##use generate_inp_01.R to generate input phen_mat.rds
##One stop code for sample selection
phen_mat <- readRDS("~/Syn_collab/phen_mat.rds")
phen_mat <- as.data.frame.matrix(phen_mat)
IBD_dist <- read.delim("~/Syn_collab/test.genome", sep = "\t", header = T)
IBD_dist <- IBD_dist[,c(1,3,10)]
IBD_dist <- IBD_dist[order(IBD_dist$PI_HAT, decreasing = F),] 

##processing phen_mat and IBD_dist
#IBD_dist1 <- IBD_dist[IBD_dist$FID1 %in% rownames(phen_mat) & IBD_dist$FID2 %in% rownames(phen_mat),]
#phen_mat <- phen_mat[rownames(phen_mat) %in% unique(c(IBD_dist1$FID1, IBD_dist1$FID2)),]

library(igraph)
IBD_dist_G <- IBD_dist
colnames(IBD_dist_G) <- c("from", "to", "weight")
graph_wgt <- igraph::graph.data.frame(IBD_dist_G, directed = F)
graph_wgt_deg <- igraph::degree(graph_wgt, v = igraph::V(graph_wgt))
graph_wgt_deg_df <- as.data.frame(graph_wgt_deg)
graph_wgt_deg_df$samp <- rownames(graph_wgt_deg_df)
##advanced parallel

##function to get degree of two nodes simultaneously
#library(igraph)
#library(data.table)
#setDTthreads(threads = 1)
get_2nodes <- function(x1, x2){
  deg1 <- graph_wgt_deg_df[match(as.character(x1), graph_wgt_deg_df$samp), 1]
  deg2 <- graph_wgt_deg_df[match(as.character(x2), graph_wgt_deg_df$samp), 1]
  degs <- cbind.data.frame("deg1" = deg1, "deg2" = deg2)
 # names(degs) <- c("deg1", "deg2")
  return(degs)
}

pick_pinf_phen <- function(ibd_dat){
  t4 <- list()
  p1 <- ibd_dat[1]
  p2 <- ibd_dat[2]
  pi_hat <- ibd_dat[3]
  phen <- colnames(phen_mat)          #[1:5]
  # p_mat <- phen_mat[,colnames(phen_mat) %in% phen]
  c1 <- phen_mat[rownames(phen_mat) %in% c(as.character(p1), as.character(p2)), colnames(phen_mat) %in% phen]
  
  if(is.null(dim(c1)) | dim(c1)[1] < 2){
    print("sample pair not found")
    print(p1)
    print(p2)
    t4 <- NULL
  }
  else if(dim(c1)[1] == 2){
    c2 <- lapply(c1, function(x) as.data.frame(t(x)))
 #   print(p1)
 #   print(p2)
    c2 <- lapply(c2, function(x){colnames(x) <- c("p1_phen", "p2_phen"); return(x)}) 
    c3 <- mapply(cbind.data.frame,"p1" = p1,"p2" = p2, "pi_hat" = pi_hat, c2, SIMPLIFY=F)
    names(c3) <- names(c2)
    ##Add degree to dataframe
    t4 <- lapply(c3, function(x)cbind.data.frame(x, get_2nodes(x[1], x[2])))
   
  }
  return(t4)
}

#para_pheno <- apply(IBD_dist[3725:3740,], 1, pick_pinf_phen)

library(parallel)
  cl <- makeCluster(25)
  
  clusterExport(cl , c("phen_mat", "graph_wgt_deg_df", "get_2nodes"), envir = .GlobalEnv)
  system.time(para_pheno <- parApply(cl, IBD_dist, 1, pick_pinf_phen))
stopCluster(cl)

t1 <- para_pheno[unlist(lapply(para_pheno, length) != 0)]
t2 <- do.call(Map, c(rbind, t1))
t11 <- Map(cbind.data.frame, t2, "phen" = names(t2))
#saveRDS(t11, file = "~/Syn_collab/para_IBD_phen_comb.rds", compress = T)
##sample selection
#function for sample selection
samp_sel_dup <- function(data_comb){
sam_sel <- ifelse(data_comb$p1_phen > data_comb$p2_phen, data_comb$p1, 
                  ifelse(data_comb$p1_phen == data_comb$p2_phen & data_comb$deg1 > data_comb$deg2,
                         data_comb$p2, 
                         ifelse(data_comb$p1_phen == data_comb$p2_phen & data_comb$deg1 < data_comb$deg2,
                                data_comb$p1, data_comb$p2)))
data_comb$samp_sel <- sam_sel
phen_ibd_comb_uniq <- data_comb[!duplicated(data_comb$samp_sel),]
phen_ibd_comb_uniq <- phen_ibd_comb_uniq[!is.na(phen_ibd_comb_uniq$samp_sel),]
return(phen_ibd_comb_uniq)
}

phen_fin_sel <- lapply(t11, function(x)samp_sel_dup(x))
saveRDS(phen_fin_sel, file = "~/Syn_collab/para_IBD_phen_comb.rds", compress = T)

##QC
# sapply(phen_fin_sel, dim)
# colSums(phen_mat)[1:5]
# sapply(phen_fin_sel,function(x)sum(x[,4]))
# sapply(phen_fin_sel,function(x)sum(x[,5]))
# 
# pos_cont <- read.delim("~/Syn_collab/selected.txt", header = T, sep = "")
# 
# apply(phen_mat[,1:5], 2, function(x)table(rownames(phen_mat)[which(x == 1)] %in% pos_cont$FID))
