##beta distribution for weighing SKAT variants based on MAFs
##Function to explicitly return SKAT weights from MAF
MAF_wt_def <- function(var_maf){
weight <- dbeta(var_maf, 1, 25)
return(weight)
}


 