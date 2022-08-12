genus_tree <- function(backbone_tree,comm_data,n) {
  
  final_tree <- list()
  i=1
  backbone_tree <- as.multiPhylo(backbone_tree)
  
  for (j in (1:length(backbone_tree))) {
    for (k in (1:n)) {
    message("making ",k," tree for posterior sample ",j)
    dummy_tree <- backbone_tree[[j]]
    
    genera<-unique(sapply(strsplit(rownames(comm_data),"_"),function(x) x[1]))
    ii<-sapply(genera,function(x,y) grep(x,y)[sample(length(grep(x,y)),1)],y=dummy_tree$tip)
    dummy_tree<-drop.tip(dummy_tree,setdiff(dummy_tree$tip.label,dummy_tree$tip[ii]))
    
    dummy_tree$tip.label<- sapply(strsplit(gsub("\"","",dummy_tree$tip.label),"\\."),function (x) x[1])
    dummy_tree <- genus.to.species.tree(dummy_tree,rownames(comm_data))
    final_tree[[i]] <- dummy_tree
    i=i+1
  }
  }
  class(final_tree) <- "multiPhylo"
  return(final_tree)
}