
#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################

############################# within gene network analysis

library(igraph)
library(stringi)
library(data.table)
library(dplyr)
library(minet)

setwd("your path here to working directly")
dir.create("Graphxtics")
dir.create("Betweeness")   #store MI

#Specify where the gene files are
Path2GeneFiles = "your path here " # original file source

Path2Results <- "./Graphxtics"
outdir=Path2Results

Path2ResultsMI <- "./Betweeness"
outdirMI=Path2ResultsMI


files <- list.files(path = Path2GeneFiles, pattern="\\.csv$")


#graphx = c()

gene<-density<-mean_dist<-transitivity<-edge_dens<-vertex_con<-edge_con<-NULL

for(file in files)
{
  data_in <- read.table(file.path(Path2GeneFiles, file), header=F, check.names=FALSE, stringsAsFactors=FALSE)
  #print(head(data))
  colnames(data_in) = c("SNP1", "SNP2","mi")
  print(head(data_in))
  data = data_in
  
  ##########################################################
  ### make adjacency matrix
  #########################################################
  # get names for row and columns
  nameVals = sort(unique(unlist(data[1:2])))
  # construct 0 matrix of correct dimensions with row and column names
  myMat = matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(data[c("SNP1", "SNP2")])] <- data[["mi"]] 
  adjmat = abs(myMat)
  #adjmat
  diag(adjmat) = 0     #make diagonal zero i.e no info btn the same snp
  print(file)
  print(adjmat)

  
  ############################################################
  ## select meaningfull edges using minet pkg and summarise graph
  ###########################################################
  # ** QUICK WAY TO RUN A SECOND ALGORITHM
  #graph_mrnet = aracne(adjmat, eps=0 )   # NOTICE THE CHANGE IN ALGORITHM

  #--------------------------------------------------------------------------
  # mrnet: Maximum Relevance Minimum Redundancy
  graph_mrnet = mrnet(adjmat)
  
  datgraph_mrnet = graph_from_adjacency_matrix(graph_mrnet, mode = "undirected", weighted = TRUE,
                                               diag = FALSE)
  #remove loops
  datgraph_mrnet = simplify(datgraph_mrnet, remove.multiple=TRUE, remove.loops=TRUE)
  
  gene[[file]] = file
  # -----------------------------   summaries of interest
  density[[file]] = graph.density(datgraph_mrnet,loop=FALSE)  #Density
  
  mean_dist[[file]] = mean_distance(datgraph_mrnet)  #Average Path Length
  
  transitivity[[file]] = transitivity(datgraph_mrnet)    #Clustering Coefficeint
  
  edge_dens[[file]] = edge_density(datgraph_mrnet, loops=F) #number of edges/no.of posible edges
  
  vertex_con[[file]] = vertex_connectivity(datgraph_mrnet) #number of edges/no.of posible edges
  edge_con[[file]] = edge_connectivity(datgraph_mrnet) #number of edges/no.of posible edges
  
  #deg <- degree(datgraph_mrnet, mode="all")
  
  #print(degree(datgraph2, mode = "all"))    #Degree: In, Out, All Centrality
  
  
  #print(degree(datgraph2, mode = "all"))    #Degree: In, Out, All Centrality
  # betweeness centrality

  # freeze everything here ----------------------------------------------------------------------------
  snpbetw_centr = betweenness(datgraph_mrnet, directed=F, weights=NA)
  snpbetw_centr = data.frame(snpbetw_centr)
  ##dim(snpbetw_centr)
  snpsbetw_centrDF <- tibble::rownames_to_column(snpbetw_centr, "SNP")
  
  outputfile_MI=paste0(outdirMI,"/",file)
  write.table(snpsbetw_centrDF,outputfile_MI,row.names=F, col.names=TRUE,quote=FALSE,sep=",")
  # freeze end -----------------------------------------------------------------------------------------
}

graphxx = data.frame(gene,density,mean_dist,transitivity,edge_dens,vertex_con,edge_con)
write.csv(graphxx, "network_analysis4000snps.csv")

#graphxx = data.frame(graphx)
graphxx = setNames(graphxx, c("gene", "density", "mean_dist", "transitivity", "edge_dens", "vertex_con", "edge_con"))

outputfile=paste0(outdir,"/",file)
write.table(graphxx,outputfile,row.names=F, col.names=TRUE,quote=FALSE,sep=",") 
# write in working directory
write.csv(graphxx, "network_analysis4000snps.csv")



print("---------------- The script has run fully ------------------------------------------")