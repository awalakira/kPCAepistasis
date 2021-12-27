#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################

# SET WORKING DIRECTORY FIRST, DEFINE CORRECT SOURCE OF FILES

################# Bayesian modelling, Antoneli et al 2020 #################################################

library(NLinteraction)

print("modelling 2 way and 3 way interactions using ns2")
setwd("your path here")

# file of kpc1 for each gene
data_genesFullx = read.csv("your path here/kpca1_genes_final.csv", header=TRUE)   
print(data_genesFullx[1:5, 1:8])
##remove .csv from column names
names(data_genesFullx) = gsub(".csv","", names(data_genesFullx), fixed=TRUE)     
print(data_genesFullx[1:5, 1:8])

# data_genesFullxmat = data_genesFullx[,3:ncol(data_genesFullx)]
data_genesFull = data_genesFullx[,3:ncol(data_genesFullx)]


#pheno = read.csv("your path here")

PHENOTYPE = data_genesFullx$PHENOTYPE

### PHENOTYPE = data_genesFullx$PHENOTYPE

X = as.matrix(data_genesFull)
Y = as.matrix(PHENOTYPE)

NLmod2 = NLint(Y=Y, X=X, C=NULL, nIter=20000, nBurn=2, thin=5, nChains=2, ns=1)

NLmod = NLmod2

##So we can now evaluate the WAIC of each model
print("WAIC for ns=1 model")
print(NLmod2$waic)

################===========================================
#Posterior inclusion probabilities

pdp = NLmod$MainPIP
genenam = as.vector(colnames(data_genesFull))
gname = data.frame(genenam)
post_inc_props = cbind(gname, pdp)

## write out gene names in a txt
#writeLines(genenam, "genes_BayesianTrial.txt")

write.csv(post_inc_props, "post_inc_probs_Iter20k_ns1.csv")
##===========================================================================

####### =======================================================================
#We now look at the matrix of two-way interaction probabilities.

intMat = NLmod$InteractionPIP
write.csv(intMat, "intMat_2way_iter20kns1.csv")
#### =========================================================================

## =========================================================================
# image
## save plot interactionProbabilities Plot
pdf(file="int_probs_2way_20k_ns1.pdf")
plotInt(NLmod = NLmod)
dev.off()

## ========================================================================

## write out gene names in a txt
writeLines(genenam, "genes_BayesianTrial.txt")

# =================================================================
## save trace plot
pdf(file="Trace_plot_20kns1.pdf")
plot(NLmod2$posterior$sigma[1,], type="l")
lines(NLmod2$posterior$sigma[2,], col=2)
dev.off()

# =================================================================




##### takes many days if you have many genes 
##########################============== trying out loop 3 way interactions================================
#inter=c()
#var1=c()
#var2=c()
#var3=c()
#for(i in 1:ncol(X)){
#  for(j in 1:ncol(X)){
#    for(k in 1:ncol(X)){
#      if(i!=j & i!=k & j!=k){
#        var1 = c(var1,colnames(X[i]))
#        var2 = c(var2,colnames(X[j]))
#        var3 = c(var3,colnames(X[k]))
#        int = InteractionProb(NLmod=NLmod2, Xsub=c(i,j,k))
#        inter = c(inter,int)
#      }
#    }
#  }
#}

#intdf_3way = data.frame(var1,var2, var3, inter)
#write.csv(intdf_3way, "intdf_3way_4000snps_20k_ns1.csv")
##########==============================================================