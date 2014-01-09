# Nitesh Turaga
# Aclust run for Shyam Biswal for 6 samples.

date()
sessionInfo()


####################################################################
# 1.
# Load required packages and datasets
####################################################################
library(Aclust)

library(minfi) # To get the annotation from minfi

setwd("~/TestRun/ShyamBiswal/")

loadDir = "~/TestRun/ShyamBiswal/circos_plot_test/"
storageDir = "~/TestRun/ShyamBiswal/"

load(file.path(loadDir,"p.rda"))
load(file.path(loadDir,"otherstuff02.rda"))
load(file.path(loadDir,"pd02.rda"))

####################################################################

# Beta values
betas = getBeta(object)

#Exposure
exposure = ifelse(test = (pd$status == "H460_parent"), yes = 1,no = 0)

####################################################################
# ANNOTATION
annot = getAnnotation(object)

# Curate annotation to fit function
annot_df = as.data.frame(annot)

# Create New data table for the Aclustering annotation 

annot_aclust = data.frame(IlmnID = as.character(rownames(annot_df)),Infinium_Design_Type = annot_df$Type)
annot_aclust$CHR = gsub(pattern = "chr",replacement="",x= annot_df$chr)
annot_aclust$Coordinate_36 = annot_df$pos
annot_aclust$Gene_Name = annot_df$UCSC_RefGene_Name
annot_aclust$UCSC_RefGene_Group = annot_df$UCSC_RefGene_Group
annot_aclust$UCSC_CpG_Islands_Name = annot_df$Islands_Name
annot_aclust$Relation_to_UCSC_CpG_Island = annot_df$Relation_to_Island
annot_original = annot
annot = annot_aclust

#Cast annot as a data table
annot = as.data.table(annot)

######################################################################


clusters.list = assign.to.clusters(betas=betas,annot=annot)
length(clusters.list)

# Covariates can be categorical or continuous, but cannot have intercept term in the,
# nor sure why it was designed this way
covariates = model.matrix(~pd$Race)
covariates = matrix(covariates[,c(2,3,4)],ncol=3)

rownames(covariates) = colnames(betas)

GEE.results.clusters = GEE.clusters(betas, clusters.list, exposure, covariates, id = colnames(betas), working.cor = "ex")
save(clusters.list,GEE.results.clusters,file = "../objs/dec11_gee.rda")
top.clusters.summary = summarize.top.clusters(betas,covariates=covariates,exposure,id = colnames(betas),GEE.results.clusters,"results.tex",annot = annot)

save.image(file = "../objs/dec11_7a.Rdata")

# 
# ?Aclust
# data(betas.7) ## upload methylation data
# exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
# covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
# rownames(covariates) <- colnames(betas.7)
# 
# data(annot)  ## load annotation created using the IlluminaHumanMethylation450k.db package on July 2013
# debug(assign.to.clusters)
# clusters.list <- assign.to.clusters(betas.7, annot)
# undebug(assign.to.clusters)
# 
# head(clusters.list)
# class(clusters.list)
# covariates
# 
# 
# GEE.results.clusters <- GEE.clusters(betas.7, clusters.list, exposure, covariates, id = colnames(betas.7), working.cor = "ex")
# top.clusters.summary <- summarize.top.clusters(betas.7, covariates, exposure, id = colnames(betas.7), GEE.results.clusters, "results.tex", annot= annot)
