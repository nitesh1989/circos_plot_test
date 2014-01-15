# RCircos Plots for Shyam Biswal.
# Nitesh Turaga
date()
sessionInfo()

#############################################
# Step 0 : Set up path and data 
#############################################

# Set path
# my.path = "~/Documents/JHMI-Research/charmData/ShyamBiswal/" # Local path
my.path = "~/TestRun/ShyamBiswal/circos_plot_test/"
setwd(my.path)
# Import libraries required
library(RCircos)
library(BiocGenerics)
library(Biobase)

# Import Peaks data from Shywam Biswals
parent = read.csv("parent_genes_distFromTSS.csv",stringsAsFactors = FALSE)
knock = read.csv("knockin_genes_distFromTSS.csv",stringsAsFactors = FALSE)


################################################
#Step 1: Curate data for RCircos package
################################################

curate.devaPeaks.to.rcircos = function(peaks.df) {
    peaks.df.x = data.frame(peaks.df[,c(3:23,1,2)])
    cols.peaks = colnames(peaks.df.x) 
    cols.peaks[2] = "START"
    cols.peaks[3] = "END"
    colnames(peaks.df.x) = cols.peaks
    return(peaks.df.x)
}

parent.x = curate.devaPeaks.to.rcircos(parent)
knock.x = curate.devaPeaks.to.rcircos(knock)


# Curate knocking pathways first
# Curate to Rcircos

collect.pathways = function(status,path) {
    files = list.files(path,full.names=T)
    new_df = data.frame()
    for (i in 1:length(files)) {     
        x = read.csv(files[i])
        pathwayName = gsub(pattern = paste0("_",status,".csv"),replacement="",x=(gsub(".+/","",files[i])))
        x = x[,c(4:7,1:3,7:22)]
        knock_df = data.frame(x,pathway = pathwayName,pathwayGene = paste0(x$name,"-",pathwayName))
        new_df = rbind(new_df,knock_df)
    }
    return(new_df)
}

knock_df = collect.pathways(status = "H460_knock","Pathways_knock")
parent_df = collect.pathways(status= "H460_parent","Pathways_parent")

collect.genes = function(row,pathway){pathwayrow$name}



BiocGenerics::table(knock_df$pathway,knock_df$name)

################################################
# Step 2: Initialize base of Circos plot
################################################
# Call RCircos Cytoband for chromosomes
data(UCSC.HG19.Human.CytoBandIdeogram)
cyto.info = UCSC.HG19.Human.CytoBandIdeogram;

# Exclude sex chromosomes
# chr.exclude = c("chrX","chrY");
chr.exclude = NULL;
# Determine tracks inside and outside
tracks.inside = 4;
tracks.outside = 6;

# Initialize RCircos Core components
RCircos.Set.Core.Components(cyto.info, chr.exclude=chr.exclude,tracks.inside=tracks.inside,tracks.outside= tracks.outside);

################################################
# Step 3: Set R Circos plot parameters
################################################

rcircos.cyto = RCircos.Get.Plot.Ideogram()
rcircos.position = RCircos.Get.Plot.Positions()
# rcircos.gene.locations = RCircos.Get.Gene.Label.Locations(genomic.data=parent.x )

params = RCircos.Get.Plot.Parameters()

# params$radius.len = 0.5;
# params$base.per.unit =1000;
# params$track.padding = 0.02;
# params$track.height = 0.2;
# 
# params$chrom.width = 0.1;
# params$chr.ideog.pos = 1.1;
# params$chr.name.pos = 1.25;
# 
# params$text.size = 0.4;

# ORIGINAL
# params$radius.len = 2.0;
# params$base.per.unit =1000;
# params$track.padding = 0.02;
# params$track.height = 0.2;
# 
# params$chrom.width = 0.2;
# params$chr.name.pos = 2.24;
# 
# params$text.size = 0.6;
# 

params$radius.len = 2.0;
params$base.per.unit =1000;
params$track.padding = 0.02;
params$track.height = 0.2;

params$chrom.width = 0.2;
params$chr.name.pos = 2.24;

params$text.size = 0.6;

RCircos.Reset.Plot.Parameters(params);

RCircos.List.Parameters() # For reference in Rout file

################################################
# Step 4: Make RCircos plots
################################################

out.file = "BiswalCircos5.png"
png(file = out.file, height = 30, width = 30,res=600,units = "in")
#plot(main="Inner Tracks: H460 Knocking-Genes; Outer Tracks: H460 Parent-genes")
RCircos.Set.Plot.Area()
# Plot chromosome ideogram
RCircos.Chromosome.Ideogram.Plot();

# Information for track 1 "inside" - Plot gene names(column 13)
RCircos.Gene.Name.Plot(gene.data= knock.x,name.col= 13,track.num=1,side="out")
RCircos.Gene.Connector.Plot(genomic.data=knock.x,track.num=2,side = "out")

#Information for track 1 "outside" - Plot gene names(column 15)
RCircos.Gene.Name.Plot(gene.data=parent.x,name.col= 13,track.num=2,side="in")
RCircos.Gene.Connector.Plot(genomic.data=parent.x,track.num=1,side="in")

# Plot Peaks like RCircos Line data

RCircos.Line.Plot(line.data=knock.x,data.col=4,track.num=3,side="out")
RCircos.Line.Plot(line.data=parent.x,data.col=4,track.num=3,side = "in")

################################################
# Step 5: Plot Pathways
################################################

# Plot in R Circos

RCircos.Gene.Name.Plot(gene.data= knock_df,name.col= which(colnames(knock_df) == "pathwayGene"),track.num=4,side="out")
# RCircos.Gene.Connector.Plot(genomic.data=knock_df,track.num=4,side = "out")

RCircos.Gene.Name.Plot(gene.data= parent_df,name.col= which(colnames(parent_df) == "pathwayGene"),track.num=4,side="in")
# RCircos.Gene.Connector.Plot(genomic.data=parent_df,track.num=4,side = "in")


dev.off()

