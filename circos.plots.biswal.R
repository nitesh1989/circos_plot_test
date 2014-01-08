# RCircos Plots for Shyam Biswal.
# Nitesh Turaga
date()
sessionInfo()

#############################################
# Step 0 : Set up path and data 
#############################################

# Set path
my.path = "~/TestRun/ShyamBiswal/"

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
tracks.inside = 5;
tracks.outside = 4;

# Initialize RCircos Core components
RCircos.Set.Core.Components(cyto.info, chr.exclude=chr.exclude,tracks.inside=tracks.inside,tracks.outside= tracks.outside);

################################################
# Step 3: Set R Circos plot parameters
################################################

rcircos.cyto = RCircos.Get.Plot.Ideogram()
rcircos.position = RCircos.Get.Plot.Positions()
# rcircos.gene.locations = RCircos.Get.Gene.Label.Locations(genomic.data=parent.x )

params = RCircos.Get.Plot.Parameters()

params$radius.len = 3.0;
params$base.per.unit =1000;
params$track.padding = 0.02;
params$track.height = 0.3;

params$chrom.width = 0.2;
params$chr.name.pos = 2.24;

params$text.size = 0.6;



RCircos.Reset.Plot.Parameters(params);

RCircos.List.Parameters() # For reference in Rout file

################################################
# Step 4: Make RCircos plots
################################################

out.file = "BiswalCircos.png"
png(file = out.file, height = 22, width = 22,res=500,units = "in")
plot(main="Inner Tracks: H460 Knocking-Genes; Outer Tracks: H460 Parent-genes")
RCircos.Set.Plot.Area()
# Plot chromosome ideogram
RCircos.Chromosome.Ideogram.Plot();

# Information for track 1 "inside" - Plot gene names(column 13)
RCircos.Gene.Name.Plot(gene.data= knock.x,name.col= 13,track.num=2,side="in")
RCircos.Gene.Connector.Plot(genomic.data=knock.x,track.num=1,side = "in")

#Information for track 1 "outside" - Plot gene names(column 15)
RCircos.Gene.Name.Plot(gene.data=parent.x,name.col= 13,track.num=2,side="out")
RCircos.Gene.Connector.Plot(genomic.data=parent.x,track.num=1,side="out")

# Plot Peaks like RCircos Line data

RCircos.Line.Plot(line.data=knock.x,data.col=4,track.num=3,side="in")
RCircos.Line.Plot(line.data=parent.x,data.col=4,track.num=3,side = "out")


dev.off()

