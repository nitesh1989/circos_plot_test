# Peak file analysis in Shyam Biswal data set

date()
sessionInfo()


setwd("~/TestRun/ShyamBiswal/")
load("objs/pd01.rda")

library(BiocGenerics)

knock = pd[pd$SAMPLE_NAME =="H460_knockin",]


knock1 = read.csv("rawData/5000to1500dnaMethFlow/542654A01_Slot9_2012-07-25_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F)
knock2 = read.csv("rawData/5000to1500dnaMethFlow/542671A01_Slot6_2012-07-24_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F)
knock3 = read.csv("rawData/5000to1500dnaMethFlow/542707A01_Slot12_2012-07-25_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F)


total_knock = rbind(knock1,knock2,knock3)

#total_knock= total_knock[! duplicated(total_knock$Name),]
#ord = order(total_knock$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,decreasing=FALSE)
#total_knock_x = total_knock[ord,]

genes_knock = BiocGenerics::intersect(BiocGenerics::intersect(knock3$Name,knock1$Name),knock2$Name)

write.csv(as.data.frame(genes_knock),file = "knockin_genes_intersect_TSS_change.csv")

parent = pd[pd$SAMPLE_NAME =="H460_parent",]


parent1 = read.csv("rawData/5000to1500dnaMethFlow/542586A01_Slot9_2012-07-24_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F,fill=T)
parent2 = read.csv("rawData/5000to1500dnaMethFlow/542737A01_Slot13_2012-07-24_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F,fill=T)
parent3 = read.csv("rawData/5000to1500dnaMethFlow/542857A01_Slot18_2012-07-24_H460_BLOCK1_DNA_Methylation_Peak_Results_Peaks_mapping.csv",stringsAsFactors=F,fill=T)


genes_parent = BiocGenerics::intersect(BiocGenerics::intersect(parent3$Name,parent1$Name),parent2$Name)

total_parent = rbind(parent1,parent2,parent3)

# total_parent = total_parent[! duplicated(total_parent$Name),]
# ord = order(total_parent$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,decreasing=FALSE)
# total_parent_x = total_parent[ord,]


write.csv(as.data.frame(genes_parent),file = "parent_genes_intersect_TSS_change.csv")


common_genes = BiocGenerics::intersect(genes_knock,genes_parent)

write.csv(as.data.frame(common_genes),file = "common_knockin_parent_TSS_change.csv")


# Remove common genes

genes_knock_no_common = BiocGenerics::setdiff(genes_knock,common_genes)
genes_parent_no_common = BiocGenerics::setdiff(genes_parent,common_genes)


# Annotate genes from set difference
annot_knock = total_knock[match(x = genes_knock_no_common,table=total_knock$Name),]
annot_knock = (as.data.frame(annot_knock))

annot_parent = total_parent[match(x = genes_parent_no_common,table = total_parent$Name),]
annot_parent = as.data.frame(annot_parent)

# Get summary

table(annot_knock$FEATURE_TRACK)
summary(annot_knock$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT)

table(annot_parent$FEATURE_TRACK)
summary(annot_parent$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT)    


library(ggplot2)


png(filename="distribution_knocking_annotated_TSS_change.png")
qplot(annot_knock$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,main = "Distance from TSS, H460_knock",xlab = "Distance from TSS in bp")
dev.off()


png(filename="distribution_parent_annotated_TSS_change.png")
qplot(annot_parent$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,main ="Distance from TSS, H460_parent",xlab = "Distance from TSS in bp")
dev.off()


write.csv(annot_knock[order(annot_knock$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,decreasing=TRUE),],file = "knockin_genes_distFromTSS_change.csv")

write.csv(annot_parent[order(annot_parent$SHORTEST_DISTANCE_FROM_FEATURE_TO_DATA_POINT,decreasing=TRUE),],file = "parent_genes_distFromTSS._change.csv")
