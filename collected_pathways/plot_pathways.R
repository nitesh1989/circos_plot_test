setwd("~/Documents/JHMI-Research/charmData/ShyamBiswal/collected_pathways/")

knock_edit = read.csv("knock_TSS_change_edit_for_latex.csv")
parent_edit = read.csv("parent_TSS_change_edit_for_latex.csv")
x = as.data.frame(table(knock_edit$pathway))
y = as.data.frame(table(parent_edit$pathway))
#png(filename="knock.png",width=1500,height=1500,res=350)
pdf(file="knock.pdf")
par(mar = c(27,7,7,7))
barplot(table(knock_edit$name,knock_edit$pathway),col=rainbow(as.factor(knock_edit$name)),ylim = c(0,40),names.arg=x$Var1,cex.names=0.7,las = 3,cex.axis=0.5)

midpoints = barplot(table(knock_edit$name,knock_edit$pathway),col = heat.colors(n=length(unique(knock_edit$name))),las=3,ylim = c(0,40))
text(x = midpoints)


library(lattice)

barchart(knock_edit$pathway~knock_edit$name|factor(knock_edit$pathway),data = knock_edit,stack =T,auto.key =T)

dev.off()

