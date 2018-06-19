########load packages
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
#load libraries
library(ape)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(miLineage)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(venneuler)



###### Read and process OTUs ######

#set working directory
"~/Users/sarahbecker/Desktop/THESIS/RDirectory/"

OTU <- read.table("/Users/sarahbecker/Desktop/THESIS/RDirectory/otu.final", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
unique(OTU$label)   # only takes value 0.03
unique(OTU$numOtus) # only takes value 859

rownames(OTU) <- OTU$Group
OTU <- OTU[, -c(1:3)]       # now subject labels are row names
OTU <- OTU[-which(rownames(OTU) %in% c("PCR1", "PCR2")), ]  # I believe these are technical controls
OTU[1:10,1:10]   # take a peek at OTU counts

taxonomy <- read.table("/Users/sarahbecker/Desktop/THESIS/RDirectory/taxonomy.final", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(taxonomy)
tax.key <- taxonomy
tax.key$Kingdom = tax.key$Phylum = tax.key$Class = tax.key$Order = tax.key$Family = tax.key$Genus = NA
for (i in 1:nrow(tax.key)) {
  this.tax <- strsplit(tax.key$Taxonomy[i], split = ";")[[1]]
  tax.key[i, 4:9] <- rev(c(this.tax, rep(NA, 6-length(this.tax))))
}
head(tax.key)

metadata <- data.frame(sample_id = rownames(OTU), subj_id = NA, group = NA, time= NA,
                       stringsAsFactors = FALSE)
metadata$time <- substr(metadata$sample_id, nchar(metadata$sample_id), nchar(metadata$sample_id))
metadata$group <- substr(metadata$sample_id, 1, 2)
metadata$subj_id <- paste("SUBJ", gsub('[:A-Z:]+', '', metadata$sample_id), sep = "")

OTU.genus <- aggregate(t(OTU), by = list(tax.key$Taxonomy), sum)
rownames(OTU.genus) <- OTU.genus[,1]
OTU.genus <- OTU.genus[,-1]
OTU.genus <- t(OTU.genus)

genus.key <- tax.key[sapply(colnames(OTU.genus), FUN = function(x) which(tax.key$Taxonomy == x)[1]), ]

rare <- which(apply(OTU.genus, 2, sum) < 2)
OTU.genus <- OTU.genus[, -rare]
genus.key <- genus.key[-rare, ]

####regenerate OTU.clean table
row.names(OTU) = OTU$Group

OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]

row.names(taxonomy) = taxonomy$OTU
tax.clean = taxonomy[row.names(taxonomy) %in% colnames(OTU.clean),]
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "Strain", "OTU"))]
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]


#(Anna descriptive analysis)

######## Stacked bar plots (time1 and time2, separated by treatment group (4))

## aggregate counts by family (otherwise there are too many to plot)
family.taxonomy <- apply(tax.key, 1, FUN = function(x) paste(x[c(9:5)], collapse = ";"))
family.tab <- aggregate(t(OTU), by = list(family.taxonomy), FUN = sum)
rownames(family.tab) <- family.tab$Group.1
family.tab <- family.tab[, -1]

## aggregate counts by treatment group and time (so we can compare treatments)
bytrt <- aggregate(t(family.tab), by = list(metadata$group, metadata$time), FUN = sum)
rownames(bytrt) <- paste(bytrt[,1], bytrt[,2], sep = "")
bytrt <- bytrt[, -c(1:2)]
bytrt[1:3,1:3]

## combine rare families (up to 1000 reads total -- actually ranges between 0 and 138)
famtot <- apply(bytrt, 2, sum)
keep <- which(famtot > 1000)    ## 44
other <- bytrt[, !c(1:ncol(bytrt)) %in% keep]  ## "other" genus
Other <- apply(other, 1, sum)

## turn into proportions
subfam <- cbind(bytrt[,keep], Other)
subfam.prop <- subfam   ## individual-level proportions
for(i in 1:nrow(subfam)){
  subfam.prop[i,] <- subfam[i,]/sum(subfam[i,])
}
subfam.prop <- t(as.matrix(subfam.prop))[, c(1,5,2,6,3,7,4,8)]

## figure out labels
family.label <- unname(sapply(rownames(subfam.prop), FUN = function(x) strsplit(x, split = ";")[[1]][5]))
family.label[20] <- "Other"
family.label <- gsub("_", " ", family.label)
subfam.prop.temp=subfam.prop


#####POSTER PLOT**********************************
#change names from column headers in metadata to treatment group names 
colnames(subfam.prop.temp)=c("LF Pre","LF Post ","HF Pre","HF Post","Saccharin Pre", "Saccharin Post", "Stevia Pre","Stevia Post")

## use subfam.prop to make stacked bar plots
jpeg("family-stackedbar5.jpeg", width = 1500, height = 800, quality = 200, res=150)
par(mar=c(9.1,5.1,4.1,15.1))
barplot(subfam.prop.temp, col = rainbow(20),
        space = c(0.2, 0.2, 0.5, 0.2, 0.5, 0.2, 0.5, 0.2),
        main = "Proportion of OTU by Family",
        xlab = "",
        ylab = "Proportion of Total Reads",
        las=2,
        legend.text = family.label,
        args.legend = list(x = "right", inset = c(-0.45, 0), bty = 'n'))
mtext("Treatment Group", side=1, adj=0.5, line=7)
dev.off()



#########generate a Shannon index
OTU.physeq = otu_table(as.matrix(OTU.clean[-c(40,41),]), taxa_are_rows=FALSE)

tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(metadata)
sample_names(meta.physeq) = sample_names(OTU.physeq)
sample


physeq.alpha = phyloseq(OTU.physeq, tax.physeq, meta.physeq)

set.seed(1)

#creates alpha diversity column
sample_data(physeq.alpha)$shannon.physeq <- estimate_richness(physeq.alpha, measures="Shannon")

#plots alpha diversity
plot_richness(physeq.alpha, "group", measures="Shannon")

metadata$shannon.vegan <- diversity(OTU, index="shannon")

#library(lme4)
#rm.shannon.all = lmer(shannon ~ group+time + (1|Animal), data=metadata)
#summary(rm.shannon.all)

#Correlations of treatment by using subfam to compare OTU at family level within a treatment group
wilcox.test(as.numeric(subfam[1,]),as.numeric(subfam[5,]),paired=TRUE)
wilcox.test(as.numeric(subfam[2,]),as.numeric(subfam[6,]),paired=TRUE)
wilcox.test(as.numeric(subfam[3,]),as.numeric(subfam[7,]),paired=TRUE)
wilcox.test(as.numeric(subfam[4,]),as.numeric(subfam[8,]),paired=TRUE)

#Bacteriodales_S47
t.test(as.numeric(subfam[1,3]),as.numeric(subfam[5,3]),paired=TRUE)

#correlations at genus level corresponding to treatment group, HF
library(xlsx)
t1dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "A"), 2]
t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), 2]
wilcox.test(t1dat, t2dat, paired = TRUE)

#loop to calculate for all genus
pvals.hf <- c() 
for (i in 1:ncol(OTU.genus)) {
  t1dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "A"), i]
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  this.p <- wilcox.test(t1dat, t2dat, paired = TRUE)$p.value
  pvals.hf <- c(pvals.hf, this.p) 
} 
#export to excel
data.frame(colnames(OTU.genus), pvals.hf)
write.xlsx(data.frame(colnames(OTU.genus), pvals.hf), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.hf.xlsx")


##difference in means
pvals.hf <- c() 
means.group1 <- c() 
means.group2 <- c() 

for (i in 1:ncol(OTU.genus)) {
  t1dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "A"), i]
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  this.p <- wilcox.test(t1dat, t2dat, paired = TRUE)
  pvals.hf <- c(pvals.hf, this.p)
  means.group1 <- c(means.group1, mean(t1dat))
  means.group2 <- c(means.group2, mean(t2dat))
} 

diff.means.hf <- means.group2 - means.group1 
diff.means
data.frame(colnames(OTU.genus), diff.means.hf)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.hf), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.hfdiffmean.xlsx")


#correlations for SC at genus level

t3dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "A" & metadata$subj_id != "SUBJ4"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t3dat, t4dat, paired = TRUE)

pvals.sc <- c() 

for (i in 1:ncol(OTU.genus)) {
  t3dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "A" & metadata$subj_id != "SUBJ4"), i]
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  this.p <- wilcox.test(t3dat, t4dat, paired = TRUE)$p.value
  pvals.sc <- c(pvals.sc, this.p) 
  
} 

data.frame(colnames(OTU.genus), pvals.sc)
write.xlsx(data.frame(colnames(OTU.genus), pvals.sc), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.sc.xlsx")

#diff in means for SC at genus level
t3dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "A" & metadata$subj_id != "SUBJ4"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t3dat, t4dat, paired = TRUE)

pvals.scdiff <- c() 
means.group3 <- c() 
means.group4 <- c()
for (i in 1:ncol(OTU.genus)) {
  t3dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "A" & metadata$subj_id != "SUBJ4"), i]
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  this.p <- wilcox.test(t3dat, t4dat, paired = TRUE)$p.value
  pvals.scdiff <- c(pvals.scdiff, this.p) 
  means.group3 <- c(means.group3, mean(t3dat))
  means.group4 <- c(means.group4, mean(t4dat))
} 
diff.means.sc <- means.group4 - means.group3 
diff.means
data.frame(colnames(OTU.genus), diff.means.sc)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.sc), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.scdiffmean.xlsx")


#correlations for TV
t5dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "A" & metadata$subj_id != "SUBJ8" & metadata$subj_id != "SUBJ2"), 2]
t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), 2]
wilcox.test(t5dat, t6dat, paired = TRUE)$p.value

pvals.tv <- c() 
for (i in 1:ncol(OTU.genus)) {
  t5dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "A" & metadata$subj_id != "SUBJ8" & metadata$subj_id != "SUBJ2"), i]
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  this.p <- wilcox.test(t5dat, t6dat, paired = TRUE)$p.value
  pvals.tv <- c(pvals.tv, this.p) 
} 

data.frame(colnames(OTU.genus), pvals.tv)
write.xlsx(data.frame(colnames(OTU.genus), pvals.tv), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.tv.xlsx")

#diff in means TV
pvals.tvdiff <- c() 
means.group5 <- c() 
means.group6 <- c()
for (i in 1:ncol(OTU.genus)) {
  t5dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "A" & metadata$subj_id != "SUBJ8" & metadata$subj_id != "SUBJ2"), i]
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  this.p <- wilcox.test(t5dat, t6dat, paired = TRUE)$p.value
  pvals.tvdiff <- c(pvals.tvdiff, this.p)
  means.group5 <- c(means.group5, mean(t5dat))
  means.group6 <- c(means.group6, mean(t6dat))
} 
diff.means.tv <- means.group6 - means.group5 
diff.means
data.frame(colnames(OTU.genus), diff.means.tv)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.tv), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.tvdiffmeans.xlsx")


#correlations for LF
t7dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "A"), 2]
t8dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "C"), 2]
wilcox.test(t7dat, t8dat, paired = TRUE)$p.value


pvals.lf <- c() 
for (i in 1:ncol(OTU.genus)) {
  t7dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "A"), i]
  t8dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "C"), i]
  this.p <- wilcox.test(t7dat, t8dat, paired = TRUE)$p.value
  pvals.lf <- c(pvals.lf, this.p) 
} 

data.frame(colnames(OTU.genus), pvals.lf)
#write.csv((colnames(OTU.genus), pvals.lf))

write.xlsx(data.frame(colnames(OTU.genus), pvals.lf), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.lf.xlsx")

#diff in mean LF
t7dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "A"), 2]
t8dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "C"), 2]
wilcox.test(t7dat, t8dat, paired = TRUE)$p.value


pvals.lfdiff <- c()
means.group7 <- c() 
means.group8 <- c()
for (i in 1:ncol(OTU.genus)) {
  t7dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "A"), i]
  t8dat <- OTU.genus[which(metadata$group == "FL" & metadata$time == "C"), i]
  this.p <- wilcox.test(t7dat, t8dat, paired = TRUE)$p.value
  pvals.lfdiff <- c(pvals.lfdiff, this.p) 
  means.group7 <- c(means.group7, mean(t7dat))
  means.group8 <- c(means.group8, mean(t8dat))
} 
diff.means.lfdiff <- means.group8 - means.group7 
diff.means
data.frame(colnames(OTU.genus), diff.means.lfdiff)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.lfdiff), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.lfdiff.xlsx")

#write.csv((colnames(OTU.genus), diff.means.lfdiff))



#correlations just at time C
#SC v STEV
t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t5dat, t4dat)

pvals.tvsctimec <- c() 
for (i in 1:ncol(OTU.genus)) {
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  this.p <- wilcox.test(t6dat, t4dat)$p.value
  pvals.tvsctimec <- c(pvals.tvsctimec, this.p) 
} 

data.frame(colnames(OTU.genus), pvals.tvsctimec)
write.xlsx(data.frame(colnames(OTU.genus), pvals.tvsctimec), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.tvvSCtimeC2.xlsx")

#diff in mean 
#SC v STEV
t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t6dat, t4dat)

pvals.tvsctimec <- c() 
means.group6 <- c() 
means.group4 <- c()
for (i in 1:ncol(OTU.genus)) {
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  this.p <- wilcox.test(t5dat, t4dat)$p.value
  pvals.tvsctimec <- c(pvals.tvsctimec, this.p)
  means.group6 <- c(means.group6, mean(t6dat))
  means.group4 <- c(means.group4, mean(t4dat))
} 
diff.means.scvstev <- means.group4 - means.group6 
diff.means
data.frame(colnames(OTU.genus), diff.means.scvstev)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.scvstev), "/Users/sarahbecker/Desktop/THESIS/RDirectory/diff.means.scvstev2.xlsx")

#STEV v HF
t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), 2]
t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), 2]

wilcox.test(t2dat, t6dat)

pvals.tvvhftimec <- c() 
for (i in 1:ncol(OTU.genus)) {
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  this.p <- wilcox.test(t6dat, t2dat)$p.value
  pvals.tvvhftimec <- c(pvals.tvvhftimec, this.p) 
} 

data.frame(colnames(OTU.genus), pvals.tvvhftimec)
write.xlsx(data.frame(colnames(OTU.genus), pvals.tvvhftimec), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.tvvhftimec.xlsx")

#diff in mean
#STEV v HF
t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), 2]
t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), 2]

wilcox.test(t2dat, t6dat)

pvals.tvvhftimec <- c() 
means.group2 <- c() 
means.group6 <- c()
for (i in 1:ncol(OTU.genus)) {
  t6dat <- OTU.genus[which(metadata$group == "TV" & metadata$time == "C"), i]
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  this.p <- wilcox.test(t6dat, t2dat)$p.value
  pvals.tvvhftimec <- c(pvals.tvvhftimec, this.p) 
  means.group2 <- c(means.group2, mean(t2dat))
  means.group6 <- c(means.group6, mean(t6dat))
} 

diff.means.stevVhf <- means.group6 - means.group2 
diff.means
data.frame(colnames(OTU.genus), pvals.tvvhftimec)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.stevVhf), "/Users/sarahbecker/Desktop/THESIS/RDirectory/diff.means.stevVhf2.xlsx")

#HF v SC
t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t2dat, t4dat)

pvals.scvhf <- c() 
for (i in 1:ncol(OTU.genus)) {
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  wilcox.test(t2dat, t4dat)$p.value
  pvals.scvhf <- c(pvals.scvhf, this.p) 
} 

data.frame(colnames(OTU.genus), pvals.scvhf)

write.xlsx(data.frame(colnames(OTU.genus), pvals.scvhf), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.genus.scvhftimec2.xlsx")

#diff in mean
#HF v SC
t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), 2]
t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), 2]
wilcox.test(t2dat, t4dat)

pvals.scvhftimec <- c() 
means.group2 <- c() 
means.group4 <- c()
for (i in 1:ncol(OTU.genus)) {
  t4dat <- OTU.genus[which(metadata$group == "SC" & metadata$time == "C"), i]
  t2dat <- OTU.genus[which(metadata$group == "HF" & metadata$time == "C" & metadata$subj_id != "SUBJ8"), i]
  this.p <- wilcox.test(t4dat, t2dat)$p.value
  pvals.scvhftimec <- c(pvals.scvhftimec, this.p)
  means.group2 <- c(means.group2, mean(t2dat))
  means.group4 <- c(means.group4, mean(t4dat))
} 

diff.means.hfsc <- means.group4 - means.group2 
diff.means
data.frame(colnames(OTU.genus), diff.means.hfsc)
data.frame(colnames(OTU.genus), pvals.scvhftimec)
write.xlsx(data.frame(colnames(OTU.genus), diff.means.hfsc), "/Users/sarahbecker/Desktop/THESIS/RDirectory/diff.means.hfsc.xlsx")
write.xlsx(data.frame(colnames(OTU.genus), pvals.scvhftimec), "/Users/sarahbecker/Desktop/THESIS/RDirectory/pvals.scvhf.xlsx")


######add column of sexes and name mice
MalesNames=c("SC1A","SC2A","SC3A","SC4A","SC5A","SC1C","SC2C","SC3C","SC4C","SC5C","TV1A","TV2A","TV3A","TV4A","TV7A","TV1C","TV2C", "TV3C","TV4C","TV7C","HF1A","HF3A","HF4A","HF9A","HF10A","HF1C","HF3C","HF4C","HF9C","HF10C","FL1A","FL5A","FL7A","FL9A","FL10A","FL1C","FL5C","FL7C","FL9C","FL10C")
FemalesNames=c("SC6A","SC7A","SC8A","SC9A","SC10A","SC6C","SC7C","SC8C","SC9C","SC10C","TV5A","TV6A","TV7A","TV9A","TV10A","TV5C","TV6C", "TV9C","TV10C","HF2A","HF5A","HF6A","HF7A","HF8A","HF2C","HF5C","HF6C","HF7C","HF8C","FL2A","FL3A","FL4A","FL6A","FL8A","FL2C","FL3C","FL4C","FL6C","FL8C")
metadata$sex =rep(NA, nrow(metadata))

metadata$sex[which(metadata$sample_id%in%MalesNames)] = "Male"
metadata$sex[which(metadata$sample_id%in%FemalesNames)] = "Female"

###########Bray Curtis plot
BC.nmds = metaMDS(OTU.clean[-c(40,41),], distance="bray", k=2, trymax=1000, halfchange=TRUE,threshold=.8)
plotpoints = BC.nmds$points

#plot males
Amalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])
Cmalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Male"),], metadata[which(metadata$time=="C" & metadata$sex=="Male"),])

#####POSTER PLOT**********************************
#a plot males
jpeg("BC-males-initial-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, ylim=c(-0.5,0.25), xlim=c(-0.6,0.5), main = "Bray-Curtis: Males Initial Beta-Diversity", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique("LF", "HF", "Saccharin", "Stevia"), col=1:length(unique(Amalepoints$group)), pch=19)


dev.off()


dev.off()
#a plot males, standard axis -0.5-0.5
jpeg("BC-males-timeA-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), main = "Bray-Curtis: Males Initial Beta-Diversity", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(Amalepoints$group), col=1:length(unique(Amalepoints$group)), pch=19)

dev.off()

#c plot males, standard axis -0.5-0.5
jpeg("BC-males-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Cmalepoints$MDS1, Cmalepoints$MDS2, col=as.factor(Cmalepoints$group), pch=19, ylim=c(-0.8,0.7), xlab="MDS1", ylab="MDS2", main="Bray-Curtis: Males Beta-Diversity at End of Treatment")
legend("topleft", unique(Cmalepoints$group), col=1:length(unique(Cmalepoints$group)), pch=19, horiz = TRUE)
dev.off()

#a plot all
jpeg("BC-all-inital-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(rbind(Amalepoints,AFemalepoints)$MDS1,rbind(Amalepoints,AFemalepoints)$MDS2, col=as.factor(rbind(Amalepoints,AFemalepoints)$group), ylim=c(-0.6,0.6), xlim=c(-0.6,0.3), pch=19, main = "Bray-Curtis: All Groups Initial", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(rbind(Amalepoints,AFemalepoints)$group), col=1:length(unique(rbind(Amalepoints,AFemalepoints)$group)), pch=19)
dev.off()

#c plot all
jpeg("BC-all-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(rbind(Cmalepoints, CFemalepoints)$MDS1,rbind(Cmalepoints,CFemalepoints)$MDS2, col=as.factor(rbind(Amalepoints,AFemalepoints)$group), ylim=c(-0.6,0.8), xlim=c(-0.6,0.8), pch=19, main = "Bray-Curtis: All Groups Post-Treatment", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(rbind(Cmalepoints,CFemalepoints)$group), col=1:length(unique(rbind(Cmalepoints,CFemalepoints)$group)), pch=19)
dev.off()

#plot females
AFemalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Female"),], metadata[which(metadata$time=="A" & metadata$sex=="Female"),])
CFemalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Female"),], metadata[which(metadata$time=="C" & metadata$sex=="Female"),])

#a plot
jpeg("BC-initial-females-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(AFemalepoints$MDS1, AFemalepoints$MDS2, col=as.factor(AFemalepoints$group), pch=19, ylim=c(-0.3,0.45), xlim=c(-0.5,0.4), main = "Bray-Curtis: Females Initial", ylab = "MDS1", xlab = "MDS2")
legend("bottomright", unique(AFemalepoints$group), col=1:length(unique(AFemalepoints$group)), pch=19)
dev.off()

#c plot
jpeg("BC-females-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(CFemalepoints$MDS1, CFemalepoints$MDS2, col=as.factor(CFemalepoints$group), pch=19, xlim=c(-0.5, 0.5), ylim=c(-0.8,0.8), xlab="MDS1", ylab="MDS2", main="Bray-Curtis: Females Post-Treatment")
legend("topleft", unique(CFemalepoints$group), col=1:length(unique(CFemalepoints$group)), pch=19, horiz = TRUE)
dev.off()

#legend("topright", "Significance: p<0.005",col="red",lwd = 3)


#Calculate distance and save as a matrix
BC.dist=vegdist(OTU.clean[-c(40,41),], distance="bray")

#
######stats anlaysis on BC values
#Run PERMANOVA on distances.
adonis(BC.dist ~ group*time, data = metadata, permutations = 1000)

OTU.temp = OTU.clean[-c(40,41),]

###males
dist.BC.males = vegdist(OTU.temp[which((metadata$sex=="Male")),], distance="bray")
adonis(dist.BC.males ~ group*time, data = metadata[which(metadata$sex=="Male"),], permutations = 1000)

##females
dist.BC.females = vegdist(OTU.temp[which((metadata$sex=="Female")),], distance="bray")
adonis(dist.BC.females ~ group*time, data = metadata[which(metadata$sex=="Female"),], permutations = 1000)

#ttests between points at time C, could do a multivariate ttest or single variable
#single variable (MDS1)
plotpoints=as.data.frame(plotpoints)
t.test(plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="HF")], plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="SC")])
t.test(plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="HF")], plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="TV")])
t.test(plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="TV")], plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="SC")])
                                                                                                                                                                                                                            
t.test(plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="HF", plotpoints$MDS1[which(metadata$time=="C" & metadata$group=="SC",)

adonis(dist.BC.timeC ~ group*time, data = metadata[which(metadata$time=="C"),], permutations = 1000)

#multivariate ttest
set.seed(1)
library(ICSNP)
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF"),], plotpoints[which(metadata$time=="C" & metadata$group=="TV"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="TV"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC"),])                                                                                                                                                                                                              

#just looking at males at time C
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF" & metadata$sex=="Male"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC"& metadata$sex=="Male"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF" & metadata$sex=="Male"),], plotpoints[which(metadata$time=="C" & metadata$group=="TV" & metadata$sex=="Male"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="TV" & metadata$sex=="Male"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC" & metadata$sex=="Male"),])                                                                                                                                                                                                              

#just looking at females at time C
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF" & metadata$sex=="Female"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC"& metadata$sex=="Female"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="HF" & metadata$sex=="Female"),], plotpoints[which(metadata$time=="C" & metadata$group=="TV" & metadata$sex=="Female"),])
HotellingsT2(plotpoints[which(metadata$time=="C" & metadata$group=="TV" & metadata$sex=="Female"),], plotpoints[which(metadata$time=="C" & metadata$group=="SC" & metadata$sex=="Female"),])                                                                                                                                                                                                              



############jaccard plot
J.nmds = metaMDS(OTU[-c(40,41),], distance="jaccard", k=2, trymax=1000)

plotpoints = J.nmds$points

#plot males
Amalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])
Cmalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Male"),], metadata[which(metadata$time=="C" & metadata$sex=="Male"),])

#a plot males
jpeg("Jacc-males-timeA-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, ylim=c(-0.5,0.6), xlim=c(-1.5,1.5), main = "Jaccard: Males Initial", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(Amalepoints$group), col=1:length(unique(Amalepoints$group)), pch=19)
dev.off()

#c plot males
jpeg("Jacc-males-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Cmalepoints$MDS1, Cmalepoints$MDS2, col=as.factor(Cmalepoints$group), pch=19, xlab="MDS1", ylab="MDS2", main="Jaccard: Males Post-Treatment")
legend("topleft", unique(Cmalepoints$group), col=1:length(unique(Cmalepoints$group)), pch=19, horiz = TRUE)
dev.off()

#plot females

AFemalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Female"),], metadata[which(metadata$time=="A" & metadata$sex=="Female"),])
CFemalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Female"),], metadata[which(metadata$time=="C" & metadata$sex=="Female"),])

#a plot
jpeg("jacc-females-initial-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(AFemalepoints$MDS1, AFemalepoints$MDS2, col=as.factor(AFemalepoints$group), pch=19, xlim=c(-1.5,0.5), main = "Jaccard: Females Initial", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(AFemalepoints$group), col=1:length(unique(AFemalepoints$group)), pch=19)
dev.off()
#c plot
jpeg("Jacc-females-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(CFemalepoints$MDS1, CFemalepoints$MDS2, col=as.factor(CFemalepoints$group), pch=19, ylim=c(-0.3,0.4), xlim=c(-1,1), xlab="MDS1", ylab="MDS2", main="Jaccard: Females Post-Treatment")
legend("bottomright", unique(CFemalepoints$group), col=1:length(unique(CFemalepoints$group)), pch=19, horiz = TRUE)
dev.off()

##stats analysis of Jaccard
#calculate distance and save as matrix
OTU.temp=OTU.clean[-c(40,41),]
J.dist=vegdist(OTU.temp[-c(74),], distance="jaccard")

J.dist=vegdist(OTU.temp, distance="jaccard")

######run PERMANOVA on distance
adonis(J.dist ~ group*time, data = metadata, permutations = 1000)

####following lines dont work
####PERMANOVA Males
adonis(J.dist$points[which(metadata$sex=="Male"),] ~ group*time, data = metadata[which(metadata$sex=="Male"),], permutations = 1000)

#####PERMANOVA females
adonis(J.dist$points[which(metadata$sex=="Female"),] ~ group*time, data = metadata[which(metadata$sex=="Female"),], permutations = 1000)



###########weighted UniFrac distance
OTU.U = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
OTU.UF = OTU.U[-c(40,41),]
tax.UF = tax_table(as.matrix(tax.clean))
meta.UF = sample_data(metadata)
sample_names(meta.UF) = sample_names(OTU.UF)
physeq = phyloseq(OTU.UF, tax.UF, meta.UF)

load("~/Desktop/THESIS/RDirectory/NJ.tree.RData")
physeq.tree = merge_phyloseq(physeq, NJ.tree)
physeq.tree

wUF.ordu = ordinate(physeq.tree, method="NMDS", distance="unifrac", weighted=TRUE)
par(mfrow=c(1,1))

plot(wUF.ordu, type="n", main="Weighted UniFrac")
plotpoints = wUF.ordu$points
#plot total
jpeg("unifrac-males-initial-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), main = "Weighted UniFrac: Males Initial", ylab = "MDS1", xlab = "MDS2")
legend("bottomright", unique(Amalepoints$group), col=1:length(unique(Amalepoints$group)), pch=19)
dev.off()

#plot males
Amalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])
Cmalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Male"),], metadata[which(metadata$time=="C" & metadata$sex=="Male"),])

#AllpointsA=cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])
#AllpointsC=cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])

#a plot males
jpeg("unifrac-males-initial-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), main = "Weighted UniFrac: Males Initial", ylab = "MDS1", xlab = "MDS2")
legend("bottomright", unique(Amalepoints$group), col=1:length(unique(Amalepoints$group)), pch=19)
dev.off()

#c plot males
jpeg("unifrac-males-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(Cmalepoints$MDS1, Cmalepoints$MDS2, col=as.factor(Cmalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), xlab="MDS1", ylab="MDS2", main="Weighted UniFrac: Males Post-Treatment")
legend("topleft", unique(Cmalepoints$group), col=1:length(unique(Cmalepoints$group)), pch=19, horiz = TRUE)
dev.off()

#plot females
AFemalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Female"),], metadata[which(metadata$time=="A" & metadata$sex=="Female"),])
CFemalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Female"),], metadata[which(metadata$time=="C" & metadata$sex=="Female"),])

#a plot
jpeg("unifrac-females-initial-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(AFemalepoints$MDS1, AFemalepoints$MDS2, col=as.factor(AFemalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), main = "Weighted UniFrac: Females Initial", ylab = "MDS1", xlab = "MDS2")
legend("bottomleft", unique(AFemalepoints$group), col=1:length(unique(AFemalepoints$group)), pch=19)
dev.off()

#c plot
jpeg("unifrac-females-final-standardaxis-poster.jpeg", width = 1500, height = 800, quality = 100, res=205)
plot(CFemalepoints$MDS1, CFemalepoints$MDS2, col=as.factor(CFemalepoints$group), pch=19, ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), xlab="MDS1", ylab="MDS2", main="Weighted UniFrac: Females Post-Treatment")
legend("bottomright", unique(CFemalepoints$group), col=1:length(unique(CFemalepoints$group)), pch=19, horiz = TRUE)
dev.off()

####stats analysis of weighted UniFrac

#Calculate distance and save as a matrix
##makes metadata match OTU temp, need to remove PCR values
OTU.temp = OTU.clean[-c(40,41),]

dist.weight=vegdist(OTU.temp[-c(74),], distance="bray")

###males
dist.weight.males = vegdist(OTU.temp[which((metadata$sex=="Male")),], distance="bray")

##females
dist.weight.females = vegdist(OTU.temp[which((metadata$sex=="Female")),], distance="bray")


#Run PERMANOVA on distances.

###total PERMANOVA
adonis(dist.weight ~ group*time, data = metadata[-c(74),], permutations = 1000)

####PERMANOVA Males
adonis(dist.weight.males ~ group*time, data = metadata[which(metadata$sex=="Male"),], permutations = 1000)

#####PERMANOVA females
adonis(dist.weight.females ~ group*time, data = metadata[which(metadata$sex=="Female"),], permutations = 1000)

#permanova just at time c between sweetener groups
sweetener=OTU.temp[c(which(metadata$group=="SC" & metadata$time=="C"), which(metadata$group=="TV" & metadata$time=="C"), which(metadata$group=="HF"& metadata$time=="C")),]
timec.wUF=vegdist(sweetener, distance="bray")

adonis(timec.wUF ~ group, data=metadata[c(which(metadata$group=="SC" & metadata$time=="C"), which(metadata$group=="TV" & metadata$time=="C"), which(metadata$group=="HF"& metadata$time=="C")),], permutations =1000)




############unweighted UniFrac
uwUF.ordu = ordinate(physeq.tree, method="NMDS", distance="unifrac", weighted=FALSE)
plotpoints = uwUF.ordu$points
#plot males
Amalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Male"),], metadata[which(metadata$time=="A" & metadata$sex=="Male"),])
Cmalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Male"),], metadata[which(metadata$time=="C" & metadata$sex=="Male"),])

#a plot males
plot(Amalepoints$MDS1, Amalepoints$MDS2, col=as.factor(Amalepoints$group), pch=19, main = "Unweighted UniFrac: Males Time A", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(Amalepoints$group), col=1:length(unique(Amalepoints$group)), pch=19)

#c plot males
plot(Cmalepoints$MDS1, Cmalepoints$MDS2, col=as.factor(Cmalepoints$group), pch=19, xlab="MDS1", ylab="MDS2", main="Unweighted UniFrac: Males Time C")
legend("topleft", unique(Cmalepoints$group), col=1:length(unique(Cmalepoints$group)), pch=19, horiz = TRUE)

#plot females
AFemalepoints = cbind(plotpoints[which(metadata$time=="A" & metadata$sex=="Female"),], metadata[which(metadata$time=="A" & metadata$sex=="Female"),])
CFemalepoints=cbind(plotpoints[which(metadata$time=="C" & metadata$sex=="Female"),], metadata[which(metadata$time=="C" & metadata$sex=="Female"),])

#a plot
plot(AFemalepoints$MDS1, AFemalepoints$MDS2, col=as.factor(AFemalepoints$group), pch=19, main = "Unweighted UniFrac: Females Time A", ylab = "MDS1", xlab = "MDS2")
legend("topright", unique(AFemalepoints$group), col=1:length(unique(AFemalepoints$group)), pch=19)

#c plot
plot(CFemalepoints$MDS1, CFemalepoints$MDS2, col=as.factor(CFemalepoints$group), pch=19, xlab="MDS1", ylab="MDS2", main="Unweighted UniFrac: Females Time C")
legend("bottomright", unique(CFemalepoints$group), col=1:length(unique(CFemalepoints$group)), pch=19, horiz = TRUE)


####stats analysis of unweighted UniFrac
set.seed(1)
dist.unweight = UniFrac(physeq.tree, weighted=FALSE, normalized=TRUE)
adonis(dist.unweight ~ group*time, data=metadata, permutations = 1000)

###PERMANOVA males
adonis(uwUF.ordu$points[which(metadata$sex=="Male"),] ~ group*time, data = metadata[which(metadata$sex=="Male"),], permutations = 1000)

adonis(uwUF.ordu$points[which(metadata$sex=="Male" & metadata$time=="C"),] ~group*time, data=metadata[which(metadata$sex=="Male" & metadata$time=="C"),], permutations = 1000)
#####PERMANOVA females
adonis(uwUF.ordu$points[which(metadata$sex=="Female"),] ~ group*time, data = metadata[which(metadata$sex=="Female"),], permutations = 1000)




###########alpha diversity
library(vegan)

## shannon for each sample -- I'll just add it to the metadata
metadata$shannon <- diversity(OTU.genus, index = "shannon")
metadata$subj_id <- paste(metadata$subj_id, metadata$group, sep = "_")


metadata.wide <- data.frame(subj_id = unique(metadata$subj_id),
                            stringsAsFactors = FALSE)
metadata.wide$group <- sapply(metadata.wide$subj_id, FUN = function(x) strsplit(x, split = "_")[[1]][2])
metadata.wide$shannon.time1 = metadata.wide$shannon.time2 = NA
for (i in 1:nrow(metadata.wide)) {
  if (any(metadata$subj_id == metadata.wide$subj_id[i] & metadata$time == "A")) {
    st1 <- metadata$shannon[which(metadata$subj_id == metadata.wide$subj_id[i] & metadata$time == "A")]
  } else {
    st1 <- NA
  }
  
  if (any(metadata$subj_id == metadata.wide$subj_id[i] & metadata$time == "C")) {
    st2 <- metadata$shannon[which(metadata$subj_id == metadata.wide$subj_id[i] & metadata$time == "C")]
  } else {
    st2 <- NA
  }
  
  metadata.wide[i, "shannon.time1"] <- st1
  metadata.wide[i, "shannon.time2"] <- st2
}
head(metadata.wide)

## t.tests across time within treatment
with(metadata.wide[metadata.wide$group == "FL", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "HF", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "SC", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "TV", ], t.test(shannon.time1, shannon.time2, paired = TRUE))

#####t test between groups
t.test(metadata.wide$shannon.time1[metadata.wide$group=="FL"], metadata.wide$shannon.time1[metadata.wide$group=="HF"])


## anova on differences
metadata.wide$shannon.diff <- metadata.wide$shannon.time2 - metadata.wide$shannon.time1
anova(lm(metadata.wide$shannon.diff ~ metadata.wide$group))


## box plots by group
jpeg("shannon-boxplot-poster3.jpeg", width = 1500, height = 800, quality = 100, res=150)
boxplot(metadata$shannon ~  metadata$time + metadata$group,
        names = rep(c("Pre","Post"), 4),
        at = c(1:2, 4:5, 7:8, 10:11),
        ylab = "Shannon Index", ylim = c(1.6, 2.8))
abline( v = c(3, 6, 9), lty = 2 )
mtext(c("LF", "HF", "Saccharin", "Stevia"), side = 1, line = 2.5,
      at = c(1.5, 4.5, 7.5, 10.5))
mtext("Treatment", side = 1, line = 3.5)
text("*", x = 7.5, y = 2.8)
segments(x0 = 6.75, x1 = 8.25, y0 = 2.75, y1 = 2.75)
text("**", x = 4.5, y = 2.8)
segments(x0 = 3.75, x1 = 5.25, y0 = 2.75, y1 = 2.75)
text("***", x = 10.5, y = 2.8)
segments(x0 = 9.75, x1 = 11.25, y0 = 2.75, y1 = 2.75)
legend(10,2, c("* p < 0.05", "** p < 0.01", "*** p < 0.001"), cex = 1,
       border = NA, bty = 'n', xjust = 0)
dev.off()

mod <- lm(metadata$shannon ~ metadata$time*metadata$group)
coef(summary(mod))


#beeswarm plot
library(beeswarm)
jpeg("shannon-beeswarm-poster.jpeg", width = 1500, height = 800, quality = 100, res=150)
beeswarm(metadata$shannon ~ factor(metadata$group) + metadata$time, 
         at = c(1:4, 6:9), pch = 1, cex = 0.75, 
         ylab = "Shannon Index", ylim = c(1.6, 2.9), 
         labels = rep(c("LF", "HF", "Saccharin", "Stevia"), 2), 
         xlab = "Treatment Group", main = "Shannon Index by Time and Treatment")
abline(v = 5, lty = 2)
sh.means <- aggregate(metadata$shannon, by = list(metadata$group, metadata$time), FUN = mean)$x
points(sh.means ~ c(1:4, 6:9), pch = 15, cex = 1.25)
mtext(c("Pre-treatment", "Post-treatment"), side = 3, line = -1.5, at = c(2.5, 7.5))
dev.off()


## t.tests across time within treatment
with(metadata.wide[metadata.wide$group == "FL", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "HF", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "SC", ], t.test(shannon.time1, shannon.time2, paired = TRUE))
with(metadata.wide[metadata.wide$group == "TV", ], t.test(shannon.time1, shannon.time2, paired = TRUE))

#####t test between groups at time1
t.test(metadata.wide$shannon.time1[metadata.wide$group=="FL"], metadata.wide$shannon.time1[metadata.wide$group=="HF"])
t.test(metadata.wide$shannon.time1[metadata.wide$group=="FL"], metadata.wide$shannon.time1[metadata.wide$group=="TV"])
t.test(metadata.wide$shannon.time1[metadata.wide$group=="FL"], metadata.wide$shannon.time1[metadata.wide$group=="SC"])
t.test(metadata.wide$shannon.time1[metadata.wide$group=="FL"], metadata.wide$shannon.time1[metadata.wide$group=="HF"])

#####t test between groups at time2
t.test(metadata.wide$shannon.time2[metadata.wide$group=="FL"], metadata.wide$shannon.time2[metadata.wide$group=="HF"])
t.test(metadata.wide$shannon.time2[metadata.wide$group=="FL"], metadata.wide$shannon.time2[metadata.wide$group=="TV"])
t.test(metadata.wide$shannon.time2[metadata.wide$group=="FL"], metadata.wide$shannon.time2[metadata.wide$group=="SC"])
t.test(metadata.wide$shannon.time2[metadata.wide$group=="FL"], metadata.wide$shannon.time2[metadata.wide$group=="HF"])


#CHAO
## chao for each sample
#chao1.out <- estimateR(OTU.genus.orig)
#metadata$chao1 <- chao1.out[1,]

## box plots by group
#jpeg("chao-boxplot-.jpeg", width = 1500, height = 800, quality = 100, res=150)

#boxplot(metadata$chao1 ~  metadata$time + metadata$group,
      #  names = rep(c("Pre", "Post"), 4),
      #  at = c(1:2, 4:5, 7:8, 10:11),
      #  ylab = "Species Richness (Chao1)")
#xlab = "Treatment Group", main = "Shannon Index by Time and Treatment")
#abline( v = c(3, 6, 9), lty = 2 )
#mtext(c("LF", "HF", "Saccharin", "Stevia"), side = 1, line = 2.5,
     # at = c(1.5, 4.5, 7.5, 10.5))
#mtext("Treatment", side = 1, line = 3.5)
#legend(10,2, c("* p < 0.05", "** p < 0.01", "*** p < 0.001"), cex = 1,
       #border = NA, bty = 'n', xjust = 0)

#dev.off()
## chao for each sample
#chao1.out <- estimateR(OTU.genus)
#metadata$chao <- chao1.out[2,]



#metadata$subj_id2 <- paste(metadata$subj_id, metadata$group, metadata$group, sep = "_")


#metadata.wide <- data.frame(subj_id = unique(metadata$subj_id),
#stringsAsFactors = FALSE)
#metadata.wide$group <- sapply(metadata.wide$subj_id, FUN = function(x) strsplit(x, split = "_")[[1]][2])
#metadata.wide$chao.time1 = metadata.wide$chao.time2 = NA
#for (i in 1:nrow(metadata.wide)) {
  #if (any(metadata$subj_id2 == metadata.wide$subj_id[i] & metadata$time == "A")) {
  #  st1 <- metadata$chao[which(metadata$subj_id2 == metadata.wide$subj_id[i] & metadata$time == "A")]
  #   st1 <- NA
 # }
 ## 
  #if (any(metadata$subj_id2 == metadata.wide$subj_id[i] & metadata$time == "C")) {
 #   st2 <- metadata$chao[which(metadata$subj_id2 == metadata.wide$subj_id[i] & metadata$time == "C")]
 # } else {
 #   st2 <- NA
 # }
  
  metadata.wide[i, "chao.time1"] <- st1
  metadata.wide[i, "chao.time2"] <- st2
}
#head(metadata.wide)


#metadata$chao1 = rep(NA, NROW(metadata))
#metadata$chao1[metadata$time=="A"] = metadata$chao[metadata$time=="A"]


#metadata$chao2 = rep(NA, NROW(metadata))
#metadata$chao2[metadata$time=="C"] = metadata$chao[metadata$time=="C"]

#metadata.wide2 <- reshape()



#ttests not working, help with code to test signficance between groups
## t.tests across time within treatment
#chao.time1 = metadata$chao1 & metadata$time=="A"]
#chao.time2= (metadata$time=="C"),]

#with(metadata.wide[metadata$group == "FL", ], t.test(chao1.time1, chao1.time2, paired = TRUE))
#with(metadata.wide[metadata$group == "HF", ], t.test(chao.time1, chao.time2, paired = TRUE))
#with(metadata.wide[metadata$group == "SC", ], t.test(chao.time1, chao.time2, paired = TRUE))
#with(metadata.wide[metadata$group == "TV", ], t.test(chao.time1, chao.time2, paired = TRUE))

## anova on differences
#metadata.wide$shannon.diff <- metadata.wide$shannon.time2 - metadata.wide$shannon.time1
#anova(lm(metadata.wide$shannon.diff ~ metadata.wide$group))


## box plots by group

#boxplot(metadata$chao1 ~  metadata$time + metadata$group,
       # names = rep(c("T1", "T2"), 4),
       # at = c(1:2, 4:5, 7:8, 10:11),
       # ylab = "Species Richness (Chao1)")
#abline( v = c(3, 6, 9), lty = 2 )
#mtext(c("FL", "HF", "SC", "TV"), side = 1, line = 2.5,
 #     at = c(1.5, 4.5, 7.5, 10.5))
#mtext("Treatment", side = 1, line = 3.5)
#mymodelFL=aov(chao1~time,data=metadata[metadata$group=="FL",])
#summary(mymodelFL)

#mymodelHF=aov(chao1~time,data=metadata[metadata$group=="HF",])
#summary(mymodelHF)

#mymodel=aov(chao1~group,data=metadata[metadata$time=="A",])
#summary(mymodel)

#mymodelSC=aov(chao1~time,data=metadata[metadata$group=="SC",])
#summary(mymodelSC)

#mymodelTV=aov(chao1~time,data=metadata[metadata$group=="TV",])
#summary(mymodelTV)

#mymodelC=aov(chao1~group,data=metadata[metadata$time=="C",])
#summary(mymodelC)

