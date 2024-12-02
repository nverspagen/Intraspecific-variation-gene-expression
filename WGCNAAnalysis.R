### Coexpression analysis for Nadja, Hannu , July 2023
sinteractive --account rlehtone --mem 8000 --time 06:00:00 --tmp 1
module load r-env
start-r


rm(list=ls())
setwd("/scratch/rlehtone/hannu/nadja_cinxia/wgcna")
library("WGCNA")
library("DESeq2")
library("vsn")
library("dendextend")
library("gtools")
library("vegan")
library(ggplot2) #To make plots
library(ggbeeswarm) #for scatter violin plots
library(ggpubr) #To organise ggplots in a grid
library(dplyr)
library(tidyr)

#Load the summarized experiment
cinxia_counts <- readRDS("smrExptOrg.rds")

#Select the counts column and round to 0 digits
count_table <- as.data.frame((cinxia_counts@assays@data@listData[[2]]))
count_table <- round(count_table, digits = 0)
# 96 samples and 14810 genes

#For WGCNA, the count table needs to be transformed
count_table <- t(count_table)

#As recommended by the WGCNA manual, check if there are missing samples or genes with 0 variance
gsg <- goodSamplesGenes(count_table, verbose=3)
gsg$allOK

#Remove the genes with missing samples or 0 variance 
if (!gsg$allOK){
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dds_stab)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
   count_table = count_table[gsg$goodSamples, gsg$goodGenes]
}

#Excluding 313 genes from the calculation due to too many missing samples or zero variance

#Transform the count table back to its original form
count_table <- as.data.frame(t(count_table))

#Read the data matrix
experiment_data = as.data.frame(cinxia_counts@colData)
experiment_data$Population <- as.factor(experiment_data$Population)
experiment_data$Temperature <- as.factor(experiment_data$Temperature)

#Merge count table and data matrix (only columns 1-4 from data matrix are needed)
count_table_m <- merge(t(count_table), experiment_data[,1:4], by.x = "row.names", by.y = "ID")
row.names(count_table_m) <- count_table_m$Row.names
count_table_m <- count_table_m[2:14498]
count_table_m <- t(count_table_m)

#get the data in the right format for the dds object
count_table_deseq <- data.frame(gene = row.names(count_table_m), count_table_m)

#Make DESeq dataset object
dds <- DESeqDataSetFromMatrix(countData= count_table_deseq, 
                              colData=experiment_data[,1:4], 
                              design=~Population + Temperature, tidy = TRUE)
                              

#Try different ways of data normalisation
data_sets <- list()
data_sets[[1]] = normTransform(dds)
data_sets[[2]] = varianceStabilizingTransformation(dds, blind = T, fitType = "mean")
data_sets[[3]] = varianceStabilizingTransformation(dds, blind = F, fitType = "mean")
data_sets[[4]] = rlog(dds, blind = T)
data_sets[[5]] = rlog(dds, blind = F)

#Plot the different normalisation methods and save
plot_names <- c("log2+1","vst_blind","vst_notblind","rlog_blind","rlog_notblind")
png(file = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/FigureS1.png", width = 18, height = 13, units = "in", res = 300)
par(mfrow=c(1,5))
for(i in 1:length(data_sets)){
a <- meanSdPlot(assay(data_sets[[i]]), plot = F)
plot(x = a$gg$data$px, y = a$gg$data$py, pch = 16, cex = 1, col = "blue", main = plot_names[i], xlab = "rank-counts", ylab = "SD", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(a$rank, a$sd, col = "red", lwd = 2)
}
dev.off()

## rlog either blind or blinds seems get best results. We use blind rlog for wgcna and data exploration
dds_stab <- data_sets[[4]]@assays@data@listData[[1]] 

#Transform the data to make a sample tree to find outlier samples
dds_stab <- t(dds_stab)

#Make the sample tree
sampleTree = hclust(dist(dds_stab), method = "average")
dend <- as.dendrogram(sampleTree)

#Make labels for the samples based on temperature and population
pop = rep("0", length(labels(dend)))
pop[grep("AF.*25", labels(dend))] = "cyan"
pop[grep("AF.*34", labels(dend))] = "darkorange"
pop[grep("CAT.*25", labels(dend))] = "cyan4"
pop[grep("CAT.*34", labels(dend))] = "darkorange4"

#Set the label colours
labels_colors(dend) <- pop

#Plot and save the sample dendrogram
png(file = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/FigureS2.png", width = 18, height = 13, units = "in", res = 300)
plot(hang.dendrogram(dend, hang = 0.1), sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
legend("topright", inset=.05, title="", c(expression("Finland, 25"*~degree*C), expression("Finland, 34"*~degree*C),expression("Spain, 25"*~degree*C), expression("Spain, 34"*~degree*C)), fill= c("cyan","darkorange","cyan4","darkorange4"), horiz=F, cex = 1.5, bty = "n")
dev.off()

#AF_142_34_117b_1 is an outlier in rlog normalized data. We remove the sample

#Transform data back
dds_stab <- t(dds_stab)

#Remove the outlier
dds_stab <- dds_stab[,grep("AF_142_34_117b_1", colnames(dds_stab), invert = T)]

#Transform data for wgcna (genes as columns)
dds_stab <- t(dds_stab)

#2.a Automatic network construction and module detection
#2.a.1 Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(dds_stab, powerVector = powers, verbose = 5, networkType = "signed")

# Plot and save the results
png(file = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/FigureS4.png", width = 18, height = 13, units = "in", res = 300)
par(mfrow = c(1,2))
cex1 = 1.5

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="p",
     pch=19, col="grey",
     ylim=c(-1,1), xaxt = "n", yaxt = "n", cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.5)
axis(side = 1, at=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20))
axis(side = 2, at=c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
title("A", adj = 0)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.80,col="black") # r2 needs to be higher than 0.8

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     pch=19, col="lightgrey",type="p",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",
     cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.5, xaxt = "n", yaxt = "n")
axis(side = 1, at=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20))
axis(side = 2, at=c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500))
title("B", adj = 0)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=100, col="black") #should be more than 100
dev.off()

# WGCNA’s authors recommend using a power that has an signed R2 above 0.80, otherwise they warn your results may be too noisy to be meaningful. If you have multiple power values with signed R2 above 0.80, then picking the one at an inflection point, in other words where the R2 values seem to have reached their saturation (Zhang and Horvath 2005). You want to a power that gives you a big enough R2 but is not excessively large. Mean connectivity should be (in the hundreds or above).

#4 satisfies both these criteria, so that is what we pick as a soft-thresholding power

#2.a.2 One-step network construction and module detection
cor <- WGCNA::cor

#Make the network
net_4 = blockwiseModules(dds_stab, power = 4,
                       TOMType = "signed", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       maxBlockSize = 20000,
                       networkType = "signed",
                       saveTOMFileBase = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/Cinxia_power_4.net",
                       verbose = 3,
                       randomSeed=1234,
                       deepSplit=4) #deepSplit: Provides a simplified control over how sensitive module detection 

#Assign colours to each module
net_4_colors <- labels2colors(net_4$colors)

#Show the number of modules and genes per module
table(net_4_colors)

# Calculate eigengenes (PC1)
MEList = moduleEigengenes(dds_stab, colors = net_4_colors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result

MEDissThres = 0.25

png(file = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/ME_tree.png", width = 22, height = 18, units = "in", res = 300)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function to merge similar modules
merged_mods = mergeCloseModules(dds_stab, net_4_colors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merged_mods$colors;

table(mergedColors)

#Add them to a dataframe and change the column name
MergedModules_genes <- as.data.frame(mergedColors)
colnames(MergedModules_genes) <- "Module"

#Add the genes
MergedModules_genes$Gene <- count_table_deseq$gene

#Save
write.csv(MergedModules_genes, "MergedModules_genes.csv")

# Eigengenes of the new merged modules:
mergedMEs = merged_mods$newMEs;

png(file = "/scratch/rlehtone/hannu/nadja_cinxia/wgcna/FigureS5.png", width = 22, height = 18, units = "in", res = 300)
plotDendroAndColors(net_4$dendrograms[[1]], cbind(net_4_colors, mergedColors),
                  c("Original modules", "Merged modules"),
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  autoColorHeight = FALSE, colorHeight = 0.2,
                  cex.colorLabels = 1.5, cex.lab=1.5, cex.axis=1.5,
                  main = "", marAll = c(0, 7, 0.5, 0))
dev.off()


##### ANOVA analysis on Modulule eigengenes (PC1) #######

# Save merged module eigengenes
MEs = mergedMEs

#Add experimental data
MEs_experiment = merge(MEs, experiment_data, by.x = "row.names", by.y = "ID")

#Remove ME from the column names and remove the grey module
me_for_anova = grep("ME",colnames(MEs_experiment), value = T)
me_for_anova = gsub("ME","",me_for_anova)
me_for_anova = me_for_anova[-which(me_for_anova =="grey")]
colnames(MEs_experiment) <- gsub("ME","", colnames(MEs_experiment))
MEs_experiment <- MEs_experiment[,!colnames(MEs_experiment)=="grey"]

#Save the merged eigengenes
write.csv(MEs_experiment, "WGCNA_eigengenes_merged.csv")

#Create a table to store the p values (column 1 = population, column 2 = temperature, column 3 = interaction, a row for each module) and variance explained (R2, column 1 = population, column 2 = temperature, column 3 = interaction, column 4 = residual, a row for each module)
anova_perm = list()
anova_perm[["p-value"]] = matrix(ncol=3, nrow = length(me_for_anova))
row.names(anova_perm[["p-value"]]) = me_for_anova
anova_perm[["variance"]] = matrix(ncol=4, nrow = length(me_for_anova))
row.names(anova_perm[["variance"]]) = me_for_anova

#Perform anova permutation using the vegan package for each module and store the p values in the anova_perm matrix
set.seed(25)
for(i in 1:length(me_for_anova)){
  a = adonis2(MEs_experiment[,me_for_anova[i]] ~ as.factor(Population) * as.factor(Temperature), method = "euclidean", permutations = 999, data = MEs_experiment)
  anova_perm[["p-value"]][i,] = a[,5][1:3]
  anova_perm[["variance"]][i,] = a[,3][1:4]

}

#Make a table for the adjusted p values, corrected for multiple testing (column 1 = population, column 2 = temperature, column 3 = interaction, a row for each module)
p_adjusted = list()
p_adjusted[["p-adjust"]] = matrix(ncol=3, nrow = length(me_for_anova))
row.names(p_adjusted[["p-adjust"]]) = me_for_anova
colnames(p_adjusted[["p-adjust"]]) = c("Population", "Temperature", "Interaction")

#Adjust p values
for(i in 1:ncol(p_adjusted[["p-adjust"]])){
  bb = unname(p.adjust(anova_perm[["p-value"]][,i], method = "BH"))
  p_adjusted[["p-adjust"]][,i] = bb
}

#Sort the table based on p values
p_adjusted[["p-adjust"]] = p_adjusted[["p-adjust"]][order(p_adjusted[["p-adjust"]][,3],p_adjusted[["p-adjust"]][,1],p_adjusted[["p-adjust"]][,2]),]

#Save adjusted p values
write.csv(p_adjusted, "Adjusted p values.csv")

#Plotting the data

#Set general ggplot theme
theme_set(theme_bw() +
            theme(aspect.ratio=1,
                  axis.title = element_text(size = 14, colour = "black"),
                  axis.text = element_text(size = 12, colour = "black", margin = margin(r = 0)),
                  plot.title = element_text(size = 12, colour = "black", hjust = 0.5),
                  strip.text.x = element_text(size = 12, colour = "black"),
                  legend.text = element_text(size = 12, colour = "black"),
                  legend.key.size = unit(0.7, "cm"),
                  legend.key.width = unit(1,"cm"),
                  legend.position="bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.3),
                  axis.ticks = element_line(size = 0.3),
                  axis.ticks.length = unit(0.1, "cm")
            )
)

#Set dodge width for all plots
d = 0.5

#Rename the merged module columns to include the number of genes
wgcna <- MEs_experiment %>% 
  rename(
    "Black (2673)" = black,
    "Blue (1389)" = blue,
    "Brown (1712)" = brown,
    "Cyan (278)" = cyan,
    "Green (1364)" = green,
    "Greenyellow (247)" = greenyellow,
    "Lightcyan (75)" = lightcyan,
    "Pink (870)" = pink,
    "Salmon (180)" = salmon,
    "Tan (200)" = tan,
    "Yellow (1243)" = yellow
  )

#from wide to long
wgcna <- gather(wgcna, "Module", "Eigengene", 2:12)

#Set the modules in the right order for plotting (first only temperature significant, then only population significant, then interaction significant, temperature + population, temperature + interaction, all, none)
wgcna$Module <- factor(wgcna$Module, levels = c("Black (2673)", "Pink (870)", "Yellow (1243)", "Salmon (180)", "Brown (1712)", "Green (1364)", "Tan (200)", "Blue (1389)", "Cyan (278)", "Greenyellow (247)", "Lightcyan (75)"))

#Set temperature as a factor
wgcna$Temperature <- as.factor(wgcna$Temperature)

#Check structure
str(wgcna)

#Calculate the mean eigengene value for each temperature and population in each module
wgcna_mean <- wgcna %>% group_by(Temperature, Population, Module) %>%
  summarize(Eigengene = mean(Eigengene))

#plot
Figure2 <- ggplot(wgcna_mean, aes(Temperature, Eigengene, group = Population, col = Population, linetype = Population, shape = Population)) +
  geom_quasirandom(data = wgcna, aes(Temperature, Eigengene), 
                   width = 0.1, size = 3, alpha=0.2, dodge.width = 0.5) +
  geom_point(position = position_dodge(d), size = 3) +
  geom_line(position = position_dodge(d), size = 1) +
  facet_wrap(~Module, ncol = 4) +
  labs(x = "", y = "Eigengene") +
  scale_color_manual(values = c("#0f5b05", "#c76800")) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks=seq(-0.3, 0.3, 0.1), labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.position = "none")

Figure2

#Save the pre-processing version of figure 2. After processing, the panels will be colour-coded based on significance and a legend will be added
ggsave("Figure2_preprocessing.png", plot = Figure2, width = 25, height = 23, units = "cm", bg = "white", dpi = 300)
