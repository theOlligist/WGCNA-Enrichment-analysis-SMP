rm(list = ls())
source("~/Dropbox/R_general_working_directory/Dictionary_of_Colors.R")

library(tidyverse); library(edgeR); library(reshape2); library(WGCNA); library(patchwork);library(jishonoiro);library(vegan);library(corrplot)
setwd("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/")
#load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/WGCNA_assemply2018.RData")
save.image("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/WGCNA_assemply2018_final_V1.RData")

one41= c("#FF5200", "#0D2B52","#A6FF47")
pal1 = c("#8c6510", "#0D2B52","#FF7399")
color_pal = c("black", "light grey", "steel blue")

# Load and process the raw count table ------------------------------------------
#filter to the dominant diatom taxa from 2018 with enough reads to meaningfully do so (Figure S#-diatoms). I left off Rhizosolenia because there just weren't enough reads.
Diatom_raw_df = read.delim("raw-count-data-metaT-SMP2018.txt") %>% 
  filter(str_detect(Taxonomy, "Thalassiosirales") | str_detect(Taxonomy, "Bacillariales")) %>% 
  mutate(tax = case_when(str_detect(Taxonomy, "Thalassiosirales") ~"Thalassiosira", 
                         str_detect(Taxonomy, "Bacillariales") ~"Pseudo_nitz")) %>% 
  select(tax, KO, everything()) %>% 
  unite("ID", c(tax,KO), sep = "-")

Diatom_raw_df %>% 
  filter(str_detect(ID, "Thalas")) %>% nrow()
nrow(Diatom_raw_df) # there are approx 26,000 transcripts
#There's only a little over 600 genes for rizosolenia but I'll roll with it.

#And count the number of KO terms. I plan to use the KO IDs as the rownames;
#therefore I'll will aggregate them for each sample. Flattening the tax diversity to FX diversity.
#KO_total = Diatom_raw_df %>% distinct(KO) %>% nrow() #count the number of distinct KO
# I can expect a max o ~5622

#because I have filtered down this dataset to be diatoms, The subject of this dataframe is the functions themselves, not the taxa.
#However, there is duplication in KO terms because this is still a mixed population of diatoms expressing similar genes in some cases.
# I will aggregate the counts from each KO in each sample. This is an important step in order to assign the KO's the rownames.

#Aggregate Table
agg_df = Diatom_raw_df %>%
  pivot_longer(cols = starts_with("SMP"), names_to = "sample", values_to = "value") %>% 
  #dcast(ID~sample) %>% column_to_rownames("ID") %>% colSums()
  mutate(sample = str_replace(sample, "SMPier\\.", ""),
         sample = str_replace(sample, regex("\\.\\d{2}\\.\\d{2}\\.2018_S\\d+"), ""),
         sample = factor(sample, levels = c("Day1.A", "Day1.B", "Day1.C", 
                                            "Day3.A", "Day3.B", "Day3.C", 
                                            "Day5.A", "Day5.B", "Day5.C", 
                                            "Day9.A", "Day9.B", "Day9.C", 
                                            "Day11.A", "Day11.B", "Day11.C"))) %>% 
  group_by(ID, sample) %>%
  summarise(count = sum(value, na.rm = TRUE)) %>%
  dcast(ID~sample)

#Assign rownames as KO
rownames(agg_df) = agg_df$ID

#Take a peek at the first three rows and first six columns
agg_df[1:3, 1:16]
#sanity check. These should be the same
colSums(agg_df[-1]) == colSums(Diatom_raw_df[-1:-2])

# Normalize via rarify ----------------------------------------------------
numbers_only = agg_df[-1]
min_reads = min(colSums(keep))
rare = rrarefy(t(round(keep)), sub)
Norm_df = as.data.frame(t(rare))
colSums(Norm_df)




# # Normalize the agg_df dataframe with edgeR -----------------------------
#(1) create a dge object and (2) calculate normalizing factors and then (3) return the cpm as a dataframe
#Ingredients for the dge object: dataframe (counts only!), sample_list, & metadata.

#Make sure these line up with the order of the columns of the agg_df
# Create a sample list 
sample_list = c("Day1","Day3","Day5","Day9","Day11")
# Create the metadata
metadata = factor(c(rep("day1",3),
                    rep("day3",3),
                    rep("day5",3),
                    rep("day9",3),
                    rep("day11",3)),
                  levels = str_to_lower(sample_list)) # levels is important for purposefully setting the sample order.


# Use this information to create the DGEList object
dge_obj = DGEList(counts = as.data.frame.matrix(agg_df[-1]), # Counts is the actual columns of the dataframe
                  genes = agg_df[1], # This is the name of the column(similar to genes)
                  group = metadata)

#calculate the normalizing factors using TMM
dge_obj = calcNormFactors(dge_obj, method = "TMM")

#generate (rounded) normalized count table for use with WGCNA
Norm_df = cpm(dge_obj) %>% #return the counts
  data.frame() %>% #convert to dataframe
  round() #round the numbers because they need to be integers for WGCNA

Norm_df[1:3, 1:15] # take a peek at the normalized, rounded count table that will be used for all downstream WGCNA analysis.
Norm_df.t = t(Norm_df) # flip it on its side. several steps require the rows to be the samples and the columns to be "genes" and a matrix. 
# t() automatically does both.
nrow(Norm_df) #~7000 genes

#barplot of abundances
#### Calculate dissimilarity between samples to get an idea of the grouping of samples based on their proportions of KOs make a tree 
sample_tree = hclust(dist(Norm_df.t), #cluster using euclidean distances from the dist function
                     method = "average") #form clusters of samples use the average method.
#Plot the results
plot(sample_tree)

#There is also flashClust
flash_tree = flashClust::flashClust(dist(Norm_df.t), #cluster using euclidean distances from the dist function
                       method = "average")
plot(flash_tree)
######_______WGCNA analysis: the manual way_________#####

# Find soft threshold -----------------------------------------------------
# Soft thresholds for which to establish links between edges. 
# Our goal here is to chose the value for power that results in an optimal soft threshold.
# explain soft and hard thresholds for tutorial (undone)

#create a vector of power values, 1-20, to test out.
powers = c(c(1:10), seq(12,20,2)) 

#Use the pick softthreshold function to model the 'fit' for the different thresholds.
sft = pickSoftThreshold(Norm_df.t, powerVector = powers)

# The function returns what it believes to be the optimal power; 
sug_power = sft$powerEstimate #its not necessary to do this because the variable 9 is already saved in sft$powerEst, but it makes it more explicit and clear to see whats happening.

#however, we can plot to see if there is an alternative.
plt1 = sft$fitIndices %>%
  ggplot(., aes(x = Power, y = -sign(slope) * SFT.R.sq, label = Power)) +
  geom_point(size = 0.0, alpha = 0) +
  geom_text() +
  labs(y = "Scale Free Topology Model Fit, signed R^2",
       x = "Soft threshold (power)") +
  #geom_hline(yintercept = .90, lwd = .2, color = "red", linetype = "dotdash", alpha = .9) +
  theme_minimal()

plt2 = sft$fitIndices %>%
  ggplot(., aes(x = Power, y = mean.k., label = Power)) +
  geom_point(size = 0.0, alpha = 0) +
  geom_text() +
  labs(y = "Mean Connectivity",
       x = "Soft threshold (power)") +
  theme_minimal()
#patchwork function / operator
plt1 + plt2

# from the above plot, 10 may be a good alternative power that produces a better soft threshold (the point before the line plateaus).
alt_power = 20



# Construct gene network: Expresssion to edges -------------------------------------------
#1. Construct Adjacency Matrix: convert gene expression to connection strength
adj_mat = adjacency(Norm_df.t, power = sug_power)
k = softConnectivity(datExpr = Norm_df.t, power = sug_power) #may or not be used value of k. softconnectivity() used for more than 5000 genes.

KO_vector = rownames(adj_mat)
#2. Convert the adjacency matrix into a topological overlap matrix; 
# The authors of WGCNA (PhD mathematicians/ biostatisticians) state that converting AM to TOM minimize noise and spurrius connections.
SMP_TOM = TOMsimilarity(adj_mat)
dissTom = 1-SMP_TOM

# Costruct modules: clusters of genes with similar expression profiles
#1. Conduct heirarchical clustering and visualize the connectedness of the genes
gene_tree = hclust(as.dist(dissTom), method = "average")
sizeGrWindow(12, 9) #set the size of the canvas for this type of plot helps to center it properly
#Plot the tree
plot(gene_tree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# Branches that form dense clusters are highly coexpressed genes (modules). Module identification amounts to trimming off
# and retaining the densely clustered genes from the tree. 
# trim some branches using the dynamictreecut function


# Characterize modules ----------------------------------------------------
#There are lots of ways to trim a gene tree: cutreeDynamic() is what I'll use.
minModuleSize = 10 #set the minimum module size to 50 because I like large modules.

#we didn't set a cutoff based on the tree. Instead, i'll go with the deepSplit option
# Its an automated way to decide how conservatively to split up the tree into modules.
dynamic_Mods = cutreeDynamic(dendro = gene_tree,
                             distM = dissTom, 
                             deepSplit = 2,#1(lump) - 4(split)
                             #cutHeight = .55,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
#dynamic_Mods is a vector that contains the module membership of each gene.
table(dynamic_Mods) #use the table function to calculate the size of each module

#convert labels to colors for plotting.
dynamic_colors = labels2colors(dynamic_Mods)
table(dynamic_colors)

#plot the dendrogram as before with their module membership below
plotDendroAndColors(gene_tree,
                    dynamic_colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Thalassiosira & Pseudo-nitzschia gene expression:\n cluster dendrogram & module color",
                    cex.main = 1)

# Module curation: merging similar modules using module eigengenes--------------------------------

# Can we merge modules with similar expression profiles?
# We use eigengenes of each module to quantify co-expression similarity between them. 
# (1) calculate the eigengenes of each module; (2) calcluate /plot their dissimilarity; (3) make merge decisions.

# 1. get the eigengenes
ME_list = moduleEigengenes(Norm_df.t, 
                           colors = dynamic_colors)
ME_s = ME_list$eigengenes
#ME_s = removeGreyME(ME_s, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))

# 2. Conduct heirarchical clustering of eigengenes
# Calcultate dissimilarity (measure distance) between eigengenes
ME_Diss = 1-cor(ME_s)

#cluster the module eigengenes
ME_tree = hclust(as.dist(ME_Diss), method = "average")
dev.off()
plot(ME_tree,
     xlab = '',
     sub = "",
     cex = .8)

#For these data, an appropriate cutoff may be at 0.42: merging the modules below the threshold
merge_thresh = 0.20
abline(h = merge_thresh,
       col = "red")
#looks good, but another strategy could have been to raise the bar to 82.

# 3. Merge modules below the threshold using mergeCloseModules()
merged = mergeCloseModules(Norm_df.t,
                           dynamic_colors,
                           cutHeight = merge_thresh,
                           verbose = 3,
                           relabel = FALSE)

#Save the eigengenes of the merged modules
merged_MEs = merged$newMEs

# Get the new merged modules (colors)
merged_colors = merged$colors

#Apply KOs to the modules
names(merged_colors) = rownames(Norm_df)

# compare the old and new; the color assignment and sizes; note what was merged and how.
table(merged_colors) #merged modules
table(dynamic_colors) #original modules from dynamic trim

# visually compare module assignments
sizeGrWindow(12, 9)
plotDendroAndColors(gene_tree,
                    #merged_colors,
                    cbind(dynamic_colors, merged_colors),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Thalassiosira & Pseudo-nitz gene expression:\nmerged module thresh .62")

#rename modules back to numbers (although I dont think this is essential and can be cut for concision)
color_Order = c("grey", standardColors(70))
module_labels = match(merged_colors, color_Order)-1

merged_colors %>% data.frame() %>% filter(. == "grey") %>% rownames_to_column("ID") %>% separate(ID, into = c("tax", "KO"), sep = "-") %>% count(tax)
unique(module_labels)
### next objective is to correlate modules with environmental data / sampling day
n_genes = nrow(Norm_df)
n_samples = ncol(Norm_df)

#Recalculste MEs with color labels or use the merged MEs
merged_MEs = orderMEs(merged_MEs)

# Correlate modules with sampling days ------------------------------------
#Module-day heatmap
merged_MEs %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols = -sample, names_to = "module", values_to = "value") %>% 
  mutate(module = str_replace(module, regex("ME"), ""),
         sample = factor(sample, levels = c("Day1.A", "Day1.B", "Day1.C", 
                                            "Day3.A", "Day3.B", "Day3.C", 
                                            "Day5.A", "Day5.B", "Day5.C", 
                                            "Day9.A", "Day9.B", "Day9.C", 
                                            "Day11.A", "Day11.B", "Day11.C"))) %>% 
  filter(module != "grey") %>% 
  #filter(module %in% prelim_modules) %>% #this filter and vector was generated post hoc.
  ggplot(., aes(x = sample, y = module, fill = value)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient2(low = "black",
                       high = jisho_picker("eosine_pink"),
                       mid = "white",
                       midpoint = 0,
                       limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "",
       y = "Module", 
       x = "Sample",
       fill = "corr")
#These represent modules containing positive corrs in all replicates

#plot correlations  between eigengenes of each module -
# because i have kept all modules that have a consensus across replicates, this may enable me to group the modules according to similarity. IE, major growth module correlates with small growth modules that span one sampling period.
cormats = signif(cor(merged_MEs, use = "p"),2) %>% data.frame() %>% 
  rownames_to_column("ID") %>% column_to_rownames("ID") %>% 
  as.matrix() %>% 
  corrplot(.)


# Correlate modules with nutrient data ------------------------------------
# read in environmental table
EnvAuxData = read.csv("~/Dropbox/Project-data/SMPier/Figures/Publication/Env/EnvAuxData-2018_edited.csv") %>% 
  dplyr::select(Phosphate, Nitrate.nitrite,Day) %>% 
  mutate(rep = rep(c("A","B","C"), 5)) %>% 
  tidyr::unite(sample, c(Day, rep), sep = ".") %>% 
  column_to_rownames("sample")

#calculate correlation correlate modules with average nutrient concentration during the sampling day
module_traitCor = cor(merged_MEs %>% select(MEturquoise, MEsteelblue, MEorangered4, MEmagenta, MEgreen, MEdarkgreen, MEblue),
                      EnvAuxData,
                      use = "p") #pairwise.complete.obs

#calculate p values for the correlations.
module_trait_pvalue = corPvalueStudent(module_traitCor, n_samples) 
#or use corrplot
library(corrplot)
corrplot::corrplot(module_traitCor,
                   method = "circle",
                   p.mat = module_trait_pvalue,
                   addCoef.col = " black",
                   sig.level = 0.01,
                   insig = "blank",
                   pch.cex = 12,
                   pch.col = "red",
                   tl.col = "black",
                   tl.cex = .8,
                   tl.srt = 0,
                   cl.ratio = 0.4, #change the thickness of the color legend
                   col = COL2(diverging = 'RdBu', n = 10))
dev.off()

#Make a dataframe containing pvalues and significant correlations annotated.
module_trait_P_DF = module_trait_pvalue %>% data.frame() %>% 
  mutate(P_sig = ifelse(Phosphate <= 0.05, "sig", "insig"),
         N_sig = ifelse(Nitrate.nitrite <= 0.05, "sig", "insig")) %>% 
  rownames_to_column("module") %>% 
  pivot_longer(cols = contains("sig"), names_to = "signif", values_to = "value") %>% 
  filter(value == "sig")
# only blue grey and black have significance.

#Color code each association by the correlatio value:
sizeGrWindow(10, 6)
# merge the correlation and pvalues for plot
text_Matrix = paste(signif(module_traitCor, 2), "\n(",
                    signif(module_trait_pvalue, 1), ")", 
                    sep = "");

#plot heatmap of correlations with pvalues
par(mar = c(7, 11, 3, 5));
labeledHeatmap(Matrix = module_traitCor,
               xLabels = names(EnvAuxData[,-6]),
               yLabels = names(merged_MEs[c(-4,-6,-7)]),
               ySymbols = names(merged_MEs[c(-4,-6,-7)]),
               colorLabels = FALSE,
               colors = blueWhiteRed(60),
               textMatrix = text_Matrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("module-trait relatinships"))
dev.off()


# Validation from tutorial III-6 ------------------------------------------
# at this point i have modules that have been aggregated according to similarity based on eigengenes. 
# the modules themselves have individual genes with similar expression patterns.
signif(cor(merged_MEs, use = "p"), 2)

#define dissimilarity measures between eigengenes that tracks the sign of the 
# correlation between the module eigengenes, and se it to cluster the eigengene
mergedME_dist = (1-t(cor(merged_MEs, method = "p")))/2
hclust_mMEs = hclust(as.dist(mergedME_dist), method = "average")
#plot the eigengene dendrograms
par(mfrow = c(1,1))
plot(hclust_mMEs)
dev.off()
#Pairwise scatter plots: these amount to plo
sizeGrWindow(8,9)
plotMEpairs(merged_MEs[-4])

# A look at the relationship between brown and blue only
merged_MEs[-4] %>% 
  select(MEplum1, MEtan) %>% 
  ggplot(., aes(x = MEplum1, y = MEtan)) +
  geom_point(shape = 1) +
  geom_smooth(se = F, color = "red", lwd = .4, lty = 2, method = "lm") +
  theme_minimal()

#Diatnostic heatmaps
val_heatmap2 = function(module_vector) {
#Function to examine the fold changes in expression patterns within different modules
sizeGrWindow(8,9)
par(mfrow=c(4,1), mar=c(1, 2, 4, 1))
which.module=module_vector[1];
plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ),
        nrgcols=30,
        clabels=F,
        rlabels=T,
        rcols=which.module,
        title=which.module )
which.module=module_vector[2];
plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ),
        nrgcols=30,
        clabels=F,
        rlabels=T,
        rcols=which.module,
        title=which.module )
which.module=module_vector[3];
plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ),
        nrgcols=30,
        clabels=F,
        rlabels=T,
        rcols=which.module,
        title=which.module )
which.module=module_vector[4];
plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ),
        nrgcols=30,
        clabels=F,
        rlabels=T,
        rcols=which.module,
        title=which.module )
}

val_heatmap2(c("blue", "darkgreen", "green", "magenta"))
val_heatmap2(c("orangered4", "steelblue", "turquoise"))



# heatmap with expression profiles
val_heatmap = function(module){
  sizeGrWindow(8,7);
  which.module=module
  ME=merged_MEs[, paste("ME",which.module, sep="")]
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ) ,
          nrgcols = 15,
          clabels = c("Day1a","Day1b","Day1c",
                      "Day3a","Day3b","Day3c",
                      "Day5a","Day5b","Day5c",
                      "Day9a","Day9b","Day9c",
                      "Day11a","Day11b","Day11c"),
          rlabels = T,
          rcols = which.module,
          main = which.module, 
          cex.main = 2)
  par(mar = c(5, 4.2, 0, 0.7));
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="sample")
}

val_heatmap("lightcyan")
val_heatmap("greenyellow")
val_heatmap("darkorange")
val_heatmap("darkgreen") 

val_heatmap_taxfilt = function(module, tax){
  #visualize the expression of genes corresponding to one taxa
  #tax uses str_detect(), so fragments of confident pieces the name is better than whole mispelled name
  #val_heatmap_taxfilt("yellow", "Thalas")-Thalassiosira gene expression for genes in the yellow module
  sizeGrWindow(8,7);
  which.module=module
  ME=merged_MEs[, paste("ME",which.module, sep="")]
  par(mfrow=c(1,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(Norm_df.t[,merged_colors==which.module ]) ) %>% data.frame() %>% filter(str_detect(rownames(.), tax)) %>% as.data.frame.matrix() ,
          nrgcols = 15,
          clabels = c("Day1a","Day1b","Day1c",
                      "Day3a","Day3b","Day3c",
                      "Day5a","Day5b","Day5c",
                      "Day9a","Day9b","Day9c",
                      "Day11a","Day11b","Day11c"),
          rlabels = T,
          rcols = which.module,
          main = which.module, 
          cex.main = 2)
}
val_heatmap("lightcyan")
val_heatmap_taxfilt("orangered4", "Thalas")
val_heatmap_taxfilt("orangered4", "Pseudo")

val_heatmap("greenyellow") 
val_heatmap_taxfilt("greenyellow", "Thalas")
val_heatmap_taxfilt("greenyellow", "Pseudo")

val_heatmap("darkgreen")
val_heatmap_taxfilt("darkgreen", "Thalas")
val_heatmap_taxfilt("darkgreen", "Pseudo")

val_heatmap("darkorange")
val_heatmap_taxfilt("darkorange", "Thalas")
val_heatmap_taxfilt("darkorange", "Pseudo")

#Number of genes inside of each module

# Note: Further examining KOs within a module -----------------------------
#Use the module_DF df to further examine the Function & logFC of the KOs within a module.
module_DF = data.frame("module" = merged_colors) %>% 
  rownames_to_column("ID")

# Pull out all transcripts associated with the darkgreen module
module_df_getter = function(MODULE){
#This function produces dataframe containing all transcripts with annotations for a given module and saves to a .csv.
  #Groom the dataframe
df_out = module_DF %>% 
  filter(module == MODULE) %>% #filter out transcripts in the darkgreen module
  left_join(.,Norm_df %>% rownames_to_column("ID"), by = "ID") %>% #join the expression profiles.
  #What follows is averaging across replicates for each day
  pivot_longer(cols = contains("Day"), names_to = "sample", values_to = "count") %>% 
  separate(sample, into = c("day", "rep"), sep = "\\.") %>% 
  dplyr::select(-rep) %>% 
  group_by(ID, module, day) %>% 
  summarise(count = round(mean(count, na.rm = T))) %>% 
  ungroup() %>% 
  dcast(ID+module~day) %>% 
  #back to it...
  separate(ID, into = c("tax", "KO"), sep = "-") #%>% #separate the ID to use the KO for joining
  #left_join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T)) %>%
  #left_join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t"))

FILENAME = str_c(MODULE,"module.csv", sep = "")
path = str_c("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/",FILENAME, sep = "") 

write.csv(df_out, file = path, row.names = F)

return(str_c(FILENAME,"saved, Dr. Oc!", sep = " "))
}

#Automate the module_getter function by iterating over the modules of interest.
#Create a vector containing modules of interest to iterate over.
module_vector = c("steelblue", "green", "darkgreen", "blue", "turquoise", "orangered4", "magenta")

#Use the map function to iterate over the vector as input to the module_df_getter() function.
map(module_vector, module_df_getter)

## Side mission- See SM_1-WGCNA_assembly_ThalasPnitz2018.R

module_composition = function(MODULE) {
  #this function takes as input a module color and prints the taxonomic compositon of the module
  module_DF %>% 
  filter(module == MODULE) %>% 
  separate(ID, c("tax","KO"), sep = "-") %>% 
  count(tax) #%>%  
  #arrange(desc(n)) #%>% 
  #ggplot(., aes(x = tax, y=n)) +
  #geom_bar(stat = "identity")
}
table(merged_colors)

#What is the compositon of all of the modules?
#use the map function to run the module_composition() function on all modules.
x = map(module_vector, module_composition)
names(x) = module_vector

# Make a table of compositions
module_composition_DF = data.frame(tax = x$steelblue[1],
           "steelblue" = x$steelblue[,2],
           "green" = x$green[,2],
           "darkgreen" = x$darkgreen[,2],
           "blue" = x$blue[,2],
           "turquoise" = x$turquoise[,2],
           "orangered4" = x$orangered4[,2],
           "magenta" = x$magenta[,2])

#Plot the module compositions
module_composition_DF %>% 
  pivot_longer(cols = -tax, names_to = "module", values_to = "value") %>% 
  ggplot(aes(x = module, y = value, fill = tax)) +
  geom_bar(stat = "identity", 
             position = "dodge",
             #position = position_dodge(0.6), width = .5, 
             size = 9) +
  #scale_y_continuous(breaks = seq(0,1500,500), minor_breaks = c(seq(0,500,50),seq(500,1500,100))) +
  scale_y_continuous(breaks = seq(0,1500,100)) +
  coord_cartesian(ylim = c(50,1500)) +
  scale_fill_manual(name = NULL,
                    breaks = c("Pseudo_nitz", "Thalassiosira"),
                    values = pal1[c(1,3)]) +
  labs(y = "No. Kegg IDs",
       x = "Module") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12))

ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/module_sizes2018.pdf", width = 6, height = 4)

## Pulling the KO ids in each of the modules - From there i can do a few things:
#1. examine the cumulative functional enrichment/vibe and
#2. juxtapose the functional enrichment/vibe of the two taxa separately


# Functional Enrichment/ Over Representation Analysis ---------------------
#Resource1: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
#Resource2: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html

#load in the clusterprofiler package for functional enrichment analysis
library(clusterProfiler)

#iterate over all data
# read in module dataframes
Enrichment_DFs = rbind(read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/bluemodule.csv"), 
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/darkgreenmodule.csv"),
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/greenmodule.csv"),
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/magentamodule.csv"),
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/orangered4module.csv"),
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/steelbluemodule.csv"),
                       read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/turquoisemodule.csv")) %>% 
  #Separate out by module
  unite(ID, tax, module, sep = "-", remove= F) %>% 
  split(.$ID) %>% 
  #Pull the vectors of KOs for each module's dataframe
  map(., ~pull(., KO)) %>% 
  #run the over representation analysis on each vector
  map(., ~enrichKEGG(gene = .,
                     organism = 'ko', #must be used with kegg orthology
                     pvalueCutoff = 0.01)) %>% 
  #convert to a (list of) dataframe(s)
  map(., ~data.frame(.)) %>% 
  #and arrange the dataframes by the Count column
  map(., ~arrange(., desc(Count)))

#Create a dot plot of modular enrichments
Enrich_table = bind_rows("Pnitz_blue" = Enrichment_DFs$`Pseudo_nitz-blue`, 
                         "Pnitz_darkgreen" = Enrichment_DFs$`Pseudo_nitz-darkgreen`, 
                         "Pnitz_green" = Enrichment_DFs$`Pseudo_nitz-green`, 
                         "Pnitz_magenta" = Enrichment_DFs$`Pseudo_nitz-magenta`,
                         "Pnitz_orangered4" = Enrichment_DFs$`Pseudo_nitz-orangered4`, 
                         "Pnitz_steelblue" = Enrichment_DFs$`Pseudo_nitz-steelblue`, 
                         "Pnitz_turquoise" = Enrichment_DFs$`Pseudo_nitz-turquoise`,
                         "Thalas_blue" = Enrichment_DFs$`Thalassiosira-blue`, 
                         "Thalas_darkgreen" = Enrichment_DFs$`Thalassiosira-darkgreen`, 
                         "Thalas_green" = Enrichment_DFs$`Thalassiosira-green`, 
                         "Thalas_magenta" = Enrichment_DFs$`Thalassiosira-magenta`,
                         "Thalas_orangered4" = Enrichment_DFs$`Thalassiosira-orangered4`, 
                         "Thalas_steelblue" = Enrichment_DFs$`Thalassiosira-steelblue`, 
                         "Thalas_turquoise" = Enrichment_DFs$`Thalassiosira-turquoise`,
                         .id = "module") %>% 
  dplyr::select(module, Description, GeneRatio, Count, p.adjust) %>%
  separate(module, into = c("tax", "module"), sep = "_") %>% 
  separate(GeneRatio, into = c("num", "denom"), sep = "/") %>% 
  mutate(GeneRatio = as.integer(num)/as.integer(denom)) %>% 
  #filter(GeneRatio > 0.025) %>% 
  dplyr::select(-num, - denom)

Human_diseases = c("Spinocerebellar ataxia", "Alzheimer disease", "Amyotrophic lateral sclerosis", "Chemical carcinogenesis - reactive oxygen species", "Epstein-Barr virus infection", "Fanconi anemia pathway", "Huntington disease", "Insulin signaling pathway", "Non-alcoholic fatty liver disease", "Parkinson disease", "Prion disease", "Pathways of neurodegeneration - multiple diseases", "Diabetic cardiomyopathy","Vibrio cholerae infection")
Amino_acid_metab = c("Alanine, aspartate and glutamate metabolism", "Cysteine and methionine metabolism")
Amino_acid_breakdown = c("Valine, leucine and isoleucine degradation")
#Nucleotide_breakdown = c("RNA degradation")
Nucleotide_metabolism = c("Purine metabolism", "Pyrimidine metabolism")
Autophagy = c("Autophagy - yeast", "Autophagy - animal")

tax_lookup = c(Pnitz = "Pseudo-nitzschia",
               Thalas = "Thalassiosira") #Havent tried this yet. Should relabel the tabs in the facet wrap
# Make a plot of the overrepresentation analysis results
Enrich_table %>% 
  #distinct(Description) %>% arrange(desc(.))
  #mutate statements are manually aggregating and renaming for ease of interpretation.
  mutate(FxTerm = case_when(str_detect(Description, regex("autoph", ignore_case = T)) ~"Autophagy",
                            str_detect(Description, regex("COVID", ignore_case = T)) ~"Ribosome",
                            str_detect(rownames(.), "(ko00062)|(ko01212)|(ko0414[02])|(ko04966)|(ko05016)|(ko05415)") ~"Lysosome binding and processing",
                            str_detect(rownames(.), "(ko00020)|(ko00220)|(ko00250)") ~"Organic Nuptake and assimilation",
                            str_detect(rownames(.), "(ko012[034]0)") ~"Organic N and P metabolism",
                            str_detect(rownames(.), "(ko04721)|(ko05110)") ~"Oxidative phosphorylation - Lysosome V-type ATPase",
                            Description %in% Human_diseases ~"Human Disease Pathways",
                            #Description %in% Nucleotide_metabolism ~"Nucleotide metabolism",
                            #Description %in% Nucleotide_breakdown ~"Nucleotide breakdown",
                            str_detect(Description, "Cell cycle") ~"Cell cycle",
                            str_detect(Description, regex("(proteolysis)|(proteosome)", ignore_case = T)) ~"Proteolysis")) %>% 
  mutate(FxTerm = ifelse(is.na(FxTerm), Description, FxTerm)) %>% 
  mutate(FxTerm = str_replace(FxTerm, "Biosynthesis of amino acids", "Amino acid biosynthesis"),
         FxTerm = str_replace(FxTerm, "Ribosome$", "Ribosomal protein - Translation"),
         FxTerm = str_replace(FxTerm, "Ribosome biogenesis .+", "Ribosome biogenesis"),
         FxTerm = str_replace(FxTerm, "Meiosis .+", "Meiosis"),
         FxTerm = str_replace(FxTerm, "Proteasome", "Proteolysis"),
         FxTerm = str_replace(FxTerm, "Alanine, aspartate and glutamate metabolism", "GS/GOGAT"),
         module = factor(module, levels = c("steelblue", "darkgreen", "green", "blue", "orangered4", "magenta", "turquoise"))) %>% 
  group_by(tax, module, FxTerm) %>% 
  summarise(GeneRatio = sum(GeneRatio, na.rm = T)) %>% 
  filter(GeneRatio > 0.025) %>% 
  #filter(FxTerm %in% hilight) %>% 
  ggplot(., aes(y = FxTerm, x = module, size = GeneRatio, fill = tax)) +
  geom_point( shape = 21) +
  scale_size_continuous(breaks = c(0.03, 0.14, 0.28, 0.56),
                        labels = c(0.03, 0.14, 0.28, 0.56)) +
  scale_fill_manual(name = "",
                    values = pal1) +
  facet_wrap(~tax,
             ncol = 2,
             labeller = labeller(tax = tax_lookup)) +
  labs(x = "",
       y = "") +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = "solid"),
        panel.grid.major.y = element_line(linetype = "dashed"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))


# Aggregate Enrichment tables and save to a .csv file for the purposes of manual assignment of functional terms. Microbes dont get COVID-19
x = bind_rows("Pnitz_blue" = Enrichment_DFs$`Pseudo_nitz-blue`, 
              "Pnitz_darkgreen" = Enrichment_DFs$`Pseudo_nitz-darkgreen`, 
              "Pnitz_green" = Enrichment_DFs$`Pseudo_nitz-green`, 
              "Pnitz_magenta" = Enrichment_DFs$`Pseudo_nitz-magenta`,
              "Pnitz_orange" = Enrichment_DFs$`Pseudo_nitz-orange`,
              "Pnitz_paleturquoise" = Enrichment_DFs$`Pseudo_nitz-paleturquoise`,
              "Thalas_blue" = Enrichment_DFs$`Thalassiosira-blue`, 
              "Thalas_darkgreen" = Enrichment_DFs$`Thalassiosira-darkgreen`, 
              "Thalas_green" = Enrichment_DFs$`Thalassiosira-green`, 
              "Thalas_magenta" = Enrichment_DFs$`Thalassiosira-magenta`,
              "Thalas_orange" = Enrichment_DFs$`Thalassiosira-orange`,
              "Thalas_paleturquoise" = Enrichment_DFs$`Thalassiosira-paleturquoise`,
              .id = "module") %>% 
  #dplyr::select(module, Description, GeneRatio, Count, p.adjust) %>%
  separate(module, into = c("tax", "module"), sep = "_") %>% 
  separate(GeneRatio, into = c("num", "denom"), sep = "/") %>% 
  mutate(GeneRatio = as.integer(num)/as.integer(denom)) 
x %>% distinct(Description) %>% arrange(Description)
mutate(Description = case_when(Description == "Coronavirus disease - COVID-19" ~ Ribosome))
write.csv(x, file = "~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/Enrichment_DF.csv")

### Gene set enrichment analysis (GSEA)
#this is just a prototype to illustrate to myself how to use the tool
#read in one table (paleturquoise module)
table = read.csv("~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/bluemodule.csv") 
#the genes in the paleturquoise module exhibited higher expression on day 9 so we are using this day for the gene expression levels.
colSums(table[4:8]) #day9 is the highest

#gseKEGG requires a sorted vector with expression and KOId as the names
gList = table %>% 
  dplyr::select(KO, Day5) %>% #grab KO and Day9 expression
  arrange(desc(Day5)) #sort

# Create a vector
vec = gList %>% 
  pull(Day5)
# Apply KOids for the names of the variables
names(vec) = gList$KO

#Run gene set enrichment analysis
x = gseKEGG(vec, organism = 'ko')


# Example: conduct enrichment analysis on the blue module
#Function to pull vector of the KOs -by default it will grab both Thalassiosira and Pseudo_nitz
VectorOfKOs = function(MODULE, TAXvec = c("Thalassiosira", "Pseudo_nitz")) {
  module_DF %>% 
    separate(ID, c("tax","KO"), sep = "-") %>%
    filter(module == MODULE,
           tax %in% TAXvec) %>% 
    pull(KO)
}

bluevec = VectorOfKOs("blue", "Thalassiosira") 
benrich = enrichKEGG(bluevec,
           organism = 'ko',
           pvalueCutoff = 0.05)

barplot(benrich,
        drop = T,
        showCategory = 10,
        font.size = 10)

#End WGCNA



# Differential Expression analysis (edgeR) -------------------------------------------------------------
# the goal here is to produce a vector of KOID that are differentially expressed over the course of the samples.
# Filter out low count rows using a conditional filter--that is, to keep rows that meet a certian condition.
keep = filterByExpr(dge_obj)
dge_obj = dge_obj[keep, ,keep.lib.sizes=FALSE]
# setup the design matrix without an intercept as day 1!
design = model.matrix(~0+group, data = dge_obj$samples)

#Set the column names for the design to match the
colnames(design) = levels(dge_obj$samples$group)

# using makeContrasts() to specify the pariwise comparisons
conts = makeContrasts(day3-day1, day5-day1, day9-day1, day11-day1, levels = str_to_lower(sample_list)) #for contrasting all samples to prebloom state

# estimate common and tagwise(pairwise) dispersion accrding the study design matrix
disp = estimateGLMCommonDisp(dge_obj, design = design)
disp = estimateGLMTagwiseDisp(disp, design = design)

# Determine differentially expressed genes
fit = glmQLFit(disp, design = design)
DEGs = glmQLFTest(fit, contrast = conts)
DEG_table = DEGs$table %>% data.frame()

# adjust p values due to so many comparissons to decrease chance of type 1 error.
DEG_table$P.adjusted = p.adjust(DEG_table$PValue, method = "fdr")
DEG_table$KO = rownames(DEG_table)
DEG_table = DEG_table %>% 
  filter(P.adjusted < 0.001) %>% 
  separate(KO, into = c("tax","KO"), sep = "-")

# The point of this analysis is to investigate diatom gene expression during this bloom so lets filter out the diatom reads
# quality filter the DEG table to those with p < 0.01 and attach functional annotations
nrow(DEG_table)
DDegs_wcat = DEG_table %>% 
  left_join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T)) %>%
  left_join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t"))  

#how many diatom DEGs?
nrow(DDegs_wcat)# There are 1626 Diatom reads with evidence of differential expression


#save(DDegs_wcat, file = "~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/DEG_table-assembly2018.rda")
#load("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/DEG_table-assembly2018.rda")

# Arrange data for plotting.
##use the below for the full contrast list
# Or everything contrasted against day 1
##use the below for the full contrast list
Diatom_Degs_wcat = DDegs_wcat %>%
  dplyr::rename("Day3" = logFC.day3...day1,
                "Day5" = logFC.day5...day1,
                "Day9" = logFC.day9...day1,
                "Day11" = logFC.day11...day1) %>%
  pivot_longer(cols = c("Day3", "Day5", "Day9", "Day11"), names_to = "contrasts", values_to = "fold_chng") %>% 
  dplyr::select(-logCPM, -F, -PValue, -A, -D, -Target, -Notes, -X) %>%
  dplyr::select(KO, Target_2, Category, B, contrasts, fold_chng, everything())

#Count the number of differentially expressed KOs
Diatom_Degs_wcat %>% filter(contrasts == "Day5") %>%
  count(KO) %>% 
  arrange(desc(n)) %>% 
  pull(n) %>% 
  sum()
# there are 1674 diatom DEGs.

#Pull the list of KOids of differentially expressed genes.
DEG_KOs = Diatom_Degs_wcat %>% distinct(KO) %>% pull()


# circle back: how many of the x KOs in the black module are differentially expressed?
##note: this will further reduce the list of targets. Is this desirable given the lack of information for each gene?
##note: remember, the goal ultimate is to determine the physiological mode during each stage of the bloom.



# DEG sandbox -------------------------------------------------------------
## The intensity of gene expression: There are different tiers of gene expression, and there must be a cutoff for what I consider major.
## how many differentially expressed genes are there with an logFC greatr than |1|? greater than |7|?
# I didn't select a higher logFC for two reasons:
# -We dont know how much of a logFC is enough to produce a significant cellular reaction.
# -Better chance of finding genes once cross referenced with the WGCNA modules.

#How many DEGs with logFCs greater than |1| on each day?
Diatom_Degs_wcat$direction = "minor" #Create 
Diatom_Degs_wcat$direction[Diatom_Degs_wcat$fold_chng > 1] = "upreg"
Diatom_Degs_wcat$direction[Diatom_Degs_wcat$fold_chng < -1] = "downreg"

#major_FCs is a table displaying the number and fraction of up and down > and < |1| 
major_FCs = Diatom_Degs_wcat %>%
  group_by(contrasts) %>%
  count(direction) %>%
  pivot_wider(names_from = direction, values_from = n) %>%
  mutate(totalDE = sum(downreg,upreg, na.rm = TRUE),
         fracDE = totalDE / sum(downreg, minor, upreg, na.rm = T),
         fracUp = upreg / sum(downreg, minor, upreg, na.rm = T),
         fracDn = downreg / sum(downreg, minor, upreg, na.rm = T),
         fracM = minor / sum(downreg, minor, upreg, na.rm = T)) %>% ungroup()
write.csv(major_FCs, file = "~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/logFC_table.csv", row.names = F)
#use permutations of the below to ask specific questions about the proportion of up down and muted DEGs
major_FCs %>% filter(contrasts != "Day9") %>% select(fracDn) %>% pull() %>% mean(., na.rm = T)

#function to plot custom pals
jisho = function(color){
  
  color = jishonoiro[[color]] %>% data.frame() %>% pull(.)
  return(color)
}

# plot this information
Diatom_Degs_wcat %>%
  #filter(KO %in% D5_subexpr$KO) %>%
  mutate(hilight = abs(fold_chng) > abs(1),
         contrasts = factor(contrasts, levels = c("Day3","Day5","Day9","Day11"))) %>%
  mutate(hilight = ifelse(abs(fold_chng) > abs(7), "super", hilight)) %>% 
  ggplot(., aes(x = contrasts, y = fold_chng)) +
  geom_boxplot(alpha = 0.1, outlier.size = 0, lwd = .2, aes(color = hilight), show.legend = F) +
  #geom_boxplot(alpha = .4, outlier.size = .9, lwd = .4,  show.legend = T, outlier.colour = "white", outlier.shape = "triangle") +
  #stat_summary(fun = "mean", shape = 23, size = .4, fill = "black") +
  geom_jitter(width = .2, shape = 21, color = "black", alpha = 0.6, show.legend = F, aes(fill = hilight)) +
  labs(y = "Log Fold Change",
       x = "Contrasts") +
  scale_fill_manual(values = jisho("american")) +
  scale_color_manual(values = jisho("american")) +
  theme_minimal()
#or as a boxplot
jishonoiro[["american"]]
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/Thalas_pnitz_logFC_jitter.pdf", width = 8, height = 8)
save.image("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/WGCNA_DEG_assemply2018(clean).RData")




# Start DEGs + WGCNA modules ---------------------------------------------------------
# Examine DE of module genes (next examine which module genes are DE)
###Use this function to plot the logFC of differentially expressed genes as a jitterplot
mod_scatterplot = function(MODULE){
  ##Examining expression profiles
  submods = data.frame(gene_id = names(merged_colors), colors = merged_colors) %>%
    filter(colors %in% MODULE)
  
  subexpr = Norm_df[submods$gene_id,] %>% data.frame()
  
  subexpr$KO = row.names(subexpr)
  subexpr_wcat = subexpr %>%
    left_join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T)) %>%
    left_join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t")) %>%
    select(-starts_with("SMPier"))
  
  
  plt = Diatom_Degs_wcat %>%
    filter(KO %in% NN_KO_vector) %>% #filter for the KOIDs found to be positively correlated with day 5 samples
    mutate(hilight = abs(fold_chng) > abs(1),
           contrasts = factor(contrasts, levels = c("Day3","Day5","Day9","Day11"))) %>%
    mutate(hilight = ifelse(abs(fold_chng) > abs(7), "super", hilight)) %>%
    ggplot(., aes(x = contrasts, y = fold_chng)) +
    #geom_boxplot(alpha = 0.1, outlier.size = 0, lwd = .2, aes(color = hilight), show.legend = F) +
    geom_jitter(width = .2, shape = 21, color = "black", alpha = 0.6, show.legend = F, aes(fill = hilight)) +
    labs(y = "Log Fold Change",
         x = "Contrasts") +
    scale_fill_manual(values = american) +
    scale_color_manual(values = american) +
    theme_minimal()
  
  return(plt)
}

###Use this function to plot the expression of the module genes s a line plot
expression_lineplot = function(MODULE){
  module_df = data.frame("module" = merged_colors) %>% 
    rownames_to_column("ID")
  
  #names(module) = MODULE
  submods = module_df %>%
    filter(module %in% MODULE)
  
  #change rownames
  row.names(module_df) = module_df$ID
  
  #grab the rows of the original normalized count table that correspond to the day5 submodules
  subexpr = Norm_df[submods$ID,] %>% data.frame()
  #length(D5_submods$gene_id) #sanity check
  #dim(subexpr) #sanity check

  subexpr %>%
    mutate(gene_id = row.names(.)) %>%
    pivot_longer(-gene_id) %>%
    mutate(module = module_df[gene_id,]$module,
           name = str_replace(name, "SMPier\\.", ""), #for values in the name column, replace "SMPier." with nothing
           name = str_replace(name, regex("\\.\\d{2}\\.\\d{2}\\.2018_S\\d+"), ""), #replace ".##.##.2018_S########..." with nothing
           #factor to specify the order of the levels
           name = factor(name, levels = c("Day1.A", "Day1.B", "Day1.C",
                                          "Day3.A", "Day3.B", "Day3.C",
                                          "Day5.A", "Day5.B", "Day5.C",
                                          "Day9.A", "Day9.B", "Day9.C",
                                          "Day11.A", "Day11.B", "Day11.C"))) %>%
    separate(name, into = c("day", "rep")) %>% 
    group_by(gene_id, day, module) %>% 
    summarise(value = mean(value, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(day = factor(day, levels = c("Day1", "Day3", "Day5", "Day9", "Day11"))) %>% 
    filter(str_detect(gene_id, "Thalas")) %>% 
    ggplot(., aes(x=day, y=value, group=gene_id)) +
    geom_line(size = .5,alpha = 0.8, aes(color = module)) + 
    #geom_point(size = .1, alpha = .3, color = "black") +
    scale_color_manual(values = MODULE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(rows = vars(module), scales = "free") +
    labs(x = "sample",
         y = "normalized expression")
}
expression_lineplot(c("blue", "darkgreen", "green", "magenta", "orange", "paleturquoise"))
ggsave(plot = last_plot(), filename = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/Redo_Jan2023/Module_expression_profiles.pdf", width = 8, height = 6)
###which modules in module_df have evidence of differential expression? Might be worth pruning further upstream.
no.moduleDEGs = function(module) {
  
  no.degs = module_df %>%
    filter(gene_id %in% pull(Diatom_Degs_wcat %>% distinct(KO)),
           colors %in% module) %>% nrow()
  mods = module_df %>%
    filter(gene_id %in% pull(Diatom_Degs_wcat %>% distinct(KO)),
           colors %in% module) %>%
    distinct(colors) %>% pull(colors)
  dat = list("number" = no.degs,
             "modules" = mods)
  return(dat)
}

#"day1 contains the colors that were highly correlated with day 1
mod_scatterplot("orange")
expression_lineplot("orange")
no.moduleDEGs(day1) 


# Circle back for DEGs in module ------------------------------------------
# the original question was which genes witin modules of interest were differentially expressed
df %>% 
  mutate(hilight = ifelse(GS > 0.25 & MM > 0.35, "sig", "not sig")) %>% 
  ggplot(., aes(x = MM, y = GS, color = hilight)) +
  geom_point(shape = 21) +
  scale_color_manual(name=NULL,
                     values = c("grey", "steel blue")) +
  labs(x = "Module membership",
       y = "Gene significance for Nitrate+nitrite") +
  theme_minimal()

# Which of these KOs are differentially expressed? And what's their function
NN_KO_vector = nitrateKOs %>% pull(KO)
Day3moduleDEGS = nitrateKOs %>% 
  filter(KO %in% DEG_KOs) %>% 
  left_join(., DEG_table, by = "KO") %>% 
  filter(P.adjusted < 0.001) %>% 
  select(-logCPM, -F, -PValue, -Notes, -X, -Phagotrophy_annotation, -P.adjusted) %>% 
  pivot_longer(cols = starts_with("logFC"), names_to = "contrast", values_to = "logFC") %>% 
  filter(str_detect(contrast, "logFC.day3")) %>% 
  select(KO, logFC, contrast, GS.Nitrate.nitrite, everything()) %>% 
  mutate(Category = ifelse(Category %in% c("TCA cycle", "Glyoxylate cycle"),"TCA-Glyoxylate cycle",Category)) 

x = Day3moduleDEGS %>% 
  select(Category) %>% rename(Category2 = Category)

x[is.na(x)] = "0"
Day3moduleDEGS = cbind(Day3moduleDEGS,x)

Day3moduleDEGS %>% 
  mutate(Category2 = ifelse(Category2 == "0", B, Category2)) %>% 
  mutate(Category2 = str_replace_all(Category2, regex("nucleotide and amino.+", ignore_case = T), "Nuc-Amino metabolism")) %>% 
  group_by(KO,Category2) %>% 
  summarise(lgFC = median(logFC, na.rm = TRUE)) %>% ungroup() %>% arrange(desc(lgFC)) %>% 
  ggplot(., aes(y = Category2, x = lgFC)) +
  geom_point()

#Function to grab a list of KOs in a module - note NO DIFFERENTIAL EXPRESSION INFO.
grab_KOs = function(mod_color){
  out = module_DF %>% 
    filter(module == mod_color) %>% 
    plyr::join(., read.csv("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/K0_geneIDs_fulllist_08192019.csv", header = T), type = "left", match = "first") %>% 
    plyr::join(., read.delim("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/Custom_KO_list_30-04-2020.txt", header = TRUE, sep = "\t"), type = "left", match = "first") %>% 
    mutate(Category = ifelse(Category %in% c("TCA cycle", "Glyoxylate cycle"),"TCA-Glyoxylate cycle",Category))

  # I want to replace the NA values in Category with 0 in order to replace them as an explicit character in an ifelse statement.
  x = out %>% 
    select(Category) %>% rename(Category2 = Category)
  x[is.na(x)] = "XXX"
  out = cbind(out,x)
  
  out = out %>% 
    mutate(Category2 = ifelse(Category2 == "XXX", B, Category2)) %>% 
    mutate(Category2 = str_replace_all(Category2, regex("nucleotide and amino.+", ignore_case = T), "Nuc-Amino metabolism"))
  return(out)
}

#Growth module - should contain pathways related to responses to replete nutrients
grab_KOs("black") %>% write_csv(file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/black_module(1-3growth).csv")
grab_KOs("lightcyan") %>% write_csv(file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/lightcyan_module(3growth).csv")
grab_KOs("black") %>% filter(Category2 == "Genetic information processing") %>% count(C) %>% arrange(desc(n))
count(Category2) %>% arrange(desc(n)) %>% 
  filter(Category2 != "<NA>") %>% 
  mutate("frac" = n/sum(n)) %>% view()

#Transition module
grab_KOs("cyan") %>% write_csv(file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/cyan_module(5transition).csv")
grab_KOs("cyan") %>% filter(Category2 == "Genetic information processing") %>% count(C) %>% arrange(desc(n))
count(Category2) %>% arrange(desc(n)) %>% 
  filter(Category2 != "<NA>") %>% 
  mutate("frac" = n/sum(n)) %>% pull(n) %>% sum()

#Death module - should contain pathways related to adaptations to no nutrients
grab_KOs("lightgreen") %>% write_csv(file = "~/Dropbox/Project-data/SMPier/Transcriptomics/Figures/assembly_2018/WGCNA/ThalasPnitz/lightgreen_module(9-11death).csv")
grab_KOs("lightgreen") %>% filter(Category2 == "Genetic information processing")%>% count(C) %>% arrange(desc(n))
  count(Category2) %>% arrange(desc(n)) %>% 
  filter(Category2 != "<NA>") %>% 
  mutate("frac" = n/sum(n)) %>% 
#END
save.image("WGCNA_DEGs-SMP2018.RData")
load("WGCNA_DEGs-SMP2018.RData")
#Sandbox

#create a dataframe containing the module and KO ids
#This dataframe can be used to bundle up KO ids associated with specific modules for downstream analysis
# using filter() %>% pull() statements


jitter_test = function(col){
# module_DF must be inherited from global env.
  KO_list = module_DF %>% filter(module == col) %>% pull(KO)
  
  plt = Diatom_Degs_wcat %>% 
    filter(KO %in% KO_list) %>% #filter for the KOIDs found to be positively correlated with day 5 samples
    mutate(hilight = abs(fold_chng) > abs(1),
           contrasts = factor(contrasts, levels = c("Day3","Day5","Day9","Day11"))) %>%
    mutate(hilight = ifelse(abs(fold_chng) > abs(7), "super", hilight)) %>%
    ggplot(., aes(x = contrasts, y = fold_chng)) +
    #geom_boxplot(alpha = 0.1, outlier.size = 0, lwd = .2, aes(color = hilight), show.legend = F) +
    geom_jitter(width = .2, shape = 21, color = "black", alpha = 0.6, show.legend = F, aes(fill = hilight)) +
    labs(y = "Log Fold Change",
         x = "Contrasts") +
    scale_fill_manual(values = c("lightgrey", "#A90636", "steelblue")) +
    scale_color_manual(values = c("lightgrey", "#A90636", "steelblue")) +
    theme_minimal()
  return(plt)
}
# Testing differential expression of genes in the black module - mostly negative logFCs in 5,9,11
jitter_test("black") 

# Testing differential expression patterns of genes in the cyan module - pronounced positive logFC in 5
jitter_test("cyan")

# Testing the differential expressio patterns of genes in the lightgreen module - pronounced positive logFC in 9 and 11
jitter_test("lightgreen")

module_df %>% filter(module == "magenta") %>% 
  filter(KO %in% NN_KO_vector) %>% pull(KO)
KO_list

length(NN_KO_vector)




