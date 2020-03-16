
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("phyloseq")
library(phyloseq)

packageVersion("phyloseq")

library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library(tidyr)


setwd("/Users/mahmoudi/Documents/Doc/MothurR/")  #set working directory
Bacteria <- read.table("Dataset/BacV5.trim.contigs.good.renamed.unique.pick.cutadapt2.unique.good.pick.dgc.0.03.abund.0.03.pick.shared", header = T , check.names = F , stringsAsFactors = F)
BacteriaTaxa <- read.table("Dataset/BacV5.trim.contigs.good.renamed.unique.pick.cutadapt2.unique.good.pick.dgc.0.03.abund.0.03.cons.pick.cons.taxonomy",header = T , check.names = F , stringsAsFactors = F)
sampleinfo <- read_excel("Sample mapping/Lib_prep_031017_edited.xlsx" , sheet = 3)
sampleinfo <- as.data.frame(sampleinfo)
Bac1 <- Bacteria%>%separate("Group" , c("Loci" , "LibM") , remove = FALSE)
sample_df <- merge(Bac1[4] , sampleinfo , by.x = "LibM" , by.y = "Lib")
Bac3 <- merge(Bac1 , sampleinfo , by.x = "LibM" , by.y = "Lib")
otu_mat <- Bac3[-c(2:5,(ncol(Bac3)-9):ncol(Bac3))]
rownames(otu_mat) <- otu_mat$LibM
otu_mat <- t(otu_mat[-1])
BacteriaTaxa$Taxonomy = gsub("\\s*\\([^\\)]+\\)","",as.character(BacteriaTaxa$Taxonomy))
otu_taxa <- BacteriaTaxa %>% separate("Taxonomy" , c("Kingdom" , "Phylum" , "Class" , "Order" , "Family" , "Genus" ,"Species" ) , ";" , remove = FALSE)
rownames(otu_taxa) <- otu_taxa$OTU
otu_taxa <- otu_taxa[-c(1,2)]
#sample_df <- sampleinfo[apply(sampleinfo,1,function(x)any(!is.na(x))),]
rownames(sample_df) <- sample_df$LibM
sample_df <- sample_df[-1]

otu <- otu_table(otu_mat , taxa_are_rows = TRUE,errorIfNULL = TRUE)
tax <- tax_table(as.matrix(otu_taxa))
samples <- sample_data(sample_df)
carbom <- phyloseq(otu,tax,samples)
carbom2 <- subset_taxa(carbom1, Phylum=="p__Firmicutes")
carbom1 <- subset_samples(carbom, Infection =="y")
carbom <- subset_taxa(carbom, !(Class %in% c("Syndiniales", "Sarcomonadea")))
carbom
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom1 = transform_sample_counts(carbom, standf)
y = plot_bar(carbom1, fill = "Kingdom")

sample_variables(carbom)
plot_bar(carbom, "Compartment", "Abundance", "Phylum", title=title)
get_taxa(physeq = carbom,"Compartment")


data("GlobalPatterns")
gpt <- subset_taxa(GlobalPatterns, Kingdom=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:300]), gpt)
plot_heatmap(gpt, sample.label="SampleType")

plot_heatmap(carbom2)
x = carbom2()
plot_richness(carbom2, x="Site", color="Site")
plot_bar(carbom2, fill="Phylum")
 plot_bar(carbom, x="Site", fill="Phylum")
ggplot()

y = as.data.frame(colSums(otu_mat))
colnames(y) = "abundance"
ggplot() + geom_bar(data = Bac3 , x="Site" , y = "y")  
 xx = plot_bar(carbom2, fill="Genus" , x = "Compartment") + geom_boxplot()
ggsave("file.pdf" , xx  , width = 20, height = 8 , dpi = 1000 )  
prevotella = subset_taxa(carbom, Genus =="g__Pseudomonas")  
plot_tree(prevotella, plot.margin = 0.5, ladderize = TRUE)  
top20OTU.names = names(sort(taxa_sums(carbom), TRUE)[1:20]) 
plot_heatmap(top20OTU.names)  
plot_heatmap(top20OTU.names, sample.label="Copartment")  


  
  
