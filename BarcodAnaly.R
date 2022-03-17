####Microbiome analysis using the QIIME2R package and others

##Installing the dependencies 

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
devtools::install_github("jbisanz/MicrobeR")
BiocManager::install("DECIPHER")
BiocManager::install("philr")
BiocManager::install("microbiome", force = TRUE)
install.packages("remotes")
remotes::install_github("microsud/microbiomeutilities")
if (!require(remotes)) install.packages("remotes")
remotes::install_github("yiluheihei/microbiomeMarker")
BiocManager::install("metagenomeFeatures")
install_github("umerijaz/microbiomeSeq") 
BiocManager::install("DESeq2")
BiocManager::install("limma")
install.packages("randomForest")
install.packages("relocate")


##Loading libraries, if they are not installed in R, look in packages-install. You need to load the packages every time you start R.

library(qiime2R)
library(tidyverse)
library(MicrobeR)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(ggplot2)
library(ggpubr)
library(plotly)
library(pheatmap)
library(metagenomeFeatures)
library(tidytree)
library(data.table)
library(vegan)
library(RColorBrewer)
library(plyr)
library(ape)
library(DT)
library(metacoder)
library(adespatial)
library(microbiomeSeq)
library(ggdendro)
library(DESeq2)
library(jpeg)
library(gridExtra)
library(knitr)
library(pander)
library(png)
library(randomForest)
library(rpart)
library(magrittr)
library(scales)
library(reshape2)
library(dplyr)
library(microDecon)
library(writexl)
library(relocate)
library(tidyr)
###tutorial: https://uw-madison-microbiome-hub.github.io/Microbiome_analysis_in-_R/

####Starting metabarcoding analysis####
# Importing the ASV abundance file - Here it can be the qza file from the dada2-table.qza

ASVs <- read_qza("otus_descontaminadas.qza")

####In silico descontamination####

##Create an object to save the spreadsheet that will be analyzed. the header keeps the table format separated by columns.
##Transform the output file from QIIME genus.qza to a .tsv file, later saved in .csv to use in this command
##The first column of the file must contain OTU_ID and OTU in sequential numbers, the taxa information is in the last column

input_decontam <- ASVs$data

# n is the total number of ASVs

###n is variable according to your data
n <- 4927
prefix <- "OTU"
suffix <- seq(1:n)
OTU <- paste(prefix,suffix,sep='')
OTU <- as.data.frame(OTU)

input2 <- cbind.data.frame(OTU,input_decontam)

input2

# move the taxa column to last column and name it as "taxa".

taxa <- rownames(input2)
rownames(input2) <- NULL
input3 <- cbind(taxa,input2)
input3

# It is necessary to provide the name of the colunm in .after

input4 <- input3 %>% relocate (taxa, .after = "XXX")
input4

### data is the object previously created with the sheet that will be read, numb.blanks is the number of blanks in the sheet
### numb.ind is the number of individuals in the sheet, rate is T or F to indicate if the sheet has a final column with taxa. T=true and F=false
results1 <- decon(data = input4, numb.blanks = 1, numb.ind = 30, taxa = T)

##command used to view the table that has been created, you can also open them directly in the R console

results1$decon.table

results1$OTUs.removed

### how to export the data out of R and save to a folder

write.csv(results1, file = "descontaminatedata.csv")

decon_data <- results1$decon.table

# Importing metadata, the same file used in the qiime2 analysis
metadata <- read.table("metadata.tsv", , sep='\t', header=T, row.names=1, comment="")

metadata <- metadata[-1,] # remove the second line that specifies the data type

# Importing tree
tree <- read_qza("rooted_tree.qza")

# Importing taxonomy
####import the taxonomy excel formatted file by import dataset.To get this file just unzip the taxonomy.qza until you find the file taxonomy.tsv. Then open excel and adjust the columns.
#import dataset -> from excel in R environment 
a_tab = as.matrix(input_taxonomy)

head(a_tab)

summary(a_tab)

str(a_tab)

is.na(a_tab2)

taxonomy_clean <- na.omit(a_tab2)

# Creating a phyloseq object 
physeq <- phyloseq(
  otu_table(ASVs$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(taxonomy_clean),
  sample_data(metadata)
)

# Checking that the files were built correctly
summarize_phyloseq(physeq)
print_ps(physeq)
summary(sample_sums(physeq))

ntaxa(physeq)

nsamples(physeq)

sample_names(physeq)[1:5]  

rank_names(physeq)  

sample_variables(physeq)  

otu_table(physeq)[1:5, 1:5]  

tax_table(physeq)[1:5, 1:4]


###plot of read/abundance distribution according to the characteristic identified in "groups
plot_read_distribution(physeq, groups = "Groups", plot.type = "density") + theme_biome_utils()
ggsave("Distribuição_reads.png")

##Generates a rarefaction file, an error appears, but according to the internet you just ignore it and the file is generated the same
physeq_rarefy <- rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)

print(physeq_rarefy)

##removing OTUs that are classified as unnasigned from the kingdom

# check if there are no OTUs present in any samples
any(taxa_sums(physeq) == 0)

ps1a <- prune_taxa(taxa_sums(physeq) > 0, physeq)

# double-check that no OTUs are present in any sample
any(taxa_sums(ps1a) == 0)

# subtract the number of OTUs in the original (ps1) from the number of OTUs in the new phyloseq object (ps1a)
# number of OTUs in the original 
ntaxa(physeq)

# number of OTUS in the new file
ntaxa(ps1a)

ntaxa(physeq) - ntaxa(ps1a)

#Check taxonomic classification information 
rank_names(physeq) 

# the table is interactive and you can search it for details.
datatable(tax_table(physeq)) 

# Double-check the taxonomy levels
rank_names(ps1a)

# Use the interactive table to check under which taxonomic classification Chloroplast and Mitrochondria are mentioned. Specify Class or Order in the code below.    
ps1a <- subset_taxa(physeq,Class!="c__Chloroplast")

ps1b <- subset_taxa(ps1a,Order!="o__Mitochondria")

ps1c <- subset_taxa(ps1b,Kingdom!="d__Eukaryota")


# Second question: Do we expect or wish to include Archea sequences? If you used primers that do not target them for Archea amplification, then it is a spurious observation.
#Keep in mind: The fact that you took these sequences with any cut-off from your otupicking algorithm means that they were relatively abundant, check the confidence of the sequences and taxonomic assignments before you decide to remove them.
#If you find high abundance/high confidence sequences that in theory shouldn't be there, maybe more things are wrong with your data. Be critical. OR, your primers got these sequences.

ps1 <- subset_taxa(ps1c,Kingdom!="d__Archaea")

sort(sample_sums(ps1))  #mark the number of reads/sample, we can see that one of the samples has only 1867 seqs, normally we would remove it, but since it is important for data analysis, we can check if this sequencing depth is enough to capture the theoretical composition, so for now we will leave it in. Just pay extra attention to where this sample ends up in the downstream analysis.


# However, in case you want to remove some samples due to low sequencing depth or some other technical reasons, you can remove them from the analysis with the following command:
# ps1.sub <-subset_samples(ps1, #phyloseq object
# sampleID !="SampleName") #metadata category name and != "SampleName"

# ps1 <- ps1.sub # this is our filtered file 

###If you are not removing anything, follow with the command. Evaluating the table interactively
datatable(tax_table(ps1))

physeq_rarefy_cleaned <- rarefy_even_depth(ps1, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)


### creating variables that will be used later

physeq.f <- format_to_besthit(ps1)

comps <- make_pairs(sample_data(physeq.f)$Groups)


####Using the last generated object with the removal of unwanted AVS 

plot_abundance_phylum <- plot_taxa_prevalence(ps1, "Phylum") 

# Bacterial names in italics
plot.map <- plot_abundance_phylum + theme(legend.text = element_text(colour = 'black', size = 10, face = 'italic'))

print(plot.map)

ggsave("abundance_phylum.png", width = 14, height = 10, dpi = 600)

physeq_fam <- microbiome::aggregate_top_taxa(ps1, "Family", top = 10)

physeq_gen <- microbiome::aggregate_top_taxa(ps1, "Genus", top = 10)

physeq.gen.rel <- microbiome::transform(physeq_gen, "compositional")

physeq.fam.rel <- microbiome::transform(physeq_fam, "compositional")

plot_composition(physeq.fam.rel,sample.sort = "Groups", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))
ggsave("family_groups.png")

plot_composition(physeq.gen.rel,sample.sort = "Groups", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Genus", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))
ggsave("family_groups.png")

taxa_barplot(summarize_taxa(ASVs$data, as.data.frame(taxonomy_clean))$Genus, metadata, "Groups")

# Make it interactive
ggplotly(taxa_barplot(summarize_taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Groups"))

# Save the plot
b.plot <- taxa_barplot(summarize_taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Groups")

ggsave("barplot_family.png", b.plot,  width = 14, height = 10, dpi = 300)

taxa_heatmap(summarize_taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Groups")


physeq_df <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                             transform.counts = "compositional")

# An additonal Sam_rep column with sample names is created for reference purposes
colnames(physeq_df)

# Box plot at family level
ggstripchart(physeq_df, "Groups", "Abundance", 
             facet.by = "Family", color = "Groups", facet = "italic",
             palette = "jco") + rremove("x.text")

ggsave("blot_family.png")

### saving a variable with the colors that will be used
mycols <- c("brown3", "steelblue", "grey50")

###Filtering what is not to appear in the image 
prevelancedf = apply(X = otu_table(ps1),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# Add taxonomy and total read count to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(ps1),
                          tax_table(ps1))
prevelancedf[1:10,]

plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

phyla2Filter = c("p__Abditibacteriota", "p__Acidobacteriota", "p__Actinobacteriota",
                 "p__Aquificota","p__Armatimonadota", "p__Fusobacteriota",
                 "p__Deinococcota", "p__Campilobacterota", "p__Chloroflexi", "p__Cyanobacteria")

# Filter entries with unidentified Phylum.
ps1_new = subset_taxa(ps1, !Phylum %in% phyla2Filter)
ps1_new

boxplot_family <- plot_taxa_boxplot(ps1_new,
                  taxonomic.level = "Phylum",
                  top.otu = 6, 
                  group = "Groups",
                  add.violin= FALSE,
                  title = "Top six family", 
                  keep.other = FALSE,
                  group.order = c("CVN", "CVP", "NSA"),
                  group.colors = mycols,
                  dot.size = 2) + theme_biome_utils()

estati_family <- boxplot_family + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.04,
  method = "wilcox.test")

estati_family + theme(axis.text = element_text(size = 17)) + theme(plot.title = element_text(size = rel(2)))

tiff("boxplot_family_estat.tiff",height = 13, width = 16, 
     units = "in", res=600)


estati_family

dev.off()


####Filtering the genera
generaFilter = c("g__Mycoplasma", "g__Streptobacillus", "g__Histophilus")

# Filter entries with unidentified Phylum.
ps1_new_g = subset_taxa(ps1, !Genus %in% generaFilter)
ps1_new_g



boxplot_genus <- plot_taxa_boxplot(ps1_new_g,
                                    taxonomic.level = "Genus",
                                    top.otu = 9, 
                                    group = "Groups",
                                    add.violin= FALSE,
                                    title = "Top six genus", 
                                    keep.other = FALSE,
                                    group.order = c("CVN", "CVP", "NSA"),
                                    group.colors = mycols,
                                    dot.size = 2) + theme_biome_utils() 


estati_genus <- boxplot_genus + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.04,
  method = "wilcox.test")

estati_genus + theme(axis.text = element_text(size = 17)) + theme(plot.title = element_text(size = rel(2)))


tiff("boxplot_genus_estat.tiff",height = 13, width = 16, 
     units = "in", res=600)

estati_genus

dev.off()


library(ggplot2)

####TOP GENERA
physeq.genus <- aggregate_taxa(ps1, "Genus")
top_four <- top_taxa(physeq.genus, 10)
top_four


physeq.family <- aggregate_taxa(ps1, "Family")
top_four_fam <- top_taxa(physeq.family, 10)
top_four_fam

mycols <- c("brown3", "steelblue", "grey50")

top_genera_fam <- plot_listed_taxa(physeq.family, top_four_fam, 
                               group= "Groups",
                               group.order = c("CNV", "CVP", "NSA"),
                               group.colors = mycols,
                               add.violin = F,
                               dot.opacity = 0.25,
                               box.opacity = 0.25,
                               panel.arrange= "wrap"
)

top_genera_fam


top_genera <- plot_listed_taxa(physeq.genus, top_four, 
                               group= "Groups",
                               group.order = c("CNV", "CVP", "NSA"),
                               group.colors = mycols,
                               add.violin = F,
                               dot.opacity = 0.25,
                               box.opacity = 0.25,
                               panel.arrange= "wrap"
)

tiff("boxplot_genera.tiff",height = 8, width = 14, 
     units = "in", res=600)

top_genera

dev.off()


boxplot <- grid.arrange(top_genera_fam, top_genera, nrow = 1)

tiff("boxplot_arrange.tiff",height = 10, width = 17, 
     units = "in", res=600)

boxplot

dev.off()

top_genera + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")


##Dominant taxa

physeq.gen <- aggregate_taxa(ps1,"Genus")

dom.tax <- dominant_taxa(ps1,level = "Genus", group="Groups")

head(dom.tax$dominant_overview)


##summary taxa

taxa_summary(ps1, "Phylum")


##Summarizing by groups
grp_abund <- get_group_abundances(ps1, 
                                  level = "Phylum", 
                                  group="Groups",
                                  transform = "compositional")


# clearing the names
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)

mean.plot <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=Groups)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("Groups", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

mean.plot

#Creating a heatmap

#You can use the input_taxonomy file or a table with the data you want to show in the Heatmap
heatmap_generos <- scale(Tabela_generos_heatmap)

#install.packages("pheatmap")
library("pheatmap")

#simple
pheatmap(heatmap_generos)


#heatmap cut into 4 groups (dendogram result)
pheatmap(heatmap_generos, cutree_rows = 4,scale = "none")

#normalized by column - Pass the original df
pheatmap(heatmap_generos, cutree_rows = 4,scale = "row")

#no column and row clustering (original entry)
pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
         cluster_rows = FALSE,cluster_cols = FALSE)

#clustering na coluna
pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
         cluster_rows = FALSE,cluster_cols = TRUE)


#install.packages("viridis")
library(viridis)

pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
         cluster_rows = TRUE,cluster_cols = FALSE,
         color=inferno(20))

pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
         cluster_rows = TRUE,cluster_cols = FALSE,
         color=viridis(15))

## coloring the groups (need to create a dataframe) ##factor are the names of your table's characteristics that will be used
grupos = data.frame(ID = factor(rep(c("CNV","CVP","NSA"), each=10)))
rownames(grupos)<-row.names(heatmap_generos)

cores = list(ID = c(CNV="#80CDC1", CVP="#FDB863",NSA="#92C5DE"))


newnames <- lapply(
  colnames(heatmap_generos),
  function(x) bquote(italic(.(x))))

heatmap_final <- pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
                          cluster_rows = TRUE,cluster_cols = FALSE,
                          annotation_row  = grupos,annotation_colors=cores[1],
                          annotation_legend = TRUE,
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          clustering_method = "complete",
                          cellwidth = 10,
                          labels_col = as.expression(newnames),
                          color=viridis(15))

tiff("heatmap_genera.tiff",height = 10, width = 28, units = 'in',res=600)

heatmap_final

dev.off()


# creating with more than one class
annotation_row = data.frame(
  grupo = factor(rep(c("CNV","CVP","NSA"), each=10)),
  adicional = c(rep(c("random1","random2"),each=25))
)

rownames(annotation_row) = row.names(heatmap_generos)

pheatmap(heatmap_generos, cutree_rows = 4,scale = "row",
         cluster_rows = TRUE,cluster_cols = FALSE,
         color=viridis(10),
         annotation_row  = annotation_row,
         annotation_legend = TRUE)


# Saving the plot
ggsave("heatmap_family.png", h.map,  width = 14, height = 10, dpi = 300)


####Alpha diversity####

###ANOVA, t-test, or general linear models with normal distribution are used when data are more or less normal
###Kruskal-Wallis, Wilcoxon rank sum test, or general linear models with another distribution are used when the data is not normal
library(ranacapa)

rarefacao <- ggrare(ps1, step = 50, color="Groups", label = "Sample", se = TRUE) 

tiff("ggrare_plot.tiff",height = 14, width = 14, 
     units = 'in',res=600)

rarefacao

dev.off()  

alpha_sample <- plot_richness(physeq_rarefy_cleaned, measures="Shannon")

print(alpha_sample)

a.div <- plot_richness(physeq_rarefy_cleaned, x="Groups", measures=c("Shannon", "simpson", "Observed"), color = "Groups") + geom_boxplot() + theme_bw()

tiff("rarefation.tiff",height = 14, width = 14, 
     units = 'in',res=600)

a.div

dev.off()

# Adding Statistical Support
alfa <- a.div + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")

tiff("alfa.tiff",height = 5, width = 8, units = 'in',res=600)

alfa

dev.off()



##Creating a .csv with the results
richness <- estimate_richness(physeq_rarefy_cleaned,measures=c("Shannon", "simpson", "Observed"))

write.csv(richness, file = "alpha_div.csv")


###graphic with only one diversity
plot_diversity_stats(ps1, group = "Groups", 
                     index = "diversity_shannon", 
                     group.order = c("CVP","CNV","NSA"),                      
                     group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE) + ylab("Shannon Diversity") + xlab("")

ggsave("shannon_groups.png", width = 14, height = 10, dpi = 300 )



plot_anova_diversity(ps1, method = c("richness", "simpson", "shannon"), grouping_column = "Groups", pValueCutoff = 0.05)


###shapiro Wilk for normality test
shannon <- estimate_richness(ps1,measures=c("Shannon"))

shapiro.test(shannon$Shannon)

hist(shannon$Shannon, main="Shannon diversity", xlab="", breaks=10)


### If the comparison data has only two variables, the t test... More uses ANOVA
###We use Kruskal-Wallis (non-parametric equivalent of ANOVA). If we only have two levels

#### kruskal-wallis
##kruskal.test(shannon$Shannon ~ campy, data=metadata)
####pairwise.wilcox.test(shannon$Shannon, metadata$campy, p.adjust.method="fdr")

#Run the ANOVA and save it as an object
aov.shannon.campy <- aov(shannon$Shannon ~ Groups, data=metadata)

#Call for the summary of this ANOVA, which will include the P-values
summary(aov.shannon.campy)

###To do all pairwise comparisons between groups and correct for multiple comparisons, we do our Tukey's ANOVA significance test
TukeyHSD(aov.shannon.campy)


####Beta diversity####

library(MESS)
library(ggrepel)

####beta non-phylogenetic diversity


physeq.ord <- ordinate(physeq_rel, "PCoA", "bray")

physeq.ord

b.div.bray <- plot_ordination(physeq_rel, physeq.ord, type= "samples", axes = 1:2, color= "Groups") + geom_point(size=4) + theme_classic() + scale_color_brewer("Groups", palette = "Set2")


tiff("bray_curtis.tiff",height = 5, width = 8, units = 'in',res=600)

b.div.bray

dev.off()

abrel_bray <- phyloseq::distance(physeq_rel, method = "bray")

abrel_bray <- as.matrix(abrel_bray)

head(abrel_bray)[,1:6]


relab_genera = transform_sample_counts(physeq_rel, function(x) x / sum(x) * 100)


samples <- data.frame(sample_data(relab_genera))

adonis(abrel_bray ~ Groups, data = samples)

permanova <- adonis2(abrel_bray ~ Groups, data = samples)

permanova

set.seed(36) #reproducible results

groups_bray.div <-adonis2(abrel_bray ~ Groups, data = samples, permutations = 999, method="bray")

groups_bray.div


###high phylogenetic diversity

# convert relative abundance
physeq_rel <- microbiome::transform(ps1, "compositional")

physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted=T)

b.div.wuni <- plot_ordination(physeq_rel, physeq.ord.wuni, type= "samples", color= "Groups") + geom_point(size=3)

b.div.wuni <- b.div.wuni + stat_ellipse() + ggtitle("Weighted Unifrac")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")

print(b.div.wuni)
ggsave("weighted_unifrac.png", width = 14, height = 10, dpi = 300)

###Permanova 

w.unifrac.dist <- UniFrac(physeq_rel, 
                          weighted = TRUE, 
                          normalized = TRUE,  
                          parallel = FALSE, 
                          fast = TRUE)

permanova <- adonis2(w.unifrac.dist ~ Groups, data = metadata)


permanova

#####Homogeneity condition check
physeq.disper <- betadisper(w.unifrac.dist, metadata$Groups)

permutest(physeq.disper, pairwise = TRUE)


###########Core microbiota for only one group
physeq.CVP <- subset_samples(ps1, Groups == "CVP")

# convert to relative abundance
physeq.CVP.rel <- microbiome::transform(physeq.CVP, "compositional")

physeq.CVP.rel2 <- prune_taxa(taxa_sums(physeq.CVP.rel) > 0, physeq.CVP.rel)

core.taxa.standard <- core_members(physeq.CVP.rel2, detection = 0.001, prevalence = 50/100)

print(core.taxa.standard)

# Extracting the taxonomy table
taxonomy_core <- as.data.frame(tax_table(physeq.CVP.rel2))

# Taxonomy table to include only major OTUs
core_taxa_id <- subset(taxonomy_core, rownames(taxonomy_core) %in% core.taxa.standard)

DT::datatable(core_taxa_id)


########A total core abundance in each sample (sum of the abundances of the core members):

core.abundance <- sample_sums(core(physeq.CVP.rel2, detection = 0.001, prevalence = 50/100))

DT::datatable(as.data.frame(core.abundance))


### AVS core visualization - Heatmap

# Core with compositions:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define the color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(physeq.CVP.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

tiff("p.core.tiff",height = 8, width = 7, units = 'in',res=600)

p.core

dev.off()

##Use the microbiomeutilities format_to_besthit function to get the best classification ASVs

physeq.CVP.rel2.f <- microbiomeutilities::format_to_besthit(physeq.CVP.rel2)

p.core <- plot_core(physeq.CVP.rel2.f, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

p.core + theme(axis.text.y = element_text(face="italic"))

tiff("p.core.tiff",height = 8, width = 7, units = 'in',res=600)

p.core

dev.off()

###Relation of prevalence of abundance
library(ggrepel)

physeq_comp <- microbiome::transform(ps1, "compositional")

# Select a group
physeq_comp_CVP <- subset_samples(physeq_comp, Groups=="CVP")

physeq_comp_CVP <- core(physeq_comp_CVP,detection = 0.0001, prevalence = 0.50) # reduce size for example

physeq_comp_CVP <- format_to_besthit(physeq_comp_CVP)

set.seed(163897)

prevalance <- plot_abund_prev(physeq_comp_CVP, 
                              label.core = TRUE,
                              color = "Phylum", # NA or "blue"
                              mean.abund.thres = 0.01,
                              mean.prev.thres = 0.99,
                              dot.opacity = 0.7,
                              label.size = 3,
                              label.opacity = 1.0,
                              nudge.label=-0.15,
                              bs.iter=999, # increase for actual analysis e.g. 999
                              size = 20, # increase to match your nsamples(asv_ps)
                              replace = TRUE,
                              label.color="#5f0f40") 
prevelance <- prevalance + 
  geom_vline(xintercept = 0.95, lty="dashed", alpha=0.7) + 
  geom_hline(yintercept = 0.01,lty="dashed", alpha=0.7) +
  scale_color_brewer(palette = "Dark2")



tiff("prevalence.tiff",height = 8, width = 7, units = 'in',res=600)

prevelance

dev.off()

###Tree

physeq_top_50 <- subset_taxa(ps1, Kingdom=="d__Bacteria")

physeq_top_50 <- prune_taxa(names(sort(taxa_sums(physeq_top_50),TRUE)[1:50]), physeq_top_50)

#phylogenetic simples
plot_tree(physeq_top_50)

# Add genera labels to the tree and bootstrap values
plot_tree(physeq_top_50, label.tips="Genus", ladderize="left")

# removing bootstrap values
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left")

# Colorize the nodes by category
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="Groups")

# Add size by abundance
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="Groups", size="abundance")
ggsave("filogenia_abundace_genus.png", width = 14, height = 10, dpi = 300)

# Convert to radial tree
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="Groups") + coord_polar(theta="y")


####Heirachial grouping to visualize the distance between samples using Weighted Unifrac and UPGMA
phy.hclust <- hclust(UniFrac(physeq_rarefy_cleaned, weighted = TRUE), method="average")

ggdendrogram(phy.hclust, rotate = TRUE, theme_dendro = TRUE)

####microbiome network

plot_net(physeq_rel, maxdist = 0.8, color = "Groups")

#Change the distance to Jaccard
plot_net(physeq_rel, maxdist = 0.8, color = "Groups", distance="jaccard")

###igraph-based network
ig <- make_network(physeq_rel, max.dist=0.8)
plot_network(ig, physeq_rel)

# Add color label 
plot_network(ig, physeq_rel, color="Groups", line_weight=0.4, label=NULL)


#Replace the (standard) Jaccard distance method with Bray-Curtis
ig <- make_network(physeq_rel, dist.fun="bray", max.dist=0.8)
plot_network(ig, physeq_rel, color="Groups", line_weight=0.4, label=NULL)
ggsave("network_groups_bray.png", width = 14, height = 10, dpi = 300)

##ramdonForest  tutorial: https://rpubs.com/michberr/randomforestmicrobe

##We then scale our reads normalizing to 15000 reads, which is slightly below our smallest library size

### Normalization ### CREATING A FUNCTION

# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down
# scale_reads <- function(ps1, n) {
#  physeq.scale <-
#    transform_sample_counts(ps1, function(x) {
#      (n * x/sum(x))
#    })
#  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
#  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
#  return(physeq.scale)
#}


## creating variables that will be used later

physeq.f <- format_to_besthit(ps1)

response <- as.factor(sample_data(physeq.f)$Groups)

predictors <- t(otu_table(physeq.f))
dim(predictors)

rf.data <- data.frame(response, predictors)


set.seed(2)
erie.classify <- randomForest(formula = response~., data=rf.data, ntree = 2000, importance=TRUE)
print(erie.classify)
varImpPlot(erie.classify, cex=1.2)


names(erie.classify)

# Make a data table with names of predictors and their importance
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)

head(imp)

# Sort the forecast levels by importance
imp.sort <- arrange(imp, dplyr::desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.20 <- imp.sort[1:20, ]

head(imp.20)

# ggplot
ggplot_ramdon <- ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important AVS for classifying NSA, CNV, CVP")

tiff("ggplot_ramdon.tiff",height = 8, width = 14, units = 'in',res=600)

ggplot_ramdon

dev.off()

otunames <- imp.20$predictors

r <- rownames(tax_table(physeq.f)) %in% otunames


varImpPlot(erie.classify, type=1, pch=19, col=1, cex=.5, main="")

var.plot <- varImpPlot(erie.classify, type=2, pch=19, col="Black", cex=.9, main="")

var.plot

importance(erie.classify, type=1)


##Top Genera

coresBasicas <- c(1:15)

##BarplotNSA


## Input the separate tables for each interest group

genera_NSA <- as.matrix(generos_NSA)

rownames(genera_NSA) <- genera_NSA[,1]

generos_NSA_fi <- genera_NSA[,-1]

box_1 <- barplot(generos_NSA_fi, 
col=c("dodgerblue2", "deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
      "#FDBF6F",
      "gray70", "khaki2",
      "maroon"))


legend("topright",
               legend = c("Pseudomonas", "Ureaplasma", "Stenotrophomonas", "UCG-010", "Mycoplasma", "UCG-005", "Rikenellaceae_RC9_gut_group", "Alistipes", "Bacteroides",
"Streptobacillus", "Histophilus", "Acinetobacter", "Campylobacter", "Bacteroidales_RF16_group", "Prevotellaceae_UCG-004"),
               fill = c("dodgerblue2", "deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
                        "#FDBF6F",
                        "gray70", "khaki2",
                        "maroon"))
        

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

######BarplotCVP

genera_CVP <- as.matrix(generos_CVP)

rownames(genera_CVP) <- genera_CVP[,1]

generos_CVP_fi <- genera_CVP[,-1]

box_2 <- barplot(generos_CVP_fi, 
        col=c("dodgerblue2", "deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
              "#FDBF6F",
              "gray70", "khaki2",
              "maroon"))

tiff("box_2.tiff",height = 13, width = 16, 
     units = "in", res=600)

box_2

dev.off()

legend("topright",
       legend = c("Pseudomonas", "Ureaplasma", "Stenotrophomonas", "UCG-010", "Mycoplasma", "UCG-005", "Rikenellaceae_RC9_gut_group", "Alistipes", "Bacteroides",
                  "Streptobacillus", "Histophilus", "Acinetobacter", "Campylobacter", "Bacteroidales_RF16_group", "Prevotellaceae_UCG-004"),
       fill = c("dodgerblue2", "deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
                "#FDBF6F",
                "gray70", "khaki2",
                "maroon"))

#####Barplot_CNV

genera_CNV <- as.matrix(generos_CNV)

rownames(genera_CNV) <- genera_CNV[,1]

generos_CNV_fi <- genera_CNV[,-1]


box_3 <- barplot(generos_CNV_fi, 
        col=c("dodgerblue2", "deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
              "#FDBF6F",
              "gray70", "khaki2",
              "maroon"))
tiff("box_3.tiff",height = 13, width = 16, 
     units = "in", res=600)

box_3

dev.off()

legend("topright",
       legend = c("Pseudomonas", "Ureaplasma", "Stenotrophomonas", "UCG-010", "Mycoplasma", "UCG-005", "Rikenellaceae_RC9_gut_group", "Alistipes", "Bacteroides",
                  "Streptobacillus", "Histophilus", "Acinetobacter", "Campylobacter", "Bacteroidales_RF16_group", "Prevotellaceae_UCG-004"),
       fill = c("dodgerblue2","deeppink1", "green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6", 
                "#FDBF6F",
                "gray70", "khaki2",
                "maroon"))


####correlation#### 
correlacao <- as.matrix(Supplementary_Table_1)

rownames(correlacao) <- correlacao[,1]

correlacao_clean <- correlacao[,-1]

coorce <- as.data.frame(correlacao_clean)

str(coorce)

sapply(coorce, is.numeric)

coorce[is.na(coorce)] <- 0

table_cor <- as.data.frame(lapply(coorce, as.numeric))

str(table_cor)

core_final <- cor(table_cor, use="pairwise.complete.obs", method = "spearman")

write.csv(core_final, file = "cor_final_spearman.csv")

library(corrr)
install.packages("Hmisc")
library(Hmisc)
install.packages("corrplot")
library(corrplot)

corrplot(core_final, sig.level = 0.01, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)



###Correlation TOP 15 genres
##Manually add the table with the genera that you want to evaluate the correlation
correlacao <- as.matrix(correlacao_top15)

rownames(correlacao) <- correlacao[,1]

correlacao_clean <- correlacao[,-1]

coorce <- as.data.frame(correlacao_clean)

str(coorce)

sapply(coorce, is.numeric)

coorce[is.na(coorce)] <- 0

table_cor <- as.data.frame(lapply(coorce, as.numeric))

str(table_cor)

core_final <- cor(table_cor, use="pairwise.complete.obs", method = "spearman")

core_final

write.csv(core_final, file = "cor_top15.csv")

corplot_paper <- corrplot(core_final, is.corr = FALSE, method = "circle", order= "hclust", type = "upper", 
         tl.col = "black", tl.srt = 45, p.mat = core_final, insig = "p-value", sig.level=-1)

tiff("corplot_teste.tiff",height = 5, width = 8, units = 'in',res=600)

corplot_paper

dev.off()
