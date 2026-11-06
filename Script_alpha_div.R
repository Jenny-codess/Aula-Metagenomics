### The Presentation
##  My task is on alpha diversity on Stool and Tissue


# ---- Verificar e carregar os pacotes ----

required_packages <- c(
  "tidyverse",
  "phyloseq",
  "vegan",
  "ggpubr",
  "cowplot",
  "patchwork",
  "ape",
  "picante",
  "compositions",
  "pheatmap",
  "DESeq2",
  "ANCOMBC",
  "Maaslin2",
  "scater",
  "pairwiseAdonis"
)

# Carregar os pacotes
lapply(
  required_packages,
  library,
  character.only = TRUE
)

# Estabelecer o diretório de trabalho
results_folder <- "C:/Users/hp/OneDrive/Desktop/Aula Metagenomica/Aula Metagenomia/results"

##################################
# -----   Importar dados   ----- #
##################################

tabASV <- read.table("feature-table-crc.tsv", header = T, check.names = F, row.names = 1)     #quantification table (matrix)
tabASV <- as.matrix(tabASV)
head(tabASV)[1:5,1:5]

tabTax <- read.delim("taxonomy_crc.txt", header = T, row.names = 1)           #taxonomy table
tabTax <- as.matrix(tabTax)
head(tabTax)[1:5,1:5]

Arv <- read_tree("tree.nwk")       #tree network

metadataTab <- read.csv("crc_metadata.csv", header = T, row.names = 1, sep = ",")      #metadata
View(metadataTab)

# Create the tissue sample dataframe
tissue_metadata <- metadataTab %>% 
  filter(str_detect(SampleType2, "Tissue"))

# Create the stool sample dataframe
stool_metadata <- metadataTab %>% 
  filter(str_detect(SampleType2, "Stool"))

##################################
# --- Importar para phyloseq --- #
##################################

ASV <- otu_table(tabASV,taxa_are_rows = TRUE)
TAX <- tax_table(tabTax)
META <- sample_data(metadataTab)
TREE <- Arv

META_tissue <- sample_data(tissue_metadata)
META_stool <- sample_data(stool_metadata)

ps_tissue <- phyloseq(ASV, TAX, META_tissue, TREE)
ps_stool <- phyloseq(ASV, TAX, META_stool, TREE)

#Check the phylosec object
ps_tissue
ps_stool


# Verify the classes of each variable
variable_classes_tissue <- sapply(META_tissue, class)
variable_classes_stool <- sapply(META_stool, class)


# Explore o phylosec object (using ps_tissue for the exploration)

nsamples(ps_tissue)     # No of samples
ntaxa(ps_tissue)        # No of ASVs
rank_names(ps_tissue)   # Taxonomy levels
sample_variables(ps_tissue)    # What metadata
sample_sums(ps_tissue)         # How many reads per sample

nsamples(ps_stool)     # No of samples
ntaxa(ps_stool)        # No of ASVs
rank_names(ps_stool)   # Taxonomy levels
sample_variables(ps_stool)    # What metadata
sample_sums(ps_stool)         # How many reads per sample

# Distribution of reads per ASV

# Create a summary table with total read count per ASV
tabelaContagem_tissue <- data.frame(ContagemTotal_tissue = taxa_sums(ps_tissue))
head(tabelaContagem_tissue)

tabelaContagem_stool <- data.frame(ContagemTotal_stool = taxa_sums(ps_stool))
head(tabelaContagem_stool)


# Use `ggplot2` to visualize the distribution of these reads with a histogram
ggplot(tabelaContagem_tissue, aes(ContagemTotal_tissue)) + geom_histogram()

ggplot(tabelaContagem_stool, aes(ContagemTotal_stool)) + geom_histogram() 


############################
# --- Alfa diversity --- #
############################

amostras_tissue <- tabASV[, colnames(tabASV) %in% rownames(tissue_metadata)]
amostras_stool <- tabASV[, colnames(tabASV) %in% rownames(stool_metadata)]

random_sample_tissue <- sample(colnames(amostras_tissue), size = 1)
random_sample_tissue

random_sample_stool <- sample(colnames(amostras_stool), size = 1)
random_sample_stool


# Calculate the total number of reads in the random sample to use as the denominator for proportions.
total_reads_tissue <- sum(amostras_tissue[, random_sample_tissue])
total_reads_stool <- sum(amostras_stool[, random_sample_stool])


# Calculate the proportion of reads assigned to each ASV in the random sample (Tissue)
proportions_tissue <- amostras_tissue[, random_sample_tissue] / total_reads_tissue

observed_richness_manual_tissue <- sum(amostras_tissue[, random_sample_tissue] > 0)
observed_richness_manual_tissue

shannon_manual_tissue <- -sum(proportions_tissue[proportions_tissue > 0] * log(proportions_tissue[proportions_tissue > 0]))
shannon_manual_tissue

singleton_tissue <- sum(amostras_tissue[, random_sample_tissue] == 1) # Number of singletons
singleton_tissue

doubleton_tissue <- sum(amostras_tissue[, random_sample_tissue] == 2) # Number of doubletons
doubleton_tissue

chao1_manual_tissue <- observed_richness_manual_tissue + (singleton_tissue^2 / (2 * max(doubleton_tissue, 1)))
chao1_manual_tissue


# Calculate the proportion of reads assigned to each ASV in the random sample (Stool)
proportions_stool <- amostras_stool[, random_sample_stool] / total_reads_stool

observed_richness_manual_stool <- sum(amostras_stool[, random_sample_stool] > 0)
observed_richness_manual_stool

shannon_manual_stool <- -sum(proportions_stool[proportions_stool > 0] * log(proportions_stool[proportions_stool > 0]))
shannon_manual_stool

singleton_stool <- sum(amostras_stool[, random_sample_stool] == 1) # Number of singletons
singleton_stool

doubleton_stool <- sum(amostras_stool[, random_sample_stool] == 2) # Number of doubletons
doubleton_stool

chao1_manual_stool <- observed_richness_manual_stool + (singleton_stool^2 / (2 * max(doubleton_stool, 1)))
chao1_manual_stool


# Estimate richness metrics available in IN 'PHYLOSEQ' package using our ps package
alpha_diversity_tissue <- estimate_richness(ps_tissue)
head(alpha_diversity_tissue)[1:4,1:9]
View(alpha_diversity_tissue)

alpha_diversity_stool <- estimate_richness(ps_stool)
head(alpha_diversity_stool)[1:4,1:9]
View(alpha_diversity_stool)


#Let's see how our manual sample results compare to our manually calculated results. Is it equal?
alpha_diversity_tissue$amostras_tissue <- rownames(alpha_diversity_tissue)
head(alpha_diversity_tissue)
subset(alpha_diversity_tissue, amostras_tissue == random_sample_tissue)

alpha_diversity_stool$amostras_stool <- rownames(alpha_diversity_stool)
head(alpha_diversity_stool)
subset(alpha_diversity_stool, amostras_stool == random_sample_stool)

#Convert our sample data to dataframe
metadata_tissue <- sample_data(ps_tissue)
metadata_stool <- sample_data(ps_stool)


alpha_div_meta_tissue <- cbind(metadata_tissue, alpha_diversity_tissue)
alpha_div_meta_stool <- cbind(metadata_stool, alpha_diversity_stool)

head(alpha_div_meta_tissue)
View(alpha_div_meta_tissue)

# Plot Shannon diversity across urbanism groups
#p_Shannon <- ggplot(alpha_div_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
#geom_boxplot() 

#p_Shannon


#The beginning of Wednesday class

# Now let's visualize the alpha diversity differences between the different classes in our project.
# Choose a Metric and a CATEGORICAL Variable to generate the image. Then, test other metrics and other variables.
# Example: ggplot(alpha_div_meta, aes(x = urbanism, y = Shannon, fill = urbanism))

# Test with Shannon variable - Tissue
p_Metrica_tissue <- ggplot(alpha_div_meta_tissue, aes(x = SampleType2, y = Shannon, fill = SampleType2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  ggtitle("Alpha Diversity: Shannon Metric vs tissue sampletypes.") +
  theme_minimal()


p_Metrica_tissue
ggsave(filename = file.path(results_folder, "shannon_tissue.png"), plot = p_Metrica_tissue)

# Tissue observed
p_Metrica_tissue2 <- ggplot(alpha_div_meta_tissue, aes(x = SampleType2, y = Observed, fill = SampleType2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  ggtitle("Alpha Diversity: Observed Metric vs tissue sampletypes.") +
  theme_minimal()


p_Metrica_tissue2
ggsave(filename = file.path(results_folder, "observed_tissue.png"), plot = p_Metrica_tissue2)

ggsave(filename = file.path(results_folder, "Alpha_div_tissue_metrics.png"), plot = p_Metrica_tissue/p_Metrica_tissue2)


# Now test with a Observed variable: - Stool
p_Metrica_stool <- ggplot(alpha_div_meta_stool, aes(x = SampleType2, y = Shannon, fill = SampleType2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  ggtitle("Alpha Diversity: Shannon Metric vs stool sampletypes.") +
  theme_minimal()

p_Metrica_stool
ggsave(filename = file.path(results_folder, "shannon_stool.png"), plot = p_Metrica_stool)


# Stool Observed
p_Metrica_stool2 <- ggplot(alpha_div_meta_stool, aes(x = SampleType2, y = Observed, fill = SampleType2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  ggtitle("Alpha Diversity: Observed Metric vs stool sampletypes.") +
  theme_minimal()


p_Metrica_stool2
p_Metrica_stool/p_Metrica_stool2

ggsave(filename = file.path(results_folder, "observed_stool.png"), plot = p_Metrica_stool2)
ggsave(filename = file.path(results_folder, "Alpha_div_stool_metrics.png"), plot = p_Metrica_stool/p_Metrica_stool2)



# With the "spicy" package, calculate alpha diversity using the Faith method, which considers phylogenetic distances. Did you observe any difference?
install.packages("picante")
library(picante)
otu_matrix_tissue <- as(otu_table(ps_tissue), "matrix")

if (taxa_are_rows(ps_tissue)) {
  otu_matrix_tissue <- t(otu_matrix_tissue)
}

phy_tree_tissue <- phy_tree(ps_tissue)

# Calculate Faith's Phylogenetic Diversity (PD)
faith_pd_tissue <- pd(otu_matrix_tissue, phy_tree_tissue)

# view the result
faith_pd_tissue

# # Add Faith's metrics to metasmetadas:
alpha_div_meta_tissue2 <- cbind(alpha_div_meta_tissue, faith_pd_tissue)
head(alpha_div_meta_tissue2)


# Visualize:
p_PD_variavel_tissue <- ggplot(alpha_div_meta_tissue2, aes(x = SampleType2, y = PD, fill=SampleType2)) +
  geom_boxplot() 
p_PD_variavel_tissue



### Rarefaction
# Will our results be different if we standardize the number of reads per sample?
# Visualize a rarefaction curve for your data:
rarecurve(t(amostras_tissue), step = 200, cex = 0.1)
rarecurve(t(amostras_stool), step = 200, cex = 0.1)

# Need to zoom in to establish an appropriate cutoff number? Replace 20000 with the maximum number of reads you want to view.
rarecurve(t(amostras_tissue), step = 200, cex = 0.6, xlim=c(0,20000))
rarecurve(t(amostras_stool), step = 200, cex = 0.6, xlim=c(0,20000))

# Rarefy your data. Replace "chosen depth" with the number of reads you want to rarefy:
ps.rare_tissue = rarefy_even_depth(ps_tissue, rngseed=1, sample.size=5000, replace=F)   # Using 5000 samples, the loss is only 3
ps.rare_stool = rarefy_even_depth(ps_stool, rngseed=1, sample.size=5000, replace=F)   # Using 5000 samples, the loss is only 4


# Calculate alpha diversity for the rarefied samples.
alpha_diversity_rare_tissue <- estimate_richness(ps.rare_tissue)
alpha_diversity_rare_stool <- estimate_richness(ps.rare_stool)

rownames(alpha_diversity_rare_tissue) <- gsub("^X", "", rownames(alpha_diversity_rare_tissue))
alpha_diversity_rare_tissue$sample <- rownames(alpha_diversity_rare_tissue)

rownames(alpha_diversity_rare_stool) <- gsub("^X", "", rownames(alpha_diversity_rare_stool))
alpha_diversity_rare_stool$sample <- rownames(alpha_diversity_rare_stool)

head(alpha_diversity_rare_tissue)
head(alpha_diversity_rare_stool)

# Combine the new metrics (with rarefaction and their metadata)
#alpha_div_rare_meta <- cbind(metadata, alpha_diversity_rare)


metadata_tissue <- sample_data(ps.rare_tissue) %>% data.frame()
alpha_div_rare_meta_tissue <- metadata_tissue %>% rownames_to_column("sample") %>%
  left_join(alpha_diversity_rare_tissue, by = "sample")

View(metadata_tissue)
View(alpha_div_rare_meta_tissue)

metadata_stool <- sample_data(ps.rare_stool) %>% data.frame()
alpha_div_rare_meta_stool <- metadata_stool %>% rownames_to_column("sample") %>%
  left_join(alpha_diversity_rare_stool, by = "sample")

View(metadata_stool)
View(alpha_div_rare_meta_stool)


# Redo your plots and observe if there are differences in the results (without statistics yet).
p_Metrica_variavel_rare_tissue <- ggplot(alpha_div_rare_meta_tissue, aes(x = SampleType2, y = Shannon, fill = SampleType)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_minimal() +
  ggtitle("With Shannon rarefaction reads per sample") +
  theme(legend.position = "none")

p_Metrica_variavel_rare_tissue
p_Metrica_variavel_rare_tissue/p_Metrica_tissue


p_Metrica_variavel_rare_stool <- ggplot(alpha_div_rare_meta_stool, aes(x = SampleType2, y = Shannon, fill = SampleType)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_minimal() +
  ggtitle("With Shannon rarefaction reads per sample") +
  theme(legend.position = "none")

p_Metrica_variavel_rare_stool
p_Metrica_variavel_rare_stool/p_Metrica_stool




#--------------------------------------------------------------------------------------------

shapiro.test(alpha_div_meta_tissue$Shannon)
shapiro.test(alpha_div_meta_stool$Shannon)
shapiro.test(alpha_div_meta$Shannon)

# Non parametric: For 2 groups, use a Wilcoxon rank-sum test (also known as the Mann-Whitney U test)
wilcox.test(Shannon ~ SampleType2, data = alpha_div_meta_tissue)
wilcox.test(Observed ~ SampleType2, data = alpha_div_meta_tissue)

# For multiple groups, use Kruskal-Wallis test
kruskal.test(Shannon ~ SampleType2, data = alpha_div_meta_stool)
kruskal.test(Observed ~ SampleType2, data = alpha_div_meta_stool)








