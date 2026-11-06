### Disciplina PPGBM/UFPA -- Bioinformática voltada para análises metagenômicas

# Ana Paula Schaan, Setembro 2026


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
setwd("/Users/xxx")
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

##################################
# --- Importar para phyloseq --- #
##################################

ASV <- otu_table(tabASV,taxa_are_rows = TRUE)
TAX <- tax_table(tabTax)
META <- sample_data(metadataTab)
TREE <- Arv

# Create the tissue sample dataframe
tissue_metadata <- metadataTab %>% 
  filter(str_detect(SampleType2, "Tissue"))

# Create the stool sample dataframe
stool_metadata <- metadataTab %>% 
  filter(str_detect(SampleType2, "Stool"))

ps <- phyloseq(ASV, TAX, META, TREE)

#Check the phylosec object
ps


# Verificar as classes de cada variável
variable_classes <- sapply(META, class)
variable_classes

# Explore o phylosec object

nsamples(ps)     # No of samples
ntaxa(ps)        # No of ASVs
rank_names(ps)   # Taxonomy levels
sample_variables(ps)    # What metadata
sample_sums(ps)         # How many reads per sample


# Distribution of reads per ASV

# Create a summary table with total read count per ASV
tabelaContagem <- data.frame(ContagemTotal = taxa_sums(ps))
head(tabelaContagem)

# Use `ggplot2` to visualize the distribution of these reads with a histogram
ggplot(tabelaContagem, aes(ContagemTotal)) + geom_histogram() 

############################
# --- Alfa diversidade --- #
############################

amostras <- tabASV[, colnames(tabASV) %in% rownames(metadataTab)]

random_sample <- sample(colnames(amostras), size = 1)
random_sample


# Calculate the total number of reads in the random sample to use as the denominator for proportions.
total_reads <- sum(amostras[, random_sample])


# Calculate the proportion of reads assigned to each ASV in the random sample
proportions <- amostras[, random_sample] / total_reads

observed_richness_manual <- sum(amostras[, random_sample] > 0)
observed_richness_manual

shannon_manual <- -sum(proportions[proportions > 0] * log(proportions[proportions > 0]))
shannon_manual

singleton <- sum(amostras[, random_sample] == 1) # Number of singletons
singleton

doubleton <- sum(amostras[, random_sample] == 2) # Number of doubletons
doubleton

chao1_manual <- observed_richness_manual + (singleton^2 / (2 * max(doubleton, 1)))
chao1_manual


# Estimate richness metrics available in IN 'PHYLOSEQ' package using our ps package
alpha_diversity <- estimate_richness(ps)
head(alpha_diversity)[1:4,1:9]
View(alpha_diversity)


#Quickly ensuring readability
rownames(alpha_diversity) <- gsub("^X", "", rownames(alpha_diversity))   #Skip if samples doesn't begin with "X"

#Let's see how our manual sample results compare to our manually calculated results. Is it equal?
alpha_diversity$amostras <- rownames(alpha_diversity)
head(alpha_diversity)

subset(alpha_diversity, amostras == random_sample)

#Convert our sample data to dataframe
metadata <- sample_data(ps)

alpha_div_meta <- cbind(metadata, alpha_diversity)

head(alpha_div_meta)
View(alpha_div_meta)

# Plot Shannon diversity across urbanism groups
#p_Shannon <- ggplot(alpha_div_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  #geom_boxplot() 

#p_Shannon


#The beginning of Wednesday class

# Now let's visualize the alpha diversity differences between the different classes in our project.
# Choose a Metric and a CATEGORICAL Variable to generate the image. Then, test other metrics and other variables.
# Example: ggplot(alpha_div_meta, aes(x = urbanism, y = Shannon, fill = urbanism))

p_Metrica <- ggplot(alpha_div_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot() +
  ggtitle("Diversity according to the Shannon Metric, comparing SampleTypes between groups.")


p_Metrica
ggsave(filename = file.path(results_folder, "p_Metrica.png"), plot = p_Metrica)


# Now test with a numeric variable:
p_Metrica_2 <- ggplot(alpha_div_meta, aes(x = SampleType, y = Chao1, fill = SampleType)) +
  geom_boxplot() +
  ggtitle("Diversity according to the Chao1 Metric, comparing SampleTypes between groups.")

p_Metrica_2
ggsave(filename = file.path(results_folder, "p_Metrica_2.png"), plot = p_Metrica_2)

p_Metrica/p_Metrica_2


#You can also test more than one variable at a time:
p_Shannon_metrica2 <- ggplot(alpha_div_meta, aes(x = SampleType, y = Shannon, fill=SampleType2)) +
  geom_boxplot() +
  facet_wrap(~ SampleType2, scales = "free_x") + # Here we establish the facets.
  theme_minimal() +
  ggtitle("Diversity according to the Shannon Matrix, Comparing SampleType and Sampleype2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_Shannon_metrica2



# With the "spicy" package, calculate alpha diversity using the Faith method, which considers phylogenetic distances. Did you observe any difference?
install.packages("picante")
library(picante)
otu_matrix <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
  otu_matrix <- t(otu_matrix)
}

phy_tree <- phy_tree(ps)

# Calculate Faith's Phylogenetic Diversity (PD)
faith_pd <- pd(otu_matrix, phy_tree)

# view the result
faith_pd

# # Add Faith's metrics to metasmetadas:
alpha_div_meta2 <- cbind(alpha_div_meta, faith_pd)
head(alpha_div_meta2)


# Visualize:
p_PD_variavel <- ggplot(alpha_div_meta2, aes(x = SampleType, y = PD, fill=SampleType)) +
  geom_boxplot() 
p_PD_variavel



### Rarefaction
# Will our results be different if we standardize the number of reads per sample?
# Visualize a rarefaction curve for your data:
rarecurve(t(amostras), step = 200, cex = 0.1)

# Need to zoom in to establish an appropriate cutoff number? Replace 20000 with the maximum number of reads you want to view.
rarecurve(t(amostras), step = 200, cex = 0.6, xlim=c(0,20000))

# Rarefy your data. Replace "chosen depth" with the number of reads you want to rarefy:
ps.rare = rarefy_even_depth(ps, rngseed=1, sample.size=5000, replace=F)
# Using 5000 samples, the loss is only 7

# Calculate alpha diversity for the rarefied samples.
alpha_diversity_rare <- estimate_richness(ps.rare)

rownames(alpha_diversity_rare) <- gsub("^X", "", rownames(alpha_diversity_rare))
alpha_diversity_rare$sample <- rownames(alpha_diversity_rare)

head(alpha_diversity_rare)
View(alpha_diversity_rare)

# Combine the new metrics (with rarefaction and their metadata)
#alpha_div_rare_meta <- cbind(metadata, alpha_diversity_rare)


metadata <- sample_data(ps.rare) %>% data.frame()
alpha_div_rare_meta <- metadata %>% rownames_to_column("sample") %>%
  left_join(alpha_diversity_rare, by = "sample")

View(metadata)
View(alpha_div_rare_meta)

# Redo your plots and observe if there are differences in the results (without statistics yet).
p_Metrica_variavel_rare <- ggplot(alpha_div_rare_meta, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot() + 
  theme_minimal() +
  ggtitle("With Shannon rarefaction reads per sample") +
  theme(legend.position = "none")

p_Metrica_variavel_rare
p_Metrica_variavel_rare/p_Metrica


#tests
# Test normality of Shannon diversity using Shapiro-Wilk test:
shapiro.test(alpha_div_rare_meta$Shannon)

# If normality, For 2 groups, use a t-test.
t.test(Shannon ~ SampleType, data = alpha_div_meta)    #these variables wont perform cos there are more than 2 groups in the variables.

# For multiple groups, a One-way ANOVA:
anova_result <- aov(Shannon ~ SampleType, data = alpha_div_meta)
summary(anova_result)

# Non parametric: For 2 groups, use a Wilcoxon rank-sum test (also known as the Mann-Whitney U test)
wilcox.test(Shannon ~ SampleType, data = alpha_div_meta)     #these variables wont perform cos there are more than 2 groups in the variables.

# For multiple groups, use Kruskal-Wallis test
kruskal.test(Shannon ~ SampleType, data = alpha_div_meta)

# For ANOVA, use Tukey’s Honest Significant Differences (Tukey HSD) for post-hoc analysis.

# Post-hoc test for ANOVA
TukeyHSD(anova_result)

# For Kruskal-Wallis, use pairwise comparisons with Wilcoxon tests (adjusting for multiple comparisons).

# Post-hoc pairwise comparisons using Wilcoxon test
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$country, p.adjust.method = "BH")

par(mar = c(5, 13, 4, 2)) 
plot(TukeyHSD(anova_result), las = 1)
## Para adicionar ao plot:

# + stat_compare_means(method="wilcox.test") ou  + stat_compare_means(method="t.test")





###############################################
# -------      Beta diversidade       ------- #
###############################################












# Perform CAP scaling (dbRDA) on Bray-Curtis distance with multiple variables
cap_bray_multi <- capscale(bray_dist ~ Variavel1 + Variavel2 + Variavel3 + Variavel4, data = sample_data_df)

summary(cap_bray_multi)

plot(cap_bray_multi, display = c("Variavel1", "Variavel2"), scaling = 1)


scores_dbRDA <- as.data.frame(scores(cap_bray_multi, display = "Variavel1"))
scores_species <- as.data.frame(scores(cap_bray_multi, display = "Variavel2"))

colnames(scores_dbRDA) <- c("CAP1", "CAP2")
colnames(scores_species) <- c("CAP1", "CAP2")

ggplot() +
  geom_point(aes(x = scores_dbRDA$CAP1, y = scores_dbRDA$CAP2, color = sample_data_df$Region)) +
  geom_segment(aes(x = 0, y = 0, xend = scores_species$CAP1, yend = scores_species$CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(aes(x = scores_species$CAP1, y = scores_species$CAP2, label = rownames(scores_species)), 
            vjust = 1, hjust = 1, color = "black") +
  theme_minimal() +
  labs(title = "dbRDA (CAP) on Bray-Curtis Distance", x = "CAP1", y = "CAP2")



###############################################
# --- Análise estatística: PERMANOVA --- #
###############################################


# Additional PERMANOVA Post-doc test
# Pairwise adonis function:
adonis_result <- adonis2(formula = bray_dist ~ sample_type, data = as(sample_data(ps.rel.g), "data.frame"))

metadata_adonis <- as(sample_data(ps.rel.g), "data.frame")
pairwise_results <- pairwise.adonis2(bray_dist ~ Location, data = metadata_adonis, method = "bray")






###################################
# --- Differential abundances --- #
###################################




# Se não funcionar o DESeq2
physeq_transformed <- transform_sample_counts(physeq_filtered_nonzero, function(x) x + 1)


# DESeq2
dds <- phyloseq_to_deseq2(physeq_filtered_nonzero, ~ Region)
dds <- DESeq(dds) 

# Resultados
res <- results(dds)
sig_deseq <- subset(res, padj < 0.05)

sig_deseq_taxa <- row.names(sig_deseq)

# Maaslin2 Analysis

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata, 
                    output = "../maaslin2_output", 
                    fixed_effects = c("Region"),
                    random_effects = c(), # Include random effects that you want to adjust for, e.g. sample IDs when you have longitudinal study
                    transform = "NONE", # Possible transformation choices:  LOG, LOGIT, AST, NONE 
                    normalization= 'TSS', # Choices: TSS, CLR, CSS, NONE, TMM
                    reference = 'Region,Norte',
                    standardize = FALSE,
                    plot_heatmap = F, 
                    plot_scatter = T)


# ANCOM-BC Analysis

ancom_result <- ancombc(phyloseq = physeq_filtered_nonzero, formula = "Region")
head(ancom_result$res$q_val)


# Abundâncias relativas em log
physeq_rel <- transform_sample_counts(physeq_filtered_nonzero, function(x) x / sum(x))
physeq_log <- transform_sample_counts(physeq_rel, function(x) log(x + 0.001))

