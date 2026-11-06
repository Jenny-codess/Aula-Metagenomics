### PPGBM/UFPA Course -- Bioinformatics focused on metagenomic analyses
# Ana Paula Schaan, September 2026

###############################################
# ---------- Composition ---------- #
###############################################

# You can cluster the ASVs to the desired level, depending on the desired taxonomic resolution.
rank_names(ps)
ps_phylum <- suppressMessages(tax_glom(ps, taxrank = "Phylum"))
ps_family <- suppressMessages(tax_glom(ps, taxrank = "Family"))
ps_genus <- suppressMessages(tax_glom(ps, taxrank = "Genus"))
# Try others as well, if you want

# Composition graph using normalized data (relative abundance)

# Generate phyloseq cluster object in relative abundance
ps_phylum.rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
# Do the same for the other taxonomic levels you chose.

#Chart
grafico_filo <- plot_bar(ps_phylum.rel, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) + 
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") +   #Replace YOUR_VARIABLE with your variable. 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + 
  theme(legend.position="bottom") +
  theme(legend.position="bottom", axis.text.x = element_text(size=2))

grafico_filo

#ggsave("plot_composition_filo.pdf",grafico_filo)
ggsave(filename = file.path(results_folder, "grafico_filo.png"), plot = grafico_filo)


# Was that too much information?
# Let's plot only the most abundant taxa! How many do you want to show?
ps_phylum.rel.tops <- ps_phylum.rel
phyla_sums <- taxa_sums(ps_phylum.rel.tops)

# Select the most abundant X's (replace 12 with the value you want to display)
top_12_phyla <- names(sort(phyla_sums, decreasing = TRUE)[1:12]) #Select here how many will be shown.

# Replace those not among the 12 with "Others"
tax_table(ps_phylum.rel.tops)[!(taxa_names(ps_phylum.rel.tops) %in% top_12_phyla), "Phylum"] <- "Others"

# Plot
grafico_filo_tops <- plot_bar(ps_phylum.rel.tops, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +  # Remove border
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + 
  theme(legend.position="bottom", axis.text.x = element_text(size=2))
grafico_filo_tops

#ggsave("plot_topX.pdf")
ggsave(filename = file.path(results_folder, "grafico_filo_tops.png"), plot = grafico_filo_tops)

# Core microbiome -- Instead of arbitrarily selecting the top 12, select taxa present in at least 50% of the samples, for example.
ps_phylum.rel.core = filter_taxa(ps_phylum.rel, function(x) sum(x > 0) > (0.5*length(x)), TRUE)

# How many taxa are part of the core microbiome?
ps_phylum.rel.core

grafico_filo_core <- plot_bar(ps_phylum.rel.core, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +  # Remove border
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + theme(legend.position="bottom") +
  theme(legend.position="bottom", axis.text.x = element_text(size=2))

grafico_filo_core

#ggsave("plot_core_filos.pdf")
ggsave(filename = file.path(results_folder, "grafico_filo_core.png"), plot = grafico_filo_core)

###############################################
# -------      Beta diversidade       ------- #
###############################################

# Cluster at the gender level for beta diversity analyses (or you can keep it at the ASV level).
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps.rel.g <- suppressMessages(tax_glom(ps.rel, "Genus"))

# Bray-Curtis
bray_dist <- phyloseq::distance(ps.rel.g, method = "bray")

# Jaccard
jaccard_dist <- phyloseq::distance(ps.rel.g, method = "jaccard")

# Unweighted UniFrac
unifrac_unw <- phyloseq::distance(ps.rel.g, method = "unifrac")

# Weighted UniFrac
unifrac_w <- phyloseq::distance(ps.rel.g, method = "wunifrac")

# View the format of the distance data.
View(as.matrix(bray_dist))
summary(bray_dist)

###############################################
# --- Sorting methods --- #
###############################################

# PCoA com Bray-Curtis
ord_bray <- ordinate(ps.rel.g, method = "MDS", distance = bray_dist)
ord_bray
# PCoA Plot
p_bray <- plot_ordination(ps.rel.g, ord_bray, color = "environmental_medium") + # Troque "environmental_medium" pela sua variável
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCoA com Bray-Curtis")

p_bray


#Sample group Ellipses
p_bray_ellipses <- plot_ordination(ps.rel.g, ord_bray, color = "environmental_medium") +
  geom_point(size = 3) +
  stat_ellipse(aes(group = environmental_medium)) +
  theme_minimal() +
  ggtitle("PCoA com Bray-Curtis e Elipses por Nature")

p_bray_ellipses


###############################################
# --- Statistical analysis: PERMANOVA --- #
###############################################

# Using Adonis to test for significant differences between groups (PERMANOVA)
adonis_result <- adonis2(bray_dist ~ SampleType, data = as(sample_data(ps.rel.g), "data.frame")) # Teste com outras métricas e outras variáveis. Qual melhor explica a variação nos seus dados?
adonis_result

# Are your results biased by sample dispersion?

# Test for beta dispersion using betadisper (group dispersion)
betadisper_bray <- betadisper(bray_dist, as.factor(sample_data(ps.rel.g)$SampleType)) # Troque "SampleType" pela sua variável
anova(betadisper_bray)  # ANOVA to test differences in dispersion

# Additional PERMANOVA Post-doc test -- Use if you are testing differences between more than 2 groups:

#Pairwise adonis function:
adonis_result <- adonis2(
  formula = bray_dist ~ SampleType, 
  data = as(sample_data(ps.rel.g), "data.frame")
)
adonis_result

# Install devtools or remotes if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

metadata_adonis <- as(sample_data(ps.rel.g), "data.frame")
pairwise_results <- pairwise.adonis2(bray_dist ~ SampleType, data = metadata_adonis, method = "bray")
pairwise_results

#Calculate dispersion.
bray_disper<-vegan::betadisper(bray_dist, type= "centroid", phyloseq::sample_data(ps.rel.g)$SampleType)

sample_data(ps.rel.g)$bray_disper<-bray_disper$distances
betadistances=as.data.frame(sample_data(ps.rel.g))
betadistances

bray_ord_p = ggplot(betadistances, aes(x=reorder(SampleType, bray_disper, FUN=mean), y=bray_disper, fill=SampleType)) + 
  geom_boxplot() + theme_classic2() + stat_compare_means()
bray_ord_p

