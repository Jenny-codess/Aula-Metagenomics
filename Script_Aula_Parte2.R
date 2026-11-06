### Disciplina PPGBM/UFPA -- Bioinformática voltada para análises metagenômicas
# Ana Paula Schaan, Setembro 2026

###############################################
# ----------      Composição       ---------- #
###############################################

# Você pode aglomerar os ASVs ao nível desejado, dependendo da resolução taxonômica desejada.
rank_names(ps)
ps_phylum <- suppressMessages(tax_glom(ps, taxrank = "Phylum"))
ps_family <- suppressMessages(tax_glom(ps, taxrank = "Family"))
ps_genus <- suppressMessages(tax_glom(ps, taxrank = "Genus"))
# Teste outros também, se quiser

# Gráfico de composição utilizando dados normalizados (abundância relativa)

# Gerar objeto phyloseq aglomerado em abundância relativa
ps_phylum.rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
# Faça o mesmo para os demais níveis taxonômicos que vc escolheu.

#Gráfico
grafico_filo <- plot_bar(ps_phylum.rel, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) + 
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") + # Substitua SUA_VARIAVEL pela sua variável 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + 
  theme(legend.position="bottom") +
  theme(legend.position="bottom", axis.text.x = element_text(size=2))

grafico_filo

ggsave("plot_composition_filo.pdf",grafico_filo)

# Foi muita informação?
# Vamos plotar somente os táxons mais abundantes! Quantos você deseja mostrar?
ps_phylum.rel.tops <- ps_phylum.rel
phyla_sums <- taxa_sums(ps_phylum.rel.tops)

# Selecionar os X mais abundantes (substitua o 12 pelo valor que vc deseja mostrar)
top_12_phyla <- names(sort(phyla_sums, decreasing = TRUE)[1:12]) #Selecione aqui quantos serão mostrados

# Substituir os que não estão entre os 12 com "Outros"
tax_table(ps_phylum.rel.tops)[!(taxa_names(ps_phylum.rel.tops) %in% top_12_phyla), "Phylum"] <- "Others"

# Plot
grafico_filo_tops <- plot_bar(ps_phylum.rel.tops, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +  # Remove border
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + 
  theme(legend.position="bottom", axis.text.x = element_text(size=2))
grafico_filo_tops

ggsave("plot_topX.pdf")

# Core microbiome -- Em vez de selecionar os top 12 arbitrariamente, selecione os táxons presentes em pelo menos 50% das amostras, por exemplo.
ps_phylum.rel.core = filter_taxa(ps_phylum.rel, function(x) sum(x > 0) > (0.5*length(x)), TRUE)

# Quantos taxons fazem parte do CORE MICROBIOME?
ps_phylum.rel.core

grafico_filo_core <- plot_bar(ps_phylum.rel.core, fill = "Phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +  # Remove border
  facet_grid( ~ SampleType, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Abundância relativa por filo") + theme(legend.position="bottom") +
  theme(legend.position="bottom", axis.text.x = element_text(size=2))

grafico_filo_core

ggsave("plot_core_filos.pdf")


###############################################
# -------      Beta diversidade       ------- #
###############################################

# Aglomerar a nível de gênero para análises de diversidade beta (ou vc pode manter a nível de ASV)
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

# Visualizar o formato dos dados de distâncias
View(as.matrix(bray_dist))
summary(bray_dist)

###############################################
# --- Métodos de ordenação --- #
###############################################

# PCoA com Bray-Curtis
ord_bray <- ordinate(ps.rel.g, method = "MDS", distance = bray_dist)
ord_bray
# Plotar PCoA
p_bray <- plot_ordination(ps.rel.g, ord_bray, color = "environmental_medium") + # Troque "urbanism" pela sua variável
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCoA com Bray-Curtis")

p_bray


#Sample group Ellipses
p_bray_ellipses <- plot_ordination(ps.rel.g, ord_bray, color = "SampleType") +
  geom_point(size = 3) +
  stat_ellipse(aes(group = SampleType)) +
  theme_minimal() +
  ggtitle("PCoA com Bray-Curtis e Elipses por Nature")

p_bray_ellipses


###############################################
# --- Análise estatística: PERMANOVA --- #
###############################################

# Usar adonis para testar diferenças significativas entre grupos (PERMANOVA)
adonis_result <- adonis2(bray_dist ~ SampleType, data = as(sample_data(ps.rel.g), "data.frame")) # Teste com outras métricas e outras variáveis. Qual melhor explica a variação nos seus dados?
adonis_result

# Será que seus resultados são enviesados pela dispersão das amostras?

# Testar a dispersão beta usando betadisper (dispersão de grupos)
betadisper_bray <- betadisper(bray_dist, as.factor(sample_data(ps.rel.g)$SampleType)) # Troque "urbanism" pela sua variável
anova(betadisper_bray)  # ANOVA para testar diferenças de dispersão

# Additional PERMANOVA Post-doc test -- Use se você tiver testando diferenças entre mais de 2 grupos:

#Pairwise adonis function:
adonis_result <- adonis2(
  formula = bray_dist ~ urbanism, 
  data = as(sample_data(ps.rel.g), "data.frame")
)
adonis_result

metadata_adonis <- as(sample_data(ps.rel.g), "data.frame")
pairwise_results <- pairwise.adonis2(bray_dist ~ urbanism, data = metadata_adonis, method = "bray")
pairwise_results

#Calcular dispersão
bray_disper<-vegan::betadisper(bray_dist, type= "centroid", phyloseq::sample_data(ps.rel.g)$urbanism)

sample_data(ps.rel.g)$bray_disper<-bray_disper$distances
betadistances=as.data.frame(sample_data(ps.rel.g))
betadistances

bray_ord_p = ggplot(betadistances, aes(x=reorder(urbanism, bray_disper, FUN=mean), y=bray_disper, fill=urbanism)) + 
  geom_boxplot() + theme_classic2() + stat_compare_means()
bray_ord_p

