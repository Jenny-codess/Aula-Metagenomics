### Disciplina PPGBM/UFPA -- Bioinformática voltada para análises metagenômicas
# Ana Paula Schaan, Setembro 2026


# ---- Verificar e carregar os pacotes ----

required_packages <- c(
  "tidyverse",
  "phyloseq", # estatisticas, processamento
  "vegan", # processamento
  "ggpubr", # plots
  "cowplot", # plots
  "patchwork", # plots
  "ape", # arvores filogenéticas
  "picante",
  "compositions",
  "pheatmap",
  "DESeq2", # abundância diferencial 
  "ANCOMBC", # abundância diferencial
  "Maaslin2", # abundância diferencial
  "scater",
  "pairwiseAdonis"
)

# Carregar os pacotes
lapply(
  required_packages,
  library,
  character.only = TRUE
)


# Estabelecer o diretório de trabalho -- substitua aqui pelo seu diretório OU clique em "Sessão/Session -> Set Working Directory -> Choose directory" e selecione a pasta ontem estão os arquivos referentes ao seu projeto
setwd("/Users/apschaan/Documents/Kiel/UKSH-CAU/Classes/KMC_workshop/") # sua pasta

##################################
# -----   Importar dados   ----- #
##################################

# Importar tabela de ASVs
tabASV <- read.table("feature-table-gmbc1.tsv", header=T, check.names=F, row.names=1) # substitua o nome do arquivo entre aspas para o nome do seu arquivo de ASVs
tabASV <- as.matrix(tabASV)
head(tabASV)[1:5,1:5]

# Importar tabela de taxonomia
tabTax <- read.csv("taxonomy_gmbc1_form.csv", header=T, row.names=1) # substitua o nome do arquivo entre aspas para o nome do seu arquivo de taxonomia
# se o seu arquivo de Tax for .txt, troque o "read.csv" por: read_table ou read.delim ou read_tsv
tabTax <- as.matrix(tabTax)
head(tabTax)[1:5,1:5]

# Importar árvore filogenética
Arv <- read_tree("tree-gmbc1.nwk")

# Importar metadados
metadataTab <- read.table("metadata_gmbc_bn10_complete.tsv", header=TRUE, row.names=1, sep="\t")
# se o seu arquivo de Metadata for .csv, troque o "read.table" por: read.csv ou read_csv ou troque o sep="\t" por sep=","

##################################
# --- Importar para phyloseq --- #
##################################

ASV <- otu_table(tabASV, taxa_are_rows = TRUE)
TAX <- tax_table(tabTax)
META <- sample_data(metadataTab)
TREE <- Arv

ps <- phyloseq(ASV, TAX, META, TREE)

# Check the phyloseq object
ps

# Verificando as classes de cada coluna (variável) -- observe quantas variáveis são caractéres, quantas são numéricas etc
variable_classes <- sapply(META, class)
tail(variable_classes)


# Explore o phyloseq object
nsamples(ps) # Número de amostras
ntaxa(ps) # Número de ASVs
rank_names(ps) # Níveis taxonômicos
sample_variables(ps) # Quais metadados estão disponíveis
sample_sums(ps) # Quantas reads por amostra vc tem?

# Distribuição de reads por amostra

# Crie uma tabela que resuma número de reads por ASV
tabelaContagem <- data.frame(ContagemTotal = taxa_sums(ps))
head(tabelaContagem)

# Use o `ggplot2` para visualizar a distribuição de reads num histograma. O que chama sua atenção?
ggplot(tabelaContagem, aes(ContagemTotal)) + geom_histogram() 

############################
# --- Alfa diversidade --- #
############################

# Vamos experimentar calcular os índices de diversidade para uma amostra manualmente.
amostras <- tabASV[, colnames(tabASV) %in% rownames(metadataTab)] 

# Pegue uma amostra aleatória
random_sample <- sample(colnames(amostras), size=1)
random_sample

# Calcule o número total de reads nessa amostra aleatória (isso será o denominador do cálculo de proporções)
total_reads <- sum(amostras[, random_sample])

# Calcule a proporção (abundância) de reads atribuidas a cada ASV nessa amostra aleatória
proportions <- amostras[, random_sample] / total_reads

# Calcule manual o número de espécies observadas (basta somar)
observed_richness_manual <- sum(amostras[, random_sample] > 0)
observed_richness_manual

# Calcule manualmente o índice de Shannon. Observe que levamos as proporções (abundâncias em consideração)
shannon_manual <- -sum(proportions[proportions > 0] * log(proportions[proportions > 0]))
shannon_manual

# Para calcular Chao1, precisamos calcular quantas singletons e doubletons temos na nossa amostra
singleton <- sum(amostras[, random_sample] == 1) # Number of singletons
singleton

doubleton <- sum(amostras[, random_sample] == 2) # Number of doubletons
doubleton

# Com essas informações, vamos calcular o índice de Shannon. Observe se ele é diferente do número de espécies observadas. Se for igual, qual o motivo?
chao1_manual <- observed_richness_manual + (singleton^2 / (2 * max(doubleton, 1)))
chao1_manual

# O bom é que não precisamos fazer esse cálculo manualmente para todas as amostras. O pacote phyloseq tem uma função chamada `estimate_richness` que já calcula automaticamente vários índices.
alpha_diversity <- estimate_richness(ps)
# Veja os resultados:
head(alpha_diversity)[1:4,1:9]
View(alpha_diversity)

# Rode o seguinte comando para arrumar o nome das amostras 
rownames(alpha_diversity) <- gsub("^X", "", rownames(alpha_diversity))

# Se vc fez o cálculo manual, compare os resultados pela função e pelos seus cálculos. Deu igual?
alpha_diversity$amostra <- rownames(alpha_diversity)
subset(alpha_diversity, amostra == random_sample)

# Vamos unir os resultados de diversidade alfa aos metadados:
metadata <- sample_data(ps)
alpha_div_meta <- cbind(metadata, alpha_diversity)

# Veja a combinação das métricas e dos seus metadados:
head(alpha_div_meta)

#### ------ FIM DA AULA DE TERÇA ------ ####

### Início da aula de quarta:
# Agora vamos visualizar as diferenças de diversidade alfa entre as diferentes classes do nosso projeto. 

# Escolha uma Metrica e uma Variavel CATEGORICA para gerar a imagem. Em seguida, teste outras métricas e outras variáveis.
# Exemplo: ggplot(alpha_div_meta, aes(x = urbanism, y = Shannon, fill = urbanism))
p_Metrica <- ggplot(alpha_div_meta, aes(x = Variavel, y = Metrica, fill = Variavel)) +
  geom_boxplot() +
  ggtitle("Diversidade pela Metrica XXX, comparando Variavel tal entre grupos")
  

p_Metrica

# Agora teste com uma variável numérica:
p_Metrica_2 <- ggplot(alpha_div_meta, aes(x = Variavel, y = Metrica, fill = Variavel)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Diversidade pela Metrica XXX, comparando Variavel tal entre grupos") +
  

p_Metrica_2


# Salve o seu plot (renomeie o arquivo com base na metrica e na variável que você utilizou):
ggsave(p_Metrica, "plot_metrica_variavel.pdf")


# Você pode também testar mais de uma variável de uma vez:
p_Shannon_metrica2 <- ggplot(alpha_div_meta, aes(x = Variavel1, y = Metrica, fill=Variavel2)) +
  geom_boxplot() +
  facet_wrap(~ Variavel2, scales = "free_x") + # Aqui estabelecemos as facetas
  theme_minimal() +
  ggtitle("Diversidade pela Metrica XXX, comparando Variavel 1 e Variavel 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_Shannon_metrica2





# Com o pacote "picante", calcule a alfa diversidade usando o método Faith, que considera distâncias filogenéticas. Você observou alguma diferenca?
install.packages("picante")
library(picante)
otu_matrix <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
  otu_matrix <- t(otu_matrix)
}

phy_tree <- phy_tree(ps)

# Calcule Faith's Phylogenetic Diversity (PD)
faith_pd <- pd(otu_matrix, phy_tree)

# Veja o resultado
faith_pd

# Una a métrica de Faith aos metadados:
alpha_div_meta2 <- cbind(alpha_div_meta, faith_pd)
head(alpha_div_meta2)

# Visualize:
p_PD_variavel <- ggplot(alpha_div_meta2, aes(x = Variavel, y = PD, fill=Variavel)) +
  geom_boxplot() 
p_PD_variavel

### Rarefação
# Será que nossos resultados serão diferentes se uniformizarmos o número de reads por amostra?
# Visualize uma curva de rarefação para seus dados:
rarecurve(t(amostras), step = 200, cex = 0.1)

# Precisa de um zoom para estabelecer um número de corte adequado? Substitua 20000 para o número máximo de reads que deseja visualizar.
rarecurve(t(subset_samples), step = 200, cex = 0.6, xlim=c(0,20000))

# Faça a rarefação dos seus dados. Substitua "chosen_depth" pelo número de reads que deseja rarefazer:
ps.rare = rarefy_even_depth(ps, rngseed=1, sample.size=chosen_depth, replace=F)

# Calcule diversidade alfa para as amostras com rarefação.
alpha_diversity_rare <- estimate_richness(ps.rare)

rownames(alpha_diversity_rare) <- gsub("^X", "", rownames(alpha_diversity_rare))
alpha_diversity_rare$sample <- rownames(alpha_diversity_rare)

head(alpha_diversity_rare)

# Combine as novas métricas (com rarefação com seus metadados)
alpha_div_rare_meta <- cbind(metadata, alpha_diversity_rare)

# Refaça seus plots e observe se há diferenças nos resultados (sem estatísticas ainda)
p_Metrica_variavel_rare <- ggplot(alpha_div_rare_meta, aes(x = Variavel, y = Metrica, fill = Variavel)) +
  geom_boxplot() + 
  theme_minimal() +
  ggtitle("Com rarefação de XXX reads por amostra") +
  theme(legend.position = "none")




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

