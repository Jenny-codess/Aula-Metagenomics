### Disciplina PPGBM/UFPA -- Bioinformática voltada para análises metagenômicas
# Ana Paula Schaan, Setembro 2026

###################################
# --- Abundância diferencial --- #
###################################

library(phyloseq)
library(DESeq2)
library(Maaslin2)
library(ANCOMBC)
library(tidyverse)

# Este script compara três métodos:
# DESeq2, MaAsLin2 e ANCOM-BC.
# variáveis categóricas com dois ou mais grupos.

sample_data(ps)

# escolha uma variável categórica presente nos seus metadados
variavel_grupo <- "urbanism"

sample_data(ps) %>% data.frame() %>% select(variavel_grupo) %>% distinct()

# Escolha um nível a ser usado como referência
grupo_referencia <- "rural"

# Nível taxonômico que será testado (pode ser "Genus", "Family", "Species" etc.)
nivel_escolhido <- "Family"


# Manter taxa presentes em pelo menos X% das amostras
cut_off <- 0.1 # Troque de acordo com sua preferência -- O que faz sentido para o seu estudo? e para o seu artigo?

# Preparo dos dados #

metadata_da <- data.frame(sample_data(ps))

variavel_grupo %in% colnames(metadata_da) # tem que retornar TRUE

# Remove samples without data for chosen variable
amostras_validas <- rownames(metadata_da)[
  !is.na(metadata_da[[variavel_grupo]]) &
    metadata_da[[variavel_grupo]] != ""
]

ps_da <- prune_samples(amostras_validas, ps)

# Níveis encontrados nos dados
niveis_grupo <- unique(
  as.character(sample_data(ps_da)[[variavel_grupo]])
)

niveis_grupo

grupos_comparacao <- setdiff(
  niveis_grupo,
  grupo_referencia
)

# Técnico: equilibrando os mapas internos
mapa_grupos <- tibble(
  grupo_original = c(
    grupo_referencia,
    grupos_comparacao
  ),
  codigo_grupo = c(
    "REF",
    paste0("G", seq_along(grupos_comparacao))
  )
) %>%
  mutate(
    comparison = if_else(
      codigo_grupo == "REF",
      NA_character_,
      grupo_original
    ),
    reference = grupo_referencia,
    contrast = if_else(
      codigo_grupo == "REF",
      NA_character_,
      paste(grupo_original, "vs", grupo_referencia)
    )
  )

mapa_grupos

grupo_observado <- as.character(
  sample_data(ps_da)[[variavel_grupo]]
)

sample_data(ps_da)$grupo_modelo <- factor(
  mapa_grupos$codigo_grupo[
    match(grupo_observado, mapa_grupos$grupo_original)
  ],
  levels = mapa_grupos$codigo_grupo
)

table(
  grupo_original = grupo_observado,
  codigo_modelo = sample_data(ps_da)$grupo_modelo
)

# Tabela de taxonomia

ps_tax <- suppressMessages(
  tax_glom(
    ps_da,
    taxrank = nivel_escolhido
  )
)

# Remover táxons sem classificação no nível escolhido
nomes_taxonomicos <- as.character(
  tax_table(ps_tax)[, nivel_escolhido]
)

taxa_classificados <-
  !is.na(nomes_taxonomicos) &
  nomes_taxonomicos != "" &
  !grepl(
    "^[a-z]__$",
    nomes_taxonomicos,
    ignore.case = TRUE
  )

ps_tax <- prune_taxa(taxa_classificados, ps_tax)

# Filtro de prevalência (aqui vamos aplicar seu cut-off)

ps_tax <- filter_taxa(
  ps_tax,
  function(x) {
    sum(x > 0) >= cut_off * length(x)
  },
  prune = TRUE
)

ps_tax <- prune_taxa(
  taxa_sums(ps_tax) > 0,
  ps_tax
)

ps_tax

########################
# Tabela de taxonomia   #
########################

taxonomia_base <- as.data.frame(
  tax_table(ps_tax),
  check.names = FALSE
) %>%
  rownames_to_column("taxon_id") %>%
  mutate(
    taxon_name = as.character(.data[[nivel_escolhido]]),
    taxon_name = str_remove(taxon_name, "^[a-z]__"),
    taxon_name = str_replace_all(taxon_name, "_", " "),
    taxonomic_rank = nivel_escolhido
  )

adicionar_taxonomia <- function(resultado, coluna_id) {
  taxon_metodo <- as.character(resultado[[coluna_id]])
  
  indice_taxon <- match(
    taxon_metodo,
    taxonomia_base$taxon_id
  )
  
  precisa_make_names <- is.na(indice_taxon)
  
  indice_taxon[precisa_make_names] <- match(
    taxon_metodo[precisa_make_names],
    make.names(taxonomia_base$taxon_id)
  )
  
  resultado$taxon_id <- taxonomia_base$taxon_id[indice_taxon]
  resultado$taxon_name <- taxonomia_base$taxon_name[indice_taxon]
  resultado$taxonomic_rank <- taxonomia_base$taxonomic_rank[indice_taxon]
  
  if (any(is.na(resultado$taxon_id))) {
    stop(
      paste(
        "Não foi possível recuperar a taxonomia de",
        sum(is.na(resultado$taxon_id)),
        "identificador(es)."
      )
    )
  }
  
  resultado %>%
    select(-all_of(coluna_id))
}

# comparações possíveis

mapa_contrastes <- mapa_grupos %>%
  filter(codigo_grupo != "REF") %>%
  select(
    codigo_grupo,
    comparison,
    reference,
    contrast
  )

mapa_contrastes

################
#    DESeq2   #
################

formula_da <- ~ grupo_modelo

dds <- phyloseq_to_deseq2(
  ps_tax,
  formula_da
)

dds <- DESeq(
  dds,
  sfType = "poscounts",
  quiet = TRUE
)

res_deseq <- map_dfr(
  mapa_contrastes$codigo_grupo,
  function(codigo_atual) {
    resultado <- results(
      dds,
      contrast = c(
        "grupo_modelo",
        codigo_atual,
        "REF"
      ),
      alpha = 0.05
    )
    
    as.data.frame(resultado) %>%
      rownames_to_column("taxon_key") %>%
      mutate(codigo_grupo = codigo_atual)
  }
) %>%
  left_join(
    mapa_contrastes,
    by = "codigo_grupo"
  ) %>%
  adicionar_taxonomia("taxon_key") %>%
  mutate(method = "DESeq2")

# Resultados
res_deseq

sig_deseq <- res_deseq %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

sig_deseq %>%
  dplyr::count(contrast, name = "numero_taxa_significativos")

################
# MaAsLin2      #
################

maaslin_data <- as(
  otu_table(ps_tax),
  "matrix"
)

if (taxa_are_rows(ps_tax)) {
  maaslin_data <- t(maaslin_data)
}

maaslin_metadata <- data.frame(
  sample_data(ps_tax)
)

maaslin_metadata <- maaslin_metadata[
  rownames(maaslin_data),
  ,
  drop = FALSE
]

maaslin_output <- file.path(
  getwd(),
  paste0(
    "maaslin2_",
    make.names(variavel_grupo),
    "_output"
  )
)

dir.create(
  maaslin_output,
  recursive = TRUE,
  showWarnings = FALSE
)

logging::logReset()
fit_maaslin <- Maaslin2(
  input_data = maaslin_data,
  input_metadata = maaslin_metadata,
  output = maaslin_output,
  fixed_effects = "grupo_modelo",
  normalization = "TSS",
  transform = "LOG",
  reference = "grupo_modelo,REF",
  standardize = FALSE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

# Resultados
res_maaslin <- fit_maaslin$results %>%
  filter(
    metadata == "grupo_modelo",
    value %in% mapa_contrastes$codigo_grupo
  ) %>%
  mutate(
    codigo_grupo = as.character(value),
    taxon_key = as.character(feature)
  ) %>%
  left_join(
    mapa_contrastes,
    by = "codigo_grupo"
  ) %>%
  adicionar_taxonomia("taxon_key") %>%
  mutate(method = "MaAsLin2")

res_maaslin

sig_maaslin <- res_maaslin %>%
  filter(
    !is.na(qval),
    qval < 0.05
  )

sig_maaslin %>%
  dplyr::count(contrast, name = "numero_taxa_significativos")

################
# ANCOM-BC      #
################

ancom_result <- ancombc(
  phyloseq = ps_tax,
  formula = "grupo_modelo",
  p_adj_method = "BH"
)

padronizar_tabela_ancom <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  
  if ("taxon" %in% colnames(x)) {
    x <- x %>%
      dplyr::rename(taxon_key = taxon)
  } else if ("taxon_id" %in% colnames(x)) {
    x <- x %>%
      dplyr::rename(taxon_key = taxon_id)
  } else {
    x <- x %>%
      rownames_to_column("taxon_key")
  }
  
  x
}

ancom_q <- padronizar_tabela_ancom(
  ancom_result$res$q_val
)

ancom_lfc <- padronizar_tabela_ancom(
  ancom_result$res$lfc
)

# Get coefficients from ANCOM and relate them to the original
mapa_coeficientes <- mapa_contrastes %>%
  mutate(
    coefficient = paste0(
      "grupo_modelo",
      codigo_grupo
    )
  ) %>%
  select(
    coefficient,
    codigo_grupo,
    comparison,
    reference,
    contrast
  )

# Can you check that your coefficients are all there?
coeficientes_ausentes <- setdiff(
  mapa_coeficientes$coefficient,
  intersect(
    colnames(ancom_q),
    colnames(ancom_lfc)
  )
)

length(coeficientes_ausentes) > 0 # should be FALSE!

ancom_q_long <- ancom_q %>%
  select(
    taxon_key,
    all_of(mapa_coeficientes$coefficient)
  ) %>%
  pivot_longer(
    cols = -taxon_key,
    names_to = "coefficient",
    values_to = "qval"
  )

ancom_lfc_long <- ancom_lfc %>%
  select(
    taxon_key,
    all_of(mapa_coeficientes$coefficient)
  ) %>%
  pivot_longer(
    cols = -taxon_key,
    names_to = "coefficient",
    values_to = "lfc"
  )

res_ancom <- ancom_q_long %>%
  left_join(
    ancom_lfc_long,
    by = c("taxon_key", "coefficient")
  ) %>%
  left_join(
    mapa_coeficientes,
    by = "coefficient"
  ) %>%
  adicionar_taxonomia("taxon_key") %>%
  mutate(method = "ANCOM-BC")

sig_ancom <- res_ancom %>%
  filter(
    !is.na(qval),
    qval < 0.05
  )

sig_ancom %>%
  dplyr::count(contrast, name = "numero_taxa_significativos")

################################
# Comparação entre os métodos  #
################################


checagem_contrastes <- dplyr::bind_rows(
  res_deseq %>%
    dplyr::distinct(
      comparison,
      reference,
      contrast
    ) %>%
    dplyr::mutate(
      method = "DESeq2"
    ),
  res_maaslin %>%
    dplyr::distinct(
      comparison,
      reference,
      contrast
    ) %>%
    dplyr::mutate(
      method = "MaAsLin2"
    ),
  res_ancom %>%
    dplyr::distinct(
      comparison,
      reference,
      contrast
    ) %>%
    dplyr::mutate(
      method = "ANCOM-BC"
    )
) %>%
  dplyr::select(
    method,
    comparison,
    reference,
    contrast
  ) %>%
  dplyr::arrange(
    comparison,
    method
  )

checagem_contrastes

garantir_taxon_id <- function(resultado) {
  if (
    "taxon_id" %in% colnames(resultado) &&
    "taxon_name" %in% colnames(resultado)
  ) {
    return(resultado)
  }
  colunas_possiveis <- c(
    "taxon_key",
    "taxon",
    "feature"
  )
  
  coluna_encontrada <- colunas_possiveis[
    colunas_possiveis %in% colnames(resultado)
  ][1]
  
  if (is.na(coluna_encontrada)) {
    stop(
      "Tax ID not found."
    )
  }
  taxon_metodo <- as.character(
    resultado[[coluna_encontrada]]
  )
  indice_taxon <- match(
    taxon_metodo,
    taxonomia_base$taxon_id
  )
  precisa_correcao <- is.na(indice_taxon)
  
  indice_taxon[precisa_correcao] <- match(
    taxon_metodo[precisa_correcao],
    make.names(taxonomia_base$taxon_id)
  )
  resultado$taxon_id <-
    taxonomia_base$taxon_id[indice_taxon]
  resultado$taxon_name <-
    taxonomia_base$taxon_name[indice_taxon]
  resultado
}

sig_deseq <- garantir_taxon_id(sig_deseq)
sig_maaslin <- garantir_taxon_id(sig_maaslin)
sig_ancom <- garantir_taxon_id(sig_ancom)

common_results <- sig_deseq %>%
  select(
    taxon_id,
    taxon_name,
    comparison,
    reference,
    contrast,
    log2FoldChange,
    padj
  ) %>%
  dplyr::rename(
    deseq_effect = log2FoldChange,
    deseq_qvalue = padj
  ) %>%
  inner_join(
    sig_maaslin %>%
      select(
        taxon_id,
        comparison,
        coef,
        qval
      ) %>%
      dplyr::rename(
        maaslin_effect = coef,
        maaslin_qvalue = qval
      ),
    by = c(
      "taxon_id",
      "comparison"
    )
  ) %>%
  inner_join(
    sig_ancom %>%
      select(
        taxon_id,
        comparison,
        lfc,
        qval
      ) %>%
      dplyr::rename(
        ancom_effect = lfc,
        ancom_qvalue = qval
      ),
    by = c(
      "taxon_id",
      "comparison"
    )
  ) %>%
  mutate(
    same_direction =
      sign(deseq_effect) == sign(maaslin_effect) &
      sign(deseq_effect) == sign(ancom_effect),
    largest_qvalue = pmax(
      deseq_qvalue,
      maaslin_qvalue,
      ancom_qvalue,
      na.rm = TRUE
    )
  )


### Consenso entre os 3 métodos! Teve algum?
consensus_results <- common_results %>%
  dplyr::filter(same_direction)

common_results %>%
  dplyr::count(
    contrast,
    same_direction,
    name = "numero_taxa"
  )


# Mostrar apenas contrastes realmente disponíveis nos resultados significativos.
contrastes_disponiveis <- bind_rows(
  sig_deseq %>% select(contrast),
  sig_maaslin %>% select(contrast),
  sig_ancom %>% select(contrast)
) %>%
  filter(!is.na(contrast)) %>%
  distinct(contrast) %>%
  arrange(contrast) %>%
  pull(contrast)

contrastes_disponiveis

# Escolha o número de uma das comparações acima
contraste_plot <- contrastes_disponiveis[1]

# Opções: "consenso", "DESeq2", "MaAsLin2" ou "ANCOM-BC"
metodo_plot <- "ANCOM-BC"

# Número máximo de táxons no gráfico, para o caso de teres centenas!
numero_maximo_taxa <- 10

resultados_plot <- switch(
  metodo_plot,
  "consenso" = consensus_results %>%
    mutate(significance = largest_qvalue),
  "DESeq2" = sig_deseq %>%
    mutate(significance = padj),
  "MaAsLin2" = sig_maaslin %>%
    mutate(significance = qval),
  "ANCOM-BC" = sig_ancom %>%
    mutate(significance = qval)
)

resultados_plot <- resultados_plot %>%
  filter(contrast == contraste_plot) %>%
  arrange(significance) %>%
  distinct(taxon_id, .keep_all = TRUE) %>%
  slice_head(n = numero_maximo_taxa)


resultados_plot %>%
  select(
    taxon_id,
    taxon_name,
    contrast,
    significance
  )

# Gráfico de abundância 
taxa_plot <- resultados_plot$taxon_id

ps_rel <- transform_sample_counts(
  ps_tax,
  function(x) {
    x / sum(x)
  }
)

ps_plot <- prune_taxa(
  taxa_plot,
  ps_rel
)

plot_data <- psmelt(ps_plot) %>%
  left_join(
    taxonomia_base %>%
      select(
        taxon_id,
        taxon_name
      ),
    by = c("OTU" = "taxon_id")
  ) %>%
  mutate(
    abundance_percent = Abundance * 100
  )

# Manter a ordem de significância 
plot_data$taxon_name <- factor(
  plot_data$taxon_name,
  levels = resultados_plot$taxon_name
)

diff_abund_plot_1 <- ggplot(
  plot_data,
  aes(
    x = .data[[variavel_grupo]],
    y = abundance_percent,
    fill = .data[[variavel_grupo]]
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  facet_wrap(
    ~ taxon_name,
    scales = "free_y"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = paste(
      metodo_plot,
      "—",
      contraste_plot
    ),
    x = variavel_grupo,
    y = "Relative abundances (%)"
  )

diff_abund_plot_1

ggsave("diff_abund_genus_variable_contrast.pdf")
