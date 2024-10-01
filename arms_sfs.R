# MERGE fraction sample ---------------------
ps_sfs <- merge_samples(ps, group = "factor_r")
df_sfs <- data.frame(sample_data(ps_sfs))
df_sfs <- df_sfs  %>%
  mutate(study = rep("arms", nrow(.)),
         site = rep(c("Algarrobo", "Las_Cruces"), each = 8 ),
         ocean_depth= rep(c("30m", "60m"), each = 4, 2),
         factor = paste(site, ocean_depth, sep = "_"),
         factor_r = paste(site, ocean_depth, replica, sep = "_")) %>%
  dplyr::select(study, site, ocean_depth, factor, factor_r)
sample_data(ps_sfs) <- df_sfs

ps_sfs@refseq <- Biostrings::DNAStringSet(asv_sequences)
otu_table(ps_sfs) <- t(otu_table(ps_sfs))

ps_sfs <- clean_zero_reads(ps_sfs, "Specie")


#RAREFY--------------------------------------
sample_sums <- sample_sums(ps_sfs) #nolint
which.min(sample_sums)
r_ps_sfs <- rarefaction(ps_sfs)

##CURVAS DE RAREFACCION--------------------------------

calculate_rarefaction_curves2 <- function(psdata, measures, depths) {
  require("plyr") # ldply
  require("reshape2") # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # Calcular equitabilidad (evenness)
    evenness_values <- diversity(t(rarified_psdata@otu_table)) / log(specnumber(t(rarified_psdata@otu_table)))
    
    # Añadir la equitabilidad a los resultados
    evenness_df <- data.frame(Sample = rownames(alpha_diversity), Evenness = evenness_values)
    alpha_diversity <- cbind(alpha_diversity, Evenness = evenness_values)
    
    # Convertir los resultados a formato largo (melted)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity),
                                   varnames = c('Sample', 'Measure'),
                                   value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # Esto habilita la adición automática de la profundidad al output por ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # Convertir Depth de factor a numérico
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# Resumen de los resultados obtenidos en las curvas de rarefacción
curve_summary_verbose <- function(rarefied_ps, r_curve_data) {
  r_curve_data_summary <- ddply(r_curve_data,
                                c("Depth", "Sample", "Measure"),
                                summarise,
                                Alpha_diversity_mean = mean(Alpha_diversity), #nolint
                                Alpha_diversity_sd = sd(Alpha_diversity)) #nolint
  r_curve_data_summary_verbose <- merge(r_curve_data_summary,
                                        data.frame(sample_data(rarefied_ps)),
                                        by.x = "Sample", by.y = "row.names")
  return(r_curve_data_summary_verbose)
}


r_curve2 <- calculate_rarefaction_curves2(r_ps_sfs,
                                        c("Observed", "Shannon", "Chao1", "Evenness"),
                                        rep(c(1:150 * 100), each = 5))

# Data summary curve rarefy
r_curve_summary <- curve_summary_verbose(r_ps_sfs, r_curve2)


##plot y:reads x:sample, agglomerated by taxonomic level----
a <- ggplot(data = data.frame(x = 1:length(sample_sums(ps_sfs)), #nolint
                              y = sort(sample_sums(ps_sfs), decreasing = TRUE)), #nolint
            aes(x = x, y = y)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = " Reads for Sample",
       x = "Samples",
       y = "Reads")
a
#plot whit sample rarefy
b <- ggplot(data = data.frame(x = 1:length(sample_sums(r_ps_sfs)), #nolint
                              y = sort(sample_sums(r_ps_sfs), decreasing = TRUE)),
            aes(x = x, y = y)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Rarefy",
       x = "Samples",
       y = "Reads")
b
plot_rarefy <- a | b
#ggsave ("Plot_papper/reads_rarefy.png", plot_rarefy)

#plot whit alpha diversity obs and shannon
c <- r_curve_summary  %>% 
  filter(Measure != "se.chao1")  %>% 
  ggplot( aes(
              x = Depth, #nolint
              y = Alpha_diversity_mean, #nolint
              ymin = Alpha_diversity_mean - Alpha_diversity_sd, #nolint
              ymax = Alpha_diversity_mean + Alpha_diversity_sd, #nolint
              colour = as.factor(factor), #nolint
              group = Sample)) + #nolint
  geom_line(linewidth = 0.2, linetype = "solid") +
  facet_wrap(facets = ~ Measure, scales = "free_y") +
  labs(x = "Reads",
       y = "Alpha diversity mean") +
  scale_colour_manual(values = c("#e46385", "#17becf", "#8263e4", "#c5e465"), 
                    labels = c("Algarrobo 30m", "Algarrobo 60m", "Las Cruces 30m", "Las Cruces 60m")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), 
        legend.title = element_blank())
c
#ggsave("plot_papper/rarefaction_curve.png", c, width = 8, height = 4)


#ASV POR PHYLUM------------------
##contar numero de assigned by tax level-------
r_ps_sfs %>%
  clean_zero_reads(.,"Phylum" ) %>%
  subset_taxa(.,      Phylum != "Unassigned") %>%
  subset_taxa(.,      Phylum != "NA") %>%
  #subset_samples(., factor == "Algarrobo_30m") %>%
  ntaxa()

##contar numero de asv por phylum-----------
count_asv_phylum <- list()

r_ps_sfs_know <- subset_taxa(r_ps_sfs, Phylum != "Unassigned")
tax_phylum <- unique(tax_table(r_ps_sfs_know)[, "Phylum"])

for( i in tax_phylum) {
  resultado <- r_ps_sfs %>%
  subset_taxa(Phylum == i)  %>%
  prune_taxa(taxa_sums(.) > 0,
             .) %>%
  ntaxa()

 count_asv_phylum[[i]] <- resultado

 }

print(count_asv_phylum)

#ELEGIR PS A OCUPAR-------------------
r_ps_oc <- r_ps_sfs
#r_ps_oc <- r_ps_sfs_know


#ABUNDNACIA POR PHYLUM------------------------
phylum_colors <- c("Annelida" = "#9467bd",# Lavanda
                   "Arthropoda" = "#ff7f0e",  # Naranja
                   "Bryozoa" = "#2ca02c",  # Verde
                   "Cnidaria" = "#bcbd22",  # Verde lima
                   "Echinodermata" = "#c5b0d5",  # Morado
                   "Mollusca" = "#d62728",  # Rojo
                   "Nematoda" = "#17becf",  # Cian
                   "Nemertea" = "#f7b6d2",  # Rosa claro
                   "Porifera" = "#1f77b4",  # Azul
                   "Phatyhelminthes" = "#62162f",  # Marrón
                   "Rotifera" = "#e377c2", 
                   "Chordata" = "#8c564b",
                   "Unassigned" = "#7f7f7f")  # Gris

###BAR PLOT ABSOLUTE RELATIVE#
taxonomic_level <- "Phylum"
sfs_rel <- r_ps_oc %>%
  tax_glom(taxrank = taxonomic_level) %>% 
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%                         # Melt to long format
  arrange(Phylum)

str(sfs_rel)
sfs_rel %>%
  group_by(factor) %>%
  summarise(total_abundance = sum(Abundance))

sfs_ps <- sfs_rel %>%
  dplyr::select(Phylum, Sample, Abundance,
         ocean_depth, site) %>%
  group_by(Phylum, ocean_depth, site) %>%
  filter(Abundance > 0.01)

bar_plot <- ggplot(sfs_ps) +
  geom_col(mapping = aes(
                         x = as.factor(ocean_depth),
                         y = Abundance,
                         fill = !!sym(taxonomic_level)),
           position = "fill",
           show.legend = TRUE) +
  theme_bw() +
  labs(x = "Ocean Depth (m)", y = "Proportion of Community") +
  facet_grid( ~ site, 
              labeller = as_labeller(c(`Algarrobo`= "Algarrobo", 
                                       `Las_Cruces` = "Las Cruces"))) +
  theme(axis.title = element_text(size = 13),  # Tamaño de los títulos de los ejes
        axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 0.5),   # Tamaño de las etiquetas del eje X
        axis.text.y = element_text(size = 9), # Tamaño de las etiquetas del eje Y
        legend.title = element_text(size = 11), # Tamaño del título de la leyenda
        legend.text = element_text(size = 9.5),
        strip.text.x = element_text(size = 12),  # Tamaño del texto en el eje x
        strip.text.y = element_text(size = 12),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = phylum_colors) + 
  guides(fill = guide_legend(override.aes = list(size = 7.5)))
bar_plot
#ggsave("Plot_papper/bar_plot_relative.png", bar_plot, width = 5, height = 5)

#PRESENCIA Y AUSENCIA-----------
r_ps_oc_pa <- normalize_pa_ps(r_ps_oc)

#ALPHA DIVERSITy-----------------
#data frame con diversity
r_df_sfs <- data.frame(sample_data(r_ps_oc))
data_r_otu_sfs <- t(data.frame(otu_table(r_ps_oc)))
data_r_richness <- estimateR(data_r_otu_sfs)
S.evenness <- diversity(data_r_otu_sfs) / log(specnumber(data_r_otu_sfs))
S.shannon <- diversity(data_r_otu_sfs, index = "shannon")
sfs_alphadiv <- cbind(r_df_sfs, t(data_r_richness), S.shannon, S.evenness)
rm(r_df_sfs, data_r_otu_sfs, data_r_richness, S.evenness, S.shannon)

##grafico
D1 <- ggplot(sfs_alphadiv, aes(x = site, y = S.obs, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Richness", x =NULL, fill = "Depth",  tag = "A") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D2 <- ggplot(sfs_alphadiv, aes(x=site, y=S.chao1, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Chao1", x =NULL, fill = "Depth",  tag = "B") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D3 <- ggplot(sfs_alphadiv, aes(x=site, y=S.evenness, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Evenness", x ="Site", fill = "Depth",  tag = "C") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D4 <- ggplot(sfs_alphadiv, aes(x=site, y=S.shannon, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Shannon", x ="Site", fill = "Depth", tag = "D") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()

plot_div <-  D1 + D2 + D3 + D4 +
  plot_layout(guides = 'collect', axis_titles = "collect") & theme(
  plot.title = element_text(size = 24, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
  axis.title = element_text(size = 22),  # Tamaño de los títulos de los ejes
  legend.title = element_text(size = 23), # Tamaño del título de la leyenda
  legend.text = element_text(size = 17),
  legend.key.size = unit(2, "lines"),
  axis.text.x = element_text(size = 17),   # Tamaño de las etiquetas del eje X
  axis.text.y = element_text(size = 17),
  axis.title.y = element_text(margin = margin(r = 5)) # Aumenta el margen derecho del título del eje Y
  
  ) & 
  guides(fill = guide_legend(override.aes = list(size = 5)))# Ajustar el tamaño de los símbolos
plot_div
#ggsave("Plot_papper/alpha_div.png", plot_div, width = 10, height = 8)


#GLM------------
###ver familias----------
hist(sfs_alphadiv$S.shannon)
hist(sfs_alphadiv$S.chao1)
hist(sfs_alphadiv$S.evenness)
hist(sfs_alphadiv$S.obs)

norm <- fitdist   (sfs_alphadiv$S.shannon, "norm")
lnormal <- fitdist(sfs_alphadiv$S.shannon, "lnorm")
exp <- fitdist    (sfs_alphadiv$S.shannon, "exp")
gamma <- fitdist  (sfs_alphadiv$S.shannon, "gamma")

# Cumulative distribution frequency, Gráficos de comparación CDF
(CDF.dens=cdfcomp(list(exp, norm, lnormal, gamma),
                  addlegend=T,main="",legendtext=c("Exp", "Normal", "lnorm", "gamma"),
                  plotstyle = "ggplot")+
    xlab("Densidad")+
    geom_line(size=0.8)+
    theme(axis.title=element_text(size=8), 
          axis.text = element_text(size=10), 
          legend.position = c(0.7,0.45),
          legend.text=element_text(size=6)))

#Gráficos de comparación QQ
(QQ.ZA.dens=qqcomp(list(exp, norm, lnormal, gamma),addlegend=F,main="",legendtext=c("Exp", "Normal", "lnorm", "gamma"),plotstyle = "ggplot")+
    theme_bw()+
    geom_jitter(size=2, height=0.2)+
    geom_line()+
    theme(axis.title=element_text(size=18), 
          axis.text = element_text(size=16), 
          title=element_blank(),
          legend.position = c(0.25,0.75),
          legend.text=element_text(size=14)))

###Homocedasticidad------------
car::leveneTest(S.shannon ~ site * ocean_depth,  data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.evenness ~ site * ocean_depth, data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.obs ~ site * ocean_depth,      data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad
car::leveneTest(S.chao1 ~ site * ocean_depth,    data = sfs_alphadiv) # Usamos Levene porque es más robusto frente a distribuciones de variables aleatroias no normales. F max se basa en el supuesto de normalidad


##GLM calculos-----
glm_shannon <- glm(S.shannon ~ site * ocean_depth,
                   family = gaussian(link = identity), data = sfs_alphadiv)

glm_evenness <- betareg(S.evenness ~ site + ocean_depth + site * ocean_depth,
                        data = sfs_alphadiv, link = "logit")

# Richness (Poisson or Negative Binomial regression)
# First check for overdispersion
poisson_model <- glm(S.obs ~ site * ocean_depth, 
                     family = poisson(link = "log"), data = sfs_alphadiv)
if (summary(poisson_model)$dispersion > 1) {
  # Overdispersion detected, use Negative Binomial
  glm_richness <- glm.nb(S.obs ~ site * ocean_depth, data = sfs_alphadiv)
} else {
  glm_richness <- poisson_model
}

#mod_glmm <- glmer(S.obs ~ site * ocean_depth + (1 | site/ocean_depth/factor_r), #Revisar si son glm o glmm
#                family = poisson(link = "log"), data = sfs_alphadiv_pa)  #paquete library(lme4)


# Chao1 (Gaussian or Gamma regression)
# Assess distribution of Chao1 values
if (skewness(sfs_alphadiv$S.chao1) > 1) {
  # Use Gamma if skewed
  glm_chao1 <- glm(S.chao1 ~ site * ocean_depth,
                   family = Gamma(link = "log"), data = sfs_alphadiv)
} else {
  # Use Gaussian if not skewed
  glm_chao1 <- glm(S.chao1 ~ site * ocean_depth,
                   family = gaussian(link = "identity"), data = sfs_alphadiv)
}

summary(glm_shannon)
summary(glm_evenness)
summary(glm_richness)
summary(glm_chao1)
summary(mod_glmm)



#NMSD--------------------
#ordination
set.seed(1)
all_nmds_pa <- ordinate(
  physeq = r_ps_oc_pa,
  method = "NMDS",
  distance = "jaccard"
)

#pasar resultados a un data frame 
all_nmds_points_pa <- as.data.frame(all_nmds_pa$points)
r_sfs_df <- data.frame(sample_data(r_ps_sfs))

all_nmds_df_pa <- all_nmds_points_pa  %>% 
mutate(site = r_sfs_df$site, 
       ocean_depth = r_sfs_df$ocean_depth, 
       type_of_sample = r_sfs_df$type_of_sample, 
       factor = paste(r_sfs_df$site, r_sfs_df$ocean_depth, sep = "_"))

                          
#plot jaccard
nmds_jaccard <- ggplot(data = all_nmds_df_pa,
                    aes(x = MDS1,
                        y = MDS2
                    )) + 
  stat_ellipse(type = "norm", linetype = 2,
               show.legend = FALSE, size = 0.6, alpha = 0.9, 
               aes(group = ocean_depth, color = ocean_depth)) +
  geom_point(size = 5, aes(shape = site, color = ocean_depth)) +
  scale_shape_manual(values = c(15, 19), 
                     labels = c(`Algarrobo` = "Algarrobo", `Las_Cruces` = "Las Cruces")) +
  scale_color_manual(values = c("#17becf", "#e46385")) +
  labs(color = "Depth", shape = "Site", x = "nMDS1", y = "nMDS2") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 17, angle = 0, vjust = 0.5, hjust = 0.5), 
        axis.title = element_text(size = 19),  # Tamaño de los títulos de los ejes
        axis.text.y = element_text(size = 17),  # Tamaño de las etiquetas del eje Y 
        strip.text.x = element_text(size = 15),  # Tamaño del texto en el eje x
        strip.text.y = element_text(size = 15),  # Tamaño del texto en el eje y
        legend.title = element_text(size = 17), # Tamaño del título de la leyenda
        legend.text = element_text(size = 15)) +# Ajustar el ancho de los elementos de la leyenda
  guides(fill = guide_legend(override.aes = list(size = 16)))
nmds_jaccard
#ggsave("Plot_papper/ndms_jaccard.png", nmds_jaccard, width = 7, height = 5)


##PERMANOVA----------------
#test adonis

df_sfs_pa <- data.frame(sample_data(r_ps_oc_pa))
dist_sfs_pa <- phyloseq::distance(r_ps_oc_pa, method = "jaccard")
test_pa <- adonis2(dist_sfs_pa ~ (site + ocean_depth + site * ocean_depth),  data=df_sfs, permutations = 1e3)
test_pa

#Turneover y Nestedness-------------
# Load necessary packages
library(phyloseq)
library(betapart)
library(pheatmap)
library(gridExtra)
library(grid)

# Assuming your phyloseq object is called `physeq`
physeq <- r_ps_oc # replace with your actual phyloseq object

# Merge samples by site and depth
merged_physeq <- merge_samples(r_ps_oc_pa, group = "factor")

#PA
physeq_pa <- normalize_pa_ps(merged_physeq)

# Extract the OTU table from the merged phyloseq object
otu_table_merged <- otu_table(physeq_pa)

# Compute beta diversity components
beta_div <- beta.pair(otu_table_merged, index.family = "sorensen")

# Extract turnover and nestedness
turnover <- beta_div$beta.sim
nestedness <- beta_div$beta.sne

# Print results
print(turnover)
print(nestedness)

# Compute beta diversity components using Jaccard index family
beta_div_jaccard <- beta.pair(otu_table_merged, index.family = "jaccard")

# Extract turnover and nestedness for Jaccard
turnover_jaccard <- as.matrix(beta_div_jaccard$beta.jtu)
nestedness_jaccard <- as.matrix(beta_div_jaccard$beta.jne)

# Print results for Jaccard
print(turnover_jaccard)
print(nestedness_jaccard)

# Eliminar los valores por encima de la diagonal
turnover_jaccard[upper.tri(turnover_jaccard)] <- NA
nestedness_jaccard[lower.tri(nestedness_jaccard)] <- NA

#heatmap 
new_labels <- c("AL30", "AL60", "LC30", "LC60")
 
# Crear una nueva matriz vacía para almacenar los valores combinados
combined_matrix <- matrix(NA, nrow = nrow(turnover_jaccard), ncol = ncol(turnover_jaccard))

# Poner turnover en la mitad inferior
combined_matrix[lower.tri(turnover_jaccard)] <- turnover_jaccard[lower.tri(turnover_jaccard)]

# Poner nestedness en la mitad superior
combined_matrix[upper.tri(nestedness_jaccard)] <- nestedness_jaccard[upper.tri(nestedness_jaccard)]

heatmap_combined <- pheatmap(combined_matrix, 
                             cluster_rows = FALSE, 
                             cluster_cols = FALSE,
                             display_numbers = TRUE, 
                             color = colorRampPalette(c("#F4FAFE", "firebrick"))(100), 
                             labels_row = new_labels,  # Etiquetas personalizadas para filas
                             labels_col = new_labels, 
                             angle_col = 0, 
                             legend = TRUE, 
                             fontsize_number = 12, 
                             fontsize_row = 12, 
                             fontsize_col = 12,
                             number_color = "black", 
                             na_col = "white", 
                             main = "")

#ggsave("Plot_papper/heatmap_combined.png", heatmap_combined, width = 10, height = 5)  # Ajusta el tamaño según lo necesites

#Beta diversity ----------------------------------
##subset por site and depth------
r_ps_AL30 <- subset_samples(r_ps_oc, factor == "Algarrobo_30m")
r_ps_AL60 <- subset_samples(r_ps_oc, factor == "Algarrobo_60m")
r_ps_LC30 <- subset_samples(r_ps_oc, factor == "Las_Cruces_30m")
r_ps_LC60 <- subset_samples(r_ps_oc, factor == "Las_Cruces_60m")

taxonomic_level <- "Specie"
r_ps_AL30 <- clean_zero_reads(r_ps_AL30, taxonomic_level)
r_ps_AL60 <- clean_zero_reads(r_ps_AL60, taxonomic_level)
r_ps_LC30 <- clean_zero_reads(r_ps_LC30, taxonomic_level)
r_ps_LC60 <- clean_zero_reads(r_ps_LC60, taxonomic_level)


##Beta diversity-------
library(ggvenn)
library(ggVennDiagram)


tax_table_AL30 <- data.frame(tax_table(r_ps_AL30))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_AL60 <- data.frame(tax_table(r_ps_AL60))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_LC30 <- data.frame(tax_table(r_ps_LC30))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_LC60 <- data.frame(tax_table(r_ps_LC60))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))


tax_table_AL30_c <- as.character(tax_table_AL30$all_taxa_level)
tax_table_AL60_c <- as.character(tax_table_AL60$all_taxa_level)
tax_table_LC30_c <- as.character(tax_table_LC30$all_taxa_level)
tax_table_LC60_c <- as.character(tax_table_LC60$all_taxa_level)

phylum_AL30_c <- as.character(tax_table_AL30$k_p)
phylum_AL60_c <- as.character(tax_table_AL60$k_p)
phylum_LC30_c <- as.character(tax_table_LC30$k_p)
phylum_LC60_c <- as.character(tax_table_LC60$k_p)

class_AL30_c <- as.character(tax_table_AL30$k_c)
class_AL60_c <- as.character(tax_table_AL60$k_c)
class_LC30_c <- as.character(tax_table_LC30$k_c)
class_LC60_c <- as.character(tax_table_LC60$k_c)

order_AL30_c <- as.character(tax_table_AL30$k_o)
order_AL60_c <- as.character(tax_table_AL60$k_o)
order_LC30_c <- as.character(tax_table_LC30$k_o)
order_LC60_c <- as.character(tax_table_LC60$k_o)

family_AL30_c <- as.character(tax_table_AL30$k_f)
family_AL60_c <- as.character(tax_table_AL60$k_f)
family_LC30_c <- as.character(tax_table_LC30$k_f)
family_LC60_c <- as.character(tax_table_LC60$k_f)


asv_beta <- list(
  AL60 = tax_table_AL60_c, AL30 = tax_table_AL30_c, 
  LC30 = tax_table_LC30_c, LC60 = tax_table_LC60_c)

phylum_beta <- list(
  AL30 = tax_table_AL30_c, AL60 = tax_table_AL60_c,
  LC30 = tax_table_LC30_c, LC60 = tax_table_LC60_c)

class_beta <- list(
  AL30 = tax_table_AL30_c, AL60 = tax_table_AL60_c,
  LC30 = tax_table_LC30_c, LC60 = tax_table_LC60_c)

order_beta <- list(
  AL30 = tax_table_AL30_c, AL60 = tax_table_AL60_c,
  LC30 = tax_table_LC30_c, LC60 = tax_table_LC60_c)

family_beta <- list(
  AL30 = tax_table_AL30_c, AL60 = tax_table_AL60_c,
  LC30 = tax_table_LC30_c, LC60 = tax_table_LC60_c)


all_venn <- ggvenn(asv_beta,
                   columns = NULL,
                   #fill_color = c("#4daf4a", "#8263e4"),
                   stroke_size = 0.6,
                   set_name_size = 9,
                   text_size = 6,
                   #auto_scale = TRUE
                   ) +
  theme_void()
all_venn
all_venn <- all_venn + expand_limits(x = c(0, 1), y = c(0, 1))
all_venn
#ggsave("Plot_papper/venn_all.png", all_venn, width = 10, height = 8)

library(VennDiagram)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


library("ggVennDiagram")
beta_asv_plot <- ggVennDiagram(asv_beta, label_alpha = 0, edge_size = .8) +
  scale_fill_gradient(low = "#F4FAFE", high = "firebrick") +
  labs(fill = "Count") +
  theme(plot.tag = element_text(size = 4), 
        plot.subtitle = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 7),            # Ajustar el tamaño de la leyenda
        legend.title = element_text(face = "bold"), 
        legend.position = "bottom",
        plot.margin = margin(15, 15, 15, 15)
  )   # Leyenda en negrita
#ggsave("Plot_papper/venn_asv.png", beta_asv_plot, width = 10, height = 8)


###------------PLOT combined beta----
plot_beta <- arrangeGrob(
  beta_asv_plot, heatmap_combined[[4]],
  ncol = 2,
  widths = c(1, 1)  # Cambia los números según el ajuste deseado
)
#ggsave("Plot_papper/plot_beta_asv2.png", plot_beta, width = 12, height = 6)



###------------subset por site and depth------
r_ps_AL <- subset_samples(r_ps_oc, site == "Algarrobo")
r_ps_LC <- subset_samples(r_ps_oc, site == "Las_Cruces")
r_ps_30 <- subset_samples(r_ps_oc, ocean_depth == "30m")
r_ps_60 <- subset_samples(r_ps_oc, ocean_depth == "60m")

taxonomic_level <- "Specie"
r_ps_AL <- clean_zero_reads(r_ps_AL, taxonomic_level)
r_ps_LC <- clean_zero_reads(r_ps_LC, taxonomic_level)
r_ps_30 <- clean_zero_reads(r_ps_30, taxonomic_level)
r_ps_60 <- clean_zero_reads(r_ps_60, taxonomic_level)


###---------------Beta diversity-------
library(ggvenn)
library(ggVennDiagram)


tax_table_AL <- data.frame(tax_table(r_ps_AL30))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_LC <- data.frame(tax_table(r_ps_AL60))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_30 <- data.frame(tax_table(r_ps_LC30))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))

tax_table_60 <- data.frame(tax_table(r_ps_LC60))  %>%
  mutate(all_taxa_level =  paste(Kingdom, Phylum, Class, Order,
                                 Family, Genus, Specie, sep = "_"),
         k_p = paste(Kingdom, Phylum, sep = "_"),
         k_c = paste(Kingdom, Phylum, Class, sep = "_"),
         k_o = paste(Kingdom, Phylum, Class, Order, sep = "_"),
         k_f = paste(Kingdom, Phylum, Class, Order, Family, sep = "_"))


tax_table_AL_c <- as.character(tax_table_AL$all_taxa_level)
tax_table_LC_c <- as.character(tax_table_LC$all_taxa_level)
tax_table_30_c <- as.character(tax_table_30$all_taxa_level)
tax_table_60_c <- as.character(tax_table_60$all_taxa_level)

phylum_AL_c <- as.character(tax_table_AL$k_p)
phylum_LC_c <- as.character(tax_table_LC$k_p)
phylum_30_c <- as.character(tax_table_30$k_p)
phylum_60_c <- as.character(tax_table_60$k_p)

class_AL_c <- as.character(tax_table_AL$k_c)
class_LC_c <- as.character(tax_table_LC$k_c)
class_30_c <- as.character(tax_table_30$k_c)
class_60_c <- as.character(tax_table_60$k_c)

order_AL_c <- as.character(tax_table_AL$k_o)
order_LC_c <- as.character(tax_table_LC$k_o)
order_30_c <- as.character(tax_table_30$k_o)
order_60_c <- as.character(tax_table_60$k_o)

family_AL_c <- as.character(tax_table_AL$k_f)
family_LC_c <- as.character(tax_table_LC$k_f)
family_30_c <- as.character(tax_table_30$k_f)
family_60_c <- as.character(tax_table_60$k_f)


asv_beta_site <- list(
  AL = tax_table_AL_c, LC = tax_table_LC_c)

asv_beta_depth <- list(
  Shallow = tax_table_30_c, Mesophotic = tax_table_60_c)

all_venn_site <- ggvenn(asv_beta_site,
                   columns = NULL,
                   #fill_color = c("#4daf4a", "#8263e4"),
                   stroke_size = 0.6,
                   set_name_size = 9,
                   text_size = 6,
                   #auto_scale = TRUE
) +
  theme_void()
all_venn_site
all_venn_site <- all_venn_site + expand_limits(x = c(0, 1), y = c(0, 1))
all_venn_site
#ggsave("Plot_papper/ven_all_site.png", all_venn_site, width = 12, height = 6)

all_venn_depth <- ggvenn(asv_beta_depth,
                        columns = NULL,
                        #fill_color = c("#4daf4a", "#8263e4"),
                        stroke_size = 0.6,
                        set_name_size = 9,
                        text_size = 6,
                        #auto_scale = TRUE
) +
  theme_void()
all_venn_depth
all_venn_depth <- all_venn_depth + expand_limits(x = c(0, 1), y = c(0, 1))
all_venn_depth
#ggsave("Plot_papper/ven_all_depth.png", all_venn_depth, width = 12, height = 6)
