###
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

ps_sfs <- clean_zero_reads(ps_sfs, "Specie")

###-----------------------RAREFY--------------------------------------
#rarefy
sample_sums <- sample_sums(ps_sfs) #nolint
which.min(sample_sums)
r_ps_sfs <- rarefaction(ps_sfs)

#
asv_transpuesta <- t(data.frame(otu_table(ps_sfs))) #cambiar filas por columnas, ahora asv en col y sample in rows
plot_curve_r <- rarecurve(asv_transpuesta, step = 100, cex = .5, las = .2,
                          xlab = "Number of reads",
                          ylab = "Number of OTU observed")


###--------CURVAS DE RAREFACCION--------------------------------

calculate_rarefaction_curves2 <- function(psdata, measures, depths) {
  require("plyr") # ldply
  require("reshape2") # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # Calcular equitabilidad (evenness)
    evenness_values <- diversity(rarified_psdata@otu_table) / log(specnumber(rarified_psdata@otu_table))
    
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
                                        rep(c(1:1500 * 100), each = 5))

# Data summary curve rarefy
r_curve_summary <- curve_summary_verbose(r_ps_sfs, r_curve2)

##

#plot y:reads x:sample, agglomerated by taxonomic level
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
ggsave ("Plot_papper/reads_rarefy.png", plot_rarefy)

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
       y = "Alpha diversity mean",
       title = "Rarefaction curve",
       colour = "Factor") +
  scale_colour_manual(values = c("#e46385", "#17becf", "#8263e4", "#c5e465"), 
                    labels = c("Algarrobo 30m", "Algarrobo 60m", "Las Cruces 30m", "Las Cruces 60m")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
c
#ggsave("plot_papper/rarefaction_curve2.png", c, width = 8, height = 4)


#contar numero de ASV unassigned--------------
ps_sfs <- prune_taxa(taxa_sums(ps_sfs) > 0,
                 ps_sfs)

count_total <- ntaxa(ps_sfs)

count_know <- ps_sfs %>%
  subset_taxa(Phylum != "Unassigned")  %>% 
  subset_taxa(Phylum != "NA")  %>% 
  ntaxa()

100 - ((count_know*100)/count_total)

#aglomerar por specie, contar nuemro de ASV y nuemeor de NA
ps <- clean_zero_reads(ps, "Specie")

###---------------ASV POR PHYLUM------------------
#contar numero de asv por tax level
r_ps_sfs %>%
  clean_zero_reads(.,"Phylum" ) %>%
  subset_taxa(.,      Phylum != "Unassigned") %>%
  subset_taxa(.,      Phylum != "NA") %>%
  ntaxa()


##contar numero de asv por phylum
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

###-----------ELEGIR PS A OCUPAR-------------------
r_ps_oc <- r_ps_sfs
#r_ps_oc <- r_ps_sfs_know


###--------ABUNDNACIA POR PHYLUM------------------------

###BAR PLOT ABSOLUTE ABUNDANCE#
taxonomic_level <- "Phylum"
sfs_ab_absol <- r_ps_oc %>%
  tax_glom(taxrank = taxonomic_level) %>% #se aliminan NA
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%                         # Melt to long format
  arrange(Phylum)

sfs_absol_ps <- sfs_ab_absol %>%
  dplyr::select(Phylum, Sample, Abundance,
         ocean_depth, site) %>%
  group_by(Phylum, ocean_depth, site) %>%
  filter(Abundance > 0.01)

bar_plot_absolute <- ggplot(sfs_absol_ps) +
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
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
        strip.text.x = element_text(size = 12),  # Tamaño del texto en el eje x
        strip.text.y = element_text(size = 12)) +  # Tamaño del texto en el eje y
  theme(axis.title = element_text(size = 13),  # Tamaño de los títulos de los ejes
        axis.text.x = element_text(size = 9),   # Tamaño de las etiquetas del eje X
        axis.text.y = element_text(size = 9), # Tamaño de las etiquetas del eje Y
        legend.title = element_text(size = 11), # Tamaño del título de la leyenda
        legend.text = element_text(size = 9.5),
        text = element_text(family = "Arial")) +
  guides(fill = guide_legend(override.aes = list(size = 7.5)))
bar_plot_absolute
#ggsave("Plot_papper/bar_plot_ab_absolute.png", bar_plot_absolute, width = 5, height = 5)


###---------presencia y ausencia-----------
r_ps_oc_pa <- normalize_pa_ps(r_ps_oc)

###-----ABUNDANCIA POR PHYLUM PA----------
Phyfrac <- as.character(tax_table(r_ps_oc_pa)[, "Phylum"])
Phyfrac <- ifelse(Phyfrac == "", "Unassigned", Phyfrac)
Phyfrac <- factor(Phyfrac)

OTUtab_a_unique = apply(otu_table(r_ps_oc_pa), MARGIN = 1, function(x) {
  tapply(x, INDEX = Phyfrac, FUN = sum, na.rm = F, simplify = TRUE)
})
OTUtab_a_unique<- as.data.frame(as.table(OTUtab_a_unique))
OTUtab_a_unique <- OTUtab_a_unique %>%
  dplyr::rename(Phylum = Var1) %>%
  mutate(site = rep(c("Algarrobo", "Las_Cruces"), each = 96), 
  ocean_depth = rep(c("30m", "60m", "30m", "60m"),each = 48)) 


bar_plot_pa <- ggplot(OTUtab_a_unique) +
  geom_col(mapping = aes(
    x = as.factor(ocean_depth),
    y = Freq,
    fill = Phylum),
    position = "fill") +
  theme_bw() +
  labs(x = "Ocean Depth (m)", y = "Proportion of Community", fill = "Phyla") +
  facet_grid(~ site, 
             labeller = as_labeller(c(`Algarrobo` = "Algarrobo", 
                                      `Las_Cruces` = "Las Cruces"))) +
  theme(axis.text.x = element_text(size= 9, angle = 0, vjust = 0.5, hjust = 0.5), 
        strip.text.x = element_text(size = 12),  # Tamaño del texto en el eje x
        strip.text.y = element_text(size = 12),  # Tamaño del texto en el eje y
        axis.title = element_text(size = 13),   # Tamaño de las etiquetas del eje X
        axis.text.y = element_text(size = 9),  # Tamaño de las etiquetas del eje Y
        legend.title = element_text(size = 11),  # Tamaño del título de la leyenda
        legend.text = element_text(size = 9.5)) +
  guides(fill = guide_legend(override.aes = list(size = 5.5)))
bar_plot_pa
#ggsave("Plot_papper/bar_plot_pa_asv.png", bar_plot_pa, width = 5, height = 5)


###--------------ALPHA DIVERSITy-----------------

#data frame con diversity
r_df_sfs <- data.frame(sample_data(r_ps_oc))
data_r_otu_sfs <- (data.frame(otu_table(r_ps_oc)))
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

plot_div <-  D1 + D2 + D3 + D4 + 
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
    axis.title = element_text(size = 22),  # Tamaño de los títulos de los ejes
    legend.title = element_text(size = 23), # Tamaño del título de la leyenda
    legend.text = element_text(size = 17),
    legend.key.size = unit(2, "lines"),
    axis.text.x = element_text(size = 17),   # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 17),
    axis.title.y = element_text(margin = margin(r = 5)) # Aumenta el margen derecho del título del eje Y
    
  ) & 
  guides(fill = guide_legend(override.aes = list(size = 5)))
plot_div
#ggsave("Plot_papper/alpha_div2.png", plot_div, width = 10, height = 8)

##alpha div pa-------
df_r <- data.frame(sample_data(r_ps_oc_pa))
data_r_otu <- data.frame(otu_table(r_ps_oc_pa))
data_r_richness <- estimateR(data_r_otu)
S.evenness <- diversity(data_r_otu) / log(specnumber(data_r_otu)) # calculate evenness index using vegan package
S.shannon <- diversity(data_r_otu, index = "shannon") # calculate Shannon index using vegan package
sfs_alphadiv_pa <- cbind(df_r, t(data_r_richness), S.shannon, S.evenness) # combine all indices in one data table
rm(df_r, data_r_otu, data_r_richness, S.evenness, S.shannon)

D1 <- ggplot(sfs_alphadiv_pa, aes(x = site, y = S.obs, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = 'Richness', y = "Alpha Diversity", x =NULL, fill = "Depth",  tag = "A") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D2 <- ggplot(sfs_alphadiv_pa, aes(x=site, y=S.chao1, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = 'Chao1', y = NULL, x =NULL, fill = "Depth",  tag = "B") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D3 <- ggplot(sfs_alphadiv_pa, aes(x=site, y=S.evenness, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = 'Evenness', y = "Alpha Diversity", x ="Site", fill = "Depth",  tag = "C") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()
D4 <- ggplot(sfs_alphadiv_pa, aes(x=site, y=S.shannon, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = 'Shannon', y = NULL, x ="Site", fill = "Depth", tag = "D") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  scale_x_discrete(labels = c("Algarrobo" = "Algarrobo", "Las_Cruces" = "Las Cruces")) +
  theme_bw()

plot_div <-  D1 + D2 + D3 + D4 +
  plot_layout(guides = 'collect', axis_titles = "collect") &
  labs(y = "Alpha diversity", tag = "") & theme(
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
#ggsave("Plot_papper/alpha_div_pa.png", plot_div, width = 10, height = 8)

plot_div <-  D1 + D2 + D3 + D4 + 
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
    axis.title = element_text(size = 22),  # Tamaño de los títulos de los ejes
    legend.title = element_text(size = 23), # Tamaño del título de la leyenda
    legend.text = element_text(size = 17),
    legend.key.size = unit(2, "lines"),
    axis.text.x = element_text(size = 17),   # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 17),
    axis.title.y = element_text(margin = margin(r = 5)) # Aumenta el margen derecho del título del eje Y
    
  ) & 
  guides(fill = guide_legend(override.aes = list(size = 5)))
plot_div
#ggsave("Plot_papper/alpha_div_pa2.png", plot_div, width = 10, height = 8)


###-------------------------GLM------------
##ver familias
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


##GLM
glm_shannon <- glm(S.shannon ~ site + ocean_depth + site * ocean_depth,
                   family = gaussian(link = identity), data = sfs_alphadiv)

glm_evenness <- betareg(S.evenness ~ site + ocean_depth + site * ocean_depth,
                        data = sfs_alphadiv, link = "logit")

# Richness (Poisson or Negative Binomial regression)
# First check for overdispersion
poisson_model <- glm(S.obs ~ ocean_depth + site + site * ocean_depth, 
                     family = poisson(link = "log"), data = sfs_alphadiv)
if (summary(poisson_model)$dispersion > 1) {
  # Overdispersion detected, use Negative Binomial
  glm_richness <- glm.nb(S.obs ~ ocean_depth + site + site * ocean_depth, data = sfs_alphadiv)
} else {
  glm_richness <- poisson_model
}

poisson_model_pa <- glm(S.obs ~ ocean_depth + site + site * ocean_depth, 
                     family = poisson(link = "log"), data = sfs_alphadiv_pa)
if (summary(poisson_model)$dispersion > 1) {
  # Overdispersion detected, use Negative Binomial
  glm_richness_pa <- glm.nb(S.obs ~ ocean_depth + site + site * ocean_depth, data = sfs_alphadiv_pa)
} else {
  glm_richness_pa <- poisson_model_pa
}



# Chao1 (Gaussian or Gamma regression)
# Assess distribution of Chao1 values
if (skewness(sfs_alphadiv$S.chao1) > 1) {
  # Use Gamma if skewed
  glm_chao1 <- glm(S.chao1 ~ ocean_depth + site + site * ocean_depth,
                   family = Gamma(link = "log"), data = sfs_alphadiv)
} else {
  # Use Gaussian if not skewed
  glm_chao1 <- glm(S.chao1 ~ ocean_depth + site + site * ocean_depth,
                   family = gaussian(link = "identity"), data = sfs_alphadiv)
}



summary(glm_shannon)
summary(glm_evenness)
summary(glm_richness)
summary(glm_richness_pa)
summary(glm_chao1)



###-------------------------NMSD--------------------
#ordination
set.seed(1)
all_pca <- ordinate(
  physeq = r_ps_oc,
  method = "NMDS",
  distance = "bray"
)

all_pca_pa <- ordinate(
  physeq = r_ps_oc_pa,
  method = "NMDS",
  distance = "jaccard"
)

#pasar resultados a un data frame 
all_pca_points <- as.data.frame(all_pca$points)
all_pca_points_pa <- as.data.frame(all_pca_pa$points)
r_sfs_df <- data.frame(sample_data(r_ps_sfs))

#agregar variables
all_pca_df <- all_pca_points  %>% 
mutate(site = r_sfs_df$site, 
       ocean_depth = r_sfs_df$ocean_depth, 
       type_of_sample = r_sfs_df$type_of_sample, 
       factor = paste(r_sfs_df$site, r_sfs_df$ocean_depth, sep = "_"))

all_pca_df_pa <- all_pca_points_pa  %>% 
mutate(site = r_sfs_df$site, 
       ocean_depth = r_sfs_df$ocean_depth, 
       type_of_sample = r_sfs_df$type_of_sample, 
       factor = paste(r_sfs_df$site, r_sfs_df$ocean_depth, sep = "_"))

##grafico
nmds_bray <- ggplot(data = all_pca_df,
                     aes(x = MDS1,
                         y = MDS2
                     )) + 
  stat_ellipse(type = "norm", linetype = 2,
               show.legend = FALSE, size = 0.6, alpha = 0.9, aes(color = ocean_depth)) +
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
nmds_bray
#ggsave("Plot_papper/ndms_bray.png", nmds_bray, width = 7, height = 5)

#plot jaccard
nmds_jaccard <- ggplot(data = all_pca_df_pa,
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

###--------------------PERMANOVA----------------
#test adonis
dist_sfs_r <- phyloseq::distance(r_ps_oc, method = "bray")
test <- adonis2(dist_sfs_r ~ (site + ocean_depth + site * ocean_depth),  data = df_sfs, permutations = 1e3)
test

df_sfs_pa <- data.frame(sample_data(r_ps_oc_pa))
dist_sfs_pa <- phyloseq::distance(r_ps_oc_pa, method = "jaccard")
test_pa <- adonis2(dist_sfs_pa ~ (site + ocean_depth + site * ocean_depth),  data=df_sfs, permutations = 1e3)
test_pa


###-----------------------SIMPER------------
##table

asv_table_r <- r_ps_oc %>% tax_glom("Order", NArm = FALSE) %>%
  otu_table() %>% as.data.frame()

asv_table_r_pa <- r_ps_oc_pa %>% tax_glom("Order", NArm = FALSE) %>%
  otu_table() %>% as.data.frame()

asv_table_r <- as.data.frame(otu_table(r_ps_oc))
tax_table_r <- as.data.frame(tax_table(r_ps_oc))
asv_table_r_pa <- as.data.frame(otu_table(r_ps_oc_pa))
dist_sfs <- vegdist(asv_table_r, method = "bray")
dist_sfs_pa <- vegdist(asv_table_r_pa, method = "jaccard")
dist_sfs_nmds <- metaMDS(dist_sfs)
dist_sfs_pa_nmds <- metaMDS(dist_sfs_pa)

##

simper_ocean <- simper(asv_table_r_pa, df_sfs$ocean_depth )
simper_ocean$`30m_60m`
summary(simper_ocean$`30m_60m`)
simper_ocean_df <- data.frame(
  ASV = simper_ocean$`30m_60m`$species,
  avegare = simper_ocean$`30m_60m`$average,
  sd = simper_ocean$`30m_60m`$sd,
  ord = simper_ocean$`30m_60m`$ord,
  p = simper_ocean$`30m_60m`$p
)
summary(simper_ocean)
simper_ocean_df <- merge(simper_ocean_df, tax_table_r, by = "row.names")
#write.csv(simper_ocean_df, "Plot_papper/simper_ocean_pa_order.csv", row.names = FALSE)


simper_site <- simper(asv_table_r_pa, df_sfs$site )
simper_site$Algarrobo_Las_Cruces
summary(simper_site$Algarrobo_Las_Cruces)
simper_site_df <- data.frame(
  ASV = simper_site$Algarrobo_Las_Cruces$species,
  avegare = simper_site$Algarrobo_Las_Cruces$average,
  sd = simper_site$Algarrobo_Las_Cruces$sd,
  ord = simper_site$Algarrobo_Las_Cruces$ord,
  p = simper_site$Algarrobo_Las_Cruces$p
)
summary(simper_site)
simper_site_df <- merge(simper_site_df, tax_table_r, by = "row.names")
#write.csv(simper_site_df, "Plot_papper/simper_site_pa.csv", row.names = FALSE)


###----Subset de filos mas ricos####
anelidos <- subset_taxa(r_ps_oc, Phylum == "Annelida")
artro    <- subset_taxa(r_ps_oc, Phylum == "Arthropoda")
briozoa  <- subset_taxa(r_ps_oc, Phylum == "Bryozoa")
cnidaria <- subset_taxa(r_ps_oc, Phylum == "Cnidaria")

anelidos <- clean_zero_reads(anelidos, "Specie")
artro    <- clean_zero_reads(artro, "Specie")
briozoa  <- clean_zero_reads(briozoa, "Specie")
cnidaria <- clean_zero_reads(cnidaria, "Specie")


calculate_alphadiv <- function(ps_subset) {
  
  r_df <- (data.frame(sample_data(ps_subset)))
  ps_otu_t <- (data.frame(otu_table(ps_subset)))
  data_r_richness <- estimateR(ps_otu_t)  
  S.evenness <- diversity(ps_otu_t) / log(specnumber(ps_otu_t))
  S.shannon <- diversity(ps_otu_t, index = "shannon")
  r_alphadiv <- cbind(r_df, t(data_r_richness), S.shannon, S.evenness)
  rm(data_r_richness, S.evenness, S.shannon)
  
  return(r_alphadiv)
  
}

anelidos_alphadiv <- calculate_alphadiv(anelidos)
artro_alphadiv    <- calculate_alphadiv(artro)
bryozoa_alphadiv  <- calculate_alphadiv(briozoa)
cnidaria_alphadiv <- calculate_alphadiv(cnidaria)


graphic_function <- function(data_alphadiv, grupo_tax) {
  
  A1 <- ggplot(data_alphadiv, aes(x=ocean_depth, y=S.obs, fill =site)) +
    geom_boxplot() +
    labs(title= grupo_tax, x= 'Ocean Depth', y= 'Richness', fill = "Site") +
    scale_fill_manual(values = c("#e46385", "#17becf"), 
                      labels = c("Algarrobo", "Las Cruces")) +
    ylim(c(0,20.5)) +
    theme_bw()
  A2 <- ggplot(data_alphadiv, aes(x=ocean_depth, y=S.chao1, fill =site)) +
    geom_boxplot() +
    labs(title= NULL, x= 'Ocean Depth', y= 'Chao1', fill = "Site") +
    scale_fill_manual(values = c("#e46385", "#17becf"),
                      labels = c("Algarrobo", "Las Cruces")) +
    ylim(c(0, 25.5)) +
    theme_bw()
  A3 <- ggplot(data_alphadiv, aes(x=ocean_depth, y=S.evenness, fill =site)) +
    geom_boxplot() +
    labs(title= NULL, x= 'Ocean Depth', y= 'Evenness', fill = "Site") +
    scale_fill_manual(values = c("#e46385", "#17becf"),
                      labels = c("Algarrobo", "Las Cruces")) +
    ylim(c(0, 1.5)) +
    theme_bw()
  A4 <- ggplot(data_alphadiv, aes(x=ocean_depth, y=S.shannon, fill =site)) +
    labs(title= NULL, x= 'Ocean Depth', y= 'Shannon', fill = "Site") +
    geom_boxplot() +
    scale_fill_manual(values = c("#e46385", "#17becf"), 
                      labels = c("Algarrobo", "Las Cruces")) +
    ylim(c(0, 2)) +
    theme_bw()
  

  return(list(A1 = A1, A2 = A2, A3 = A3, A4 = A4))
}

plot_anelidos <- graphic_function(anelidos_alphadiv, "Annelida")
plot_artro <- graphic_function(artro_alphadiv, "Arthropoda")
plot_bryozoa <- graphic_function(bryozoa_alphadiv, "Bryozoa")
plot_cnidaria <- graphic_function(cnidaria_alphadiv, "Cnidaria")

plot_div_group_tax <- plot_anelidos$A1 + plot_artro$A1 + plot_bryozoa$A1 + plot_cnidaria$A1 +
  plot_anelidos$A2 + plot_artro$A2 + plot_bryozoa$A2 + plot_cnidaria$A2 +
  plot_anelidos$A3 + plot_artro$A3 + plot_bryozoa$A3 + plot_cnidaria$A3 +
  plot_anelidos$A4 + plot_artro$A4 + plot_bryozoa$A4 + plot_cnidaria$A4 +
  plot_layout(guides = 'collect', axis_titles = "collect") & 
  theme(
    plot.title = element_text(size = 12.5, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
    axis.title = element_text(size = 12.5),  # Tamaño de los títulos de los ejes
    legend.title = element_text(size = 13.5), # Tamaño del título de la leyenda
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(2, "lines"),
    axis.text.x = element_text(size = 7.5),   # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 7.5)) &
  guides(fill = guide_legend(override.aes = list(size = 6)))# Ajustar el tamaño de los símbolos

ggsave("Plot_papper/group_tax_alpha_div.png", plot_div_group_tax, width = 10, height = 8)


calculate_glm <- function(subset_data_div) {
    
  glm_obs <- glm(S.obs ~ site + ocean_depth + site * ocean_depth,
                 family = gaussian(link = identity),
                 data = subset_data_div)
  glm_shannon <- glm(S.shannon ~ site + ocean_depth + site * ocean_depth,
                     family = gaussian(link = identity),
                     data = subset_data_div)
  glm_evenness <- glm(S.evenness ~ site + ocean_depth + site * ocean_depth,
                      family = gaussian(link = identity),
                      data = subset_data_div)
  glm_chao <- glm(S.chao1 ~ site + ocean_depth + site * ocean_depth,
                  family = gaussian(link = identity),
                  data = subset_data_div)
  
 return(list(glm_obs, glm_shannon, glm_evenness, glm_chao))

  }

glm_anelidos <- calculate_glm(anelidos_alphadiv)
glm_artro <- calculate_glm(artro_alphadiv)
glm_bryozoa <- calculate_glm(bryozoa_alphadiv)
glm_cnidaria <- calculate_glm(cnidaria_alphadiv)

summary(glm_anelidos[[1]])
summary(glm_anelidos[[2]])
summary(glm_anelidos[[3]])
summary(glm_anelidos[[4]])
summary(glm_artro[[1]])
summary(glm_artro[[2]])
summary(glm_artro[[3]])
summary(glm_artro[[4]])
summary(glm_bryozoa[[1]])
summary(glm_bryozoa[[2]])
summary(glm_bryozoa[[3]])
summary(glm_bryozoa[[4]])
summary(glm_cnidaria[[1]])
summary(glm_cnidaria[[2]])
summary(glm_cnidaria[[3]])
summary(glm_cnidaria[[4]])


nmds_bray_function <- function(ps_subset, tax_subset) {
  all_pca <- ordinate(
  physeq = ps_subset,
  method = "NMDS",
  distance = "bray")
  
  #pasar resultados a un data frame 
  all_pca_points <- as.data.frame(all_pca$points)
  df_subset <- data.frame(sample_data(ps_subset))
  
  #agregar variables
  all_pca_df <- all_pca_points  %>% 
    mutate(site = df_subset$site, 
           ocean_depth = df_subset$ocean_depth, 
           factor = paste(df_subset$site, df_subset$ocean_depth, sep = "_"))
  
  
  nmds_bray <- ggplot(data = all_pca_df,
                      aes(x = MDS1,
                          y = MDS2
                      )) + 
    stat_ellipse(type = "norm", linetype = 2,
                 show.legend = FALSE, size = 0.6, alpha = 0.9, aes(color = ocean_depth)) +
    geom_point(size = 5, aes(shape = site, color = ocean_depth)) +
    scale_shape_manual(values = c(15, 19), 
                       labels = c(`Algarrobo` = "Algarrobo", `Las_Cruces` = "Las Cruces")) +
    scale_color_manual(values = c("#17becf", "#e46385")) +
    labs(color = "Ocean Depth", shape = "Site", x = "nMDS1", y = "nMDS2", title = tax_subset) +
    theme_bw() 
  #test adonis
  dist_sfs_r <- phyloseq::distance(ps_subset, method = "bray")
  test <- adonis2(dist_sfs_r ~ (site + ocean_depth + site * ocean_depth),  data = df_subset, permutations = 1e3)
  test
  
  return(list (test = test, plot = nmds_bray))
  
}

nmds_anelidos <- nmds_bray_function(anelidos, "Annelida")
nmds_artro    <- nmds_bray_function(artro, "Arthropoda")
nmds_bryozoa  <- nmds_bray_function(briozoa, "Bryozoa")
nmds_cnidaria <- nmds_bray_function(cnidaria, "Cnidaria")

plot_nmds_subset <- nmds_anelidos$plot + nmds_artro$plot + nmds_cnidaria$plot +
  plot_layout(guides = 'collect', axis_titles = "collect") & 
  theme(
    plot.title = element_text(size = 12.5, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
    axis.title = element_text(size = 12.5),  # Tamaño de los títulos de los ejes
    legend.title = element_text(size = 13.5), # Tamaño del título de la leyenda
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(2, "lines"),
    axis.text.x = element_text(size = 7.5),   # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 7.5)) &
  guides(fill = guide_legend(override.aes = list(size = 6)))# Ajustar el tamaño de los símbolos
#ggsave("Plot_papper/plot_nmds_bray_subset.png", plot_nmds_subset, width = 8, height = 3)

##PA

anelidos_pa <- subset_taxa(r_ps_oc_pa, Phylum == "Annelida")
artro_pa    <- subset_taxa(r_ps_oc_pa, Phylum == "Arthropoda")
briozoa_pa  <- subset_taxa(r_ps_oc_pa, Phylum == "Bryozoa")
cnidaria_pa <- subset_taxa(r_ps_oc_pa, Phylum == "Cnidaria")

nmds_jaccard_function <- function(ps_subset_pa, tax_subset) {
  all_pca <- ordinate(
    physeq = ps_subset_pa,
    method = "NMDS",
    distance = "jaccard")
  
  #pasar resultados a un data frame 
  all_pca_points <- as.data.frame(all_pca$points)
  df_subset <- data.frame(sample_data(ps_subset_pa))
  
  #agregar variables
  all_pca_df <- all_pca_points  %>% 
    mutate(site = df_subset$site, 
           ocean_depth = df_subset$ocean_depth, 
           factor = paste(df_subset$site, df_subset$ocean_depth, sep = "_"))
  
  
  nmds_jaccard <- ggplot(data = all_pca_df,
                      aes(x = MDS1,
                          y = MDS2
                      )) + 
    stat_ellipse(type = "norm", linetype = 2,
                 show.legend = FALSE, size = 0.6, alpha = 0.9, aes(color = ocean_depth)) +
    geom_point(size = 5, aes(shape = site, color = ocean_depth)) +
    scale_shape_manual(values = c(15, 19), 
                       labels = c(`Algarrobo` = "Algarrobo", `Las_Cruces` = "Las Cruces")) +
    scale_color_manual(values = c("#17becf", "#e46385")) +
    labs(color = "Ocean Depth", shape = "Site", x = "nMDS1", y = "nMDS2", title = tax_subset) +
    theme_bw()
  
  #test adonis
  dist_sfs_r <- phyloseq::distance(ps_subset, method = "jaccard")
  test <- adonis2(dist_sfs_r ~ (site + ocean_depth + site * ocean_depth),  data = df_subset, permutations = 1e3)
  test
  
  return(list (test = test, plot = nmds_jaccard))
  
}


nmds_anelidos_pa <- nmds_bray_function(anelidos_pa, "Annelida")
nmds_artro_pa    <- nmds_bray_function(artro_pa, "Arthropoda")
nmds_bryozoa_pa  <- nmds_bray_function(briozoa_pa, "Bryozoa")
nmds_cnidaria_pa <- nmds_bray_function(cnidaria_pa, "Cnidaria")

plot_nmds_subset_pa <- nmds_anelidos_pa$plot + nmds_artro_pa$plot + nmds_cnidaria_pa$plot +
  plot_layout(guides = 'collect', axis_titles = "collect") & 
  theme(
    plot.title = element_text(size = 12.5, hjust = 0.5, vjust = 0.5),  # Tamaño del título general
    axis.title = element_text(size = 12.5),  # Tamaño de los títulos de los ejes
    legend.title = element_text(size = 13.5), # Tamaño del título de la leyenda
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(2, "lines"),
    axis.text.x = element_text(size = 7.5),   # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 7.5)) &
  guides(fill = guide_legend(override.aes = list(size = 6)))# Ajustar el tamaño de los símbolos
#ggsave("Plot_papper/plot_nmds_jaccard_subset.png", plot_nmds_subset_pa, width = 8, height = 3)

###---------------Beta diversity-------------
# Load necessary packages
library(phyloseq)
library(betapart)

# Assuming your phyloseq object is called `physeq`
physeq <- r_ps_oc # replace with your actual phyloseq object

# Merge samples by site
merged_physeq <- merge_samples(physeq, group = "factor")

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
turnover_jaccard <- beta_div_jaccard$beta.jtu
nestedness_jaccard <- beta_div_jaccard$beta.jne

# Print results for Jaccard
print(turnover_jaccard)
print(nestedness_jaccard)

