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
