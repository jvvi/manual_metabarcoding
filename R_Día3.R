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
D1 <- ggplot(sfs_alphadiv, aes(x = ocean_depth, y = S.obs, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Richness", x =NULL, fill = "Depth",  tag = "A") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  theme(legend.position = "none") +
  theme_bw()
D2 <- ggplot(sfs_alphadiv, aes(x=ocean_depth, y=S.chao1, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Chao1", x =NULL, fill = "Depth",  tag = "B") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  theme(legend.position = "none") +
  theme_bw()
D3 <- ggplot(sfs_alphadiv, aes(x=ocean_depth, y=S.evenness, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Evenness", x ="Site", fill = "Depth",  tag = "C") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  theme(legend.position = "none") +
  theme_bw()
D4 <- ggplot(sfs_alphadiv, aes(x=ocean_depth, y=S.shannon, fill =ocean_depth)) +
  geom_boxplot() +
  labs(title = NULL, y = "Shannon", x ="Site", fill = "Depth", tag = "D") +
  scale_fill_manual(
    values = c("30m" = "#17becf", "60m" = "#e46385")) +
  theme(legend.position = "none") +
  theme_bw()

plot_div <-  D1 + D2 + D3 + D4 +
  plot_layout(axis_titles = "collect") & theme(
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

#ELEGIR PS A OCUPAR-------------------
r_ps_oc <- r_ps_sfs

#PRESENCIA Y AUSENCIA-----------
#presence absence of the objeto phyloseq
#Args:
#>phyloseq objectrarefy or not
#
#Return:
#>phyloseq object whit presence and absence
normalize_pa_ps <- function(ps) {
  comm <- as(phyloseq::otu_table(ps), Class = "matrix")
  comm <- vegan::decostand(comm, method = "pa", MARGIN = 1)
  ps_pa <- ps
  trows <- phyloseq::taxa_are_rows(ps)
  phyloseq::otu_table(ps_pa) <- phyloseq::otu_table(comm, taxa_are_rows = trows)
  return(ps_pa)
}

r_ps_oc_pa <- normalize_pa_ps(r_ps_oc)

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
  geom_point(size = 5, aes(color = ocean_depth)) +
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
new_labels <- c("AL30", "AL60")
 
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

taxonomic_level <- "Specie"

r_ps_AL30 <- clean_zero_reads(r_ps_AL30, taxonomic_level)
r_ps_AL60 <- clean_zero_reads(r_ps_AL60, taxonomic_level)


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

tax_table_AL30_c <- as.character(tax_table_AL30$all_taxa_level)
tax_table_AL60_c <- as.character(tax_table_AL60$all_taxa_level)

asv_beta <- list(
  AL60 = tax_table_AL60_c, AL30 = tax_table_AL30_c)


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
