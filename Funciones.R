#Librerias
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(decontam)
library(vegan)
library(FSA)
library(ggVennDiagram)
library(patchwork)
library(stringr)
library(fitdistrplus)
library(extrafont)
library(car)
library(MASS)
library(betareg)
library(e1071)
#library(multtest)
#library(gridExtra)
#library(grid)
#library(rlang)
#library(ape)

# Importar y cargar fuentes
font_import()
loadfonts(device = "win")

###############FUNCTIONS################


#create a phyloseq compatible otu table format
#Args:
#>file_otu (in the rows the ASVs and in the columns the samples)
#
#Returns:
#>OTU table
otu <- function(file_otu) {
  otus <- as.matrix(file_otu)
  otus_format <- otu_table(otus, taxa_are_rows = TRUE)
  return(otus_format)
}


#create a phyloseq compatible taxonomy table format
#Args:
#>file_taxa (in the rows the ASVs and in the columns the taxonomic levels)
#
#Returns:
#>Taxonomy table
taxa <- function(file_taxa) {
  taxas <- (as.matrix(file_taxa))
  taxa_format <- (tax_table(taxas))
  return(taxa_format)
}


#create a phyloseq compatible samplee data format
#Args:
#>file_data (in the rows the samples and in the columns the variables)
#
#Returns:
#>Sample table (sample information)
data <- function(file_df) {
  df <- data.frame(file_df)
  df_format <- sample_data(df)
  return(df_format)
}


#agrupar por tax and clean tax con reads = 0
#Args:
#>object phyloseq (pa or not)
#>taxonoic level (remember that this has the same name as the taxonomy table)
#default
#>NArm = FALSE
#
#Returns:
#>object phyloseq whit reads > 0 in the specified taxonomic level
clean_zero_reads <- function(ps, taxonomic_level) {
  tax_glom_ps <- tax_glom(ps, taxonomic_level, NArm = FALSE) #agrupa por taxonomic leves eg. genus #nolint
  clean_zero_ps <- prune_taxa(taxa_sums(tax_glom_ps) > 0,
                              tax_glom_ps) #remove taxonommic_level=0
  return(clean_zero_ps)
}


#rarefaction
#Args:
#>phyloseq object without reads=0 (reuturn fnction clean_zero_reads)
#Default: sample size is the min number reads in the otu table
#
#Returns:
#>Object phyloseq rarefy
rarefaction <- function(clean_zero_ps) {
  sample_sums <- sample_sums(clean_zero_ps) #nolint
  min_reads_sample <- which.min(sample_sums) #nolint
  min_reads <- sample_sums[min_reads_sample] #nolint
  min_reads <- as.numeric(min_reads) #nolint
  rarefied_ps <- rarefy_even_depth(clean_zero_ps,
                                   sample.size = min_reads,
                                   replace = TRUE)
  return(rarefied_ps)
}


#curves rarefaction, data, alpha div and plot
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require("plyr") # ldply
  require("reshape2") # melt

  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)

    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)

    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)

    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity),
                                   varnames = c('Sample', 'Measure'),
                                   value.name = 'Alpha_diversity')

    molten_alpha_diversity
  }

  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))

  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

  rarefaction_curve_data
}


#summary of the result obtein in the rarefaction curves
#Args:
#>phyloseq object rarefy (can be unrarefactioned, review other functions)
#>rarefaction curve (retun function rarefaction_curve)
#
#Returns:
#>data frame whit depth, sample, and measue (alpha or observed)
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


#' Calculate relative abundance of sequences by taxonomic level
#'
#' @param ps A phyloseq object (with or without rarefaction, without pa)
#' @param taxonomic_level The taxonomic level at which to calculate relative abundance
#'
#' @return A phyloseq object agglomerated by taxonomic level and transformed to relative abundance
#'
#' @examples
#' relative_ab(r_ps, taxonomic_level)
#'
relative_ab <- function(ps, taxonomic_level) {
  lev_tax_abundance <- ps %>%
    tax_glom(taxrank = taxonomic_level) %>%
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>%                         # Melt to long format
    arrange(!!sym(taxonomic_level))

  return(lev_tax_abundance)
}


#' Calculate absolute abundance of sequences by taxonomic level
#'
#' @param ps A phyloseq object (with or without rarefaction, without pa)
#' @param taxonomic_level The taxonomic level at which to calculate relative abundance
#'
#' @return A phyloseq object agglomerated by taxonomic level and transformed to relative abundance
#'
#' @examples
#' relative_ab(r_ps, taxonomic_level)
#'
absolute_ab <- function(ps, taxonomic_level) {
  lev_tax_abundance <- ps %>%
    tax_glom(taxrank = taxonomic_level) %>%
    psmelt() %>%                         # Melt to long format
    arrange(!!sym(taxonomic_level))

  return(lev_tax_abundance)
}


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

#clean NA y ""
#Args:
#>obect phyloseq rarery or not pa or not
#>taxonomic_level
#
#Return:
#>phyloseq object without NA and ""
clean_NA <- function(ps_pa, taxonomic_level) {
  ps_pa_level_clean <- subset_taxa(ps_pa, taxonomic_level!="NA") #nolint
  ps_pa_level_clean <- subset_taxa(ps_pa_level_clean, taxonomic_level!="unassigned") #nolint
  return(ps_pa_level_clean)
}


#contar numero de asv vacias ("") en base a tax table
#Arg:
#>phylosec object (rarefy or not)
#
#return:
#>numero de filas en taxa tabla(numero de asv)
nrow_subset <- function(ps_subset) {
  nrow_subset <- ps_subset %>%
  tax_table() %>%                           # Obtener la tabla de taxonomía
  as.data.frame() %>%                       # Convertir a data frame
  filter(Specie == "") %>%                 # Filtrar por especies vacías
  nrow()

  return(nrow_subset)
}


