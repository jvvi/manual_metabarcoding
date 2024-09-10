#directory
setwd("C:/Users/javie/Documents/Scripts/agua_arms/agua_arms")

project_path <- "C:/Users/javie/Documents/Scripts/agua_arms/agua_arms"
project_path

#open files
asv_1 <- read.table("asv_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1 <- read.table("tax_table.txt", header = TRUE, sep = "\t", row.names = 1)
tax_1[tax_1 == ""] <- "Unassigned"
df <- read.table("df_table.txt", header = TRUE, sep = "\t", row.names = 1)

df <- df %>% 
mutate(factor =  paste(site, ocean_depth, sep = "_")) %>%
mutate(factor_r =  paste(site, ocean_depth, replica, sep = "_"))

##table insect###
#load("insect/longDF_table.RData")
#tax_2 <- aggregate(longDF[3:12], longDF["taxID"], head, 1)
#tax_2<- tax_2 %>% 
#  dplyr::select(1,5:11) %>%
#  rename(CO1_seq_number = taxID, Kingdom = kingdom, Phylum = phylum, Class= class, Order= order, 
#                        Family = family, Genus = genus, Specie = species)
#tax_2 <- tax_2 %>%
#  mutate_all(~ ifelse(. %in% c("Metazoa", "Viridiplantae"), "Eukarya", .))

#asv_2 <- aggregate(longDF[13:ncol(longDF)], longDF["taxID"], sum)
#colnames(asv_2) <- gsub("-", ".", colnames(asv_2))
#asv_2 <- rename(asv_2, CO1_seq_number = taxID)
#
##merge 
#setdiff(colnames(asv_1), colnames(asv_2))
#intersect(colnames(tax_1), colnames(tax_2))
#asv <- merge(x = asv_1, y = asv_2, all = TRUE)
#tax <- merge(x = tax_1, y = tax_2, all = TRUE)
#rownames(asv) <- asv$CO1_seq_number; asv <- asv[,-1]
#rownames(tax) <- tax$CO1_seq_number; tax <- tax[,-1]

#####create object phylseq#####
otu_ps <- otu(asv_1)
tax_ps <- taxa(tax_1)
data_ps <- data(df)
ps <- phyloseq(otu_ps, tax_ps, data_ps)
ps <- subset_samples(ps, study == "arms")


####remove contaminarion based on the negative control####
#data preparation
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#visualization
contam <- ggplot(data = df, aes(x = sample_name, y = LibrarySize,
        color = factor)) + geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 6))
contam
#ggsave("plot_papper/contam.png", contam, width = 6, height = 3)

#identification of negative controls
sample_data(ps)$is.neg <- sample_data(ps)$type == "control"

#detection contaminants
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
filter(contamdf.prev05, contaminant == TRUE)

#selection seq contaminants
contam <- head(which(contamdf.prev05$contaminant))
contam <- rownames(contamdf.prev05[c(contam),])

#remove contaminants
otu_table(ps) <- otu_table(ps)[!(row.names(otu_table(ps)) %in% contam),]
ps <- subset_samples(ps, type != "control")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)


#contar numero de secuencias
ps <- prune_taxa(taxa_sums(ps) > 0,
                 ps)

count_total <- ntaxa(ps)

count_know <- ps %>%
  subset_taxa(Phylum != "Unassigned")  %>% 
  subset_taxa(Phylum != "NA")  %>% 
  ntaxa()

100 - ((count_know*100)/count_total)


#### data clean####
#limpiar datos de sp no marinas, y peces
ps = subset_taxa(ps, Specie !="Eurytemora foveola")#AGUA DULCE
ps = subset_taxa(ps, Specie !="Eurytemora herdmani")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Fridericia")
ps = subset_taxa(ps, Genus !="Dussartcyclops")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Neoergasilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Sinocalanus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Bryodrilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Rhysida")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Potamothrix")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Baikalodrilus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Fridericia")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Lumbriculus")#AGUA DULCE
ps = subset_taxa(ps, Genus !="Lamprodrilus")#AGUA DULCE
ps = subset_taxa(ps, Family !="Diaptomidae")
ps = subset_taxa(ps, Family !="Paramelitidae")
ps = subset_taxa(ps, Family !="Potamidae")
ps = subset_taxa(ps, Family !="Ampullariidae")
ps = subset_taxa(ps, Family !="Orobdellidae")
ps = subset_taxa(ps, Family !="Trombiculidae")
ps = subset_taxa(ps, Family !="Buthidae")
ps = subset_taxa(ps, Family !="Phytoseiidae")
ps = subset_taxa(ps, Family !="Unionicolidae") 
ps = subset_taxa(ps, Order !="Crassiclitellata") #gusano terrestre
ps = subset_taxa(ps, Order !="Lumbriculida")#AGUA DULCE
ps = subset_taxa(ps, Order !="Anura")
ps = subset_taxa(ps, Order !="Artiodactyla") #angulados
ps = subset_taxa(ps, Order !="Chiroptera") #murcielago
ps = subset_taxa(ps, Order !="Eulipotyphla")
ps = subset_taxa(ps, Order !="Physariida")
ps = subset_taxa(ps, Order !="Mesostigmata")
ps = subset_taxa(ps, Order !="Araneae")
ps = subset_taxa(ps, Order !="Sarcoptiformes")
ps = subset_taxa(ps, Order !="Ixodida")
ps = subset_taxa(ps, Order !="Opiliones")
ps = subset_taxa(ps, Order !="Squamata")
ps = subset_taxa(ps, Class !="Eumycetozoa") #nolint
ps = subset_taxa(ps, Class !="Insecta") #nolint
ps = subset_taxa(ps, Class !="Collembola")
ps = subset_taxa(ps, Class !="Diplopoda")
ps = subset_taxa(ps, Class !="Actinopteri")
ps = subset_taxa(ps, Phylum !="Ascomycota")
ps = subset_taxa(ps, Phylum !="Basidiomycota")
ps = subset_taxa(ps, Phylum !="Mucoromycota")
ps = subset_taxa(ps, Phylum !="Chordata")
ps = subset_taxa(ps, Phylum !="Chlorophyta")
ps = subset_taxa(ps, Phylum !="Rhodophyta")
ps = subset_taxa(ps, Phylum !="Discosea")
ps = subset_taxa(ps, Phylum !="Bacillariophyta")
ps = subset_taxa(ps, Kingdom !="Bacteria")


#contar numero de secuencias
ps <- prune_taxa(taxa_sums(ps) > 0,
                      ps)

count_total <- ntaxa(ps)

count_know <- ps %>%
  subset_taxa(Phylum != "Unassigned")  %>% 
  subset_taxa(Phylum != "NA")  %>% 
  ntaxa()

100 - ((count_know*100)/count_total)

#aglomerar por specie, contar nuemro de ASV y nuemeor de NA
ps <- clean_zero_reads(ps, "Specie")
