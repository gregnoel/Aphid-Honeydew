##Mise en place du set directory tjs utile lorsqu'on veut importer des jeux de donn?es 
setwd(dir = "F://Analyse_Microbio/Honeydew/01032023")
R.version
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MicrobiotaProcess")
library(BiocManager)
library(phyloseq)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("apeglm")

library(DESeq2)
library(ape)
library(ggplot2)
library(vegan)
library(metacoder)
###Pr?paration des graphiques par rapport aux OTUs en barplot
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(microbiome)
library(magrittr) # need to run every time you start R and want to use %>%
library(devtools)
library(picante)
library(phyloseq)

#############################################################
#Preparation des donn?es et construction d'un objet phyloseq#
#############################################################

#Importation de la matrice de communaute bacterinne / sample
OTUTable<-read.csv2("F://Analyse_Microbio/Honeydew/01032023/feature-table_Honeydew01032023.csv",row.names=1)
str(OTUTable)
#Firstcolumn<-OTUTable[,1]
OTUTable[,1:78]<-sapply(OTUTable[,1:78], function(x) as.numeric(x))

#OTUTable<-as.data.frame(sapply(OTUTable[,2:85], function(x) as.numeric(x)), make.names = T)
#OTUTable<-cbind(Firstcolumn,OTUTable)


OTUTable2<-as.matrix(OTUTable)
OTU = otu_table(OTUTable2, taxa_are_rows = TRUE)

#str(OTUTable2)
OTU_taxonomy<-read.csv2("F://Analyse_Microbio/Honeydew/01032023/taxonomy_honeydew01032023_2.csv", row.names=1,na.strings = "")

TAX = phyloseq::tax_table(as.matrix(OTU_taxonomy))
#Ici on peut direct cr?er un objet phyloseq en se basant juste sur la table d'OTU et de taxo
physeq = phyloseq(OTU, TAX)
#plot_bar(physeq, fill = "Phylum",facet_grid = Domain)
#?plot_bar

#On peut ajouter l'arbre et les donn?es envi dans l'objet phyloseq (optionnel)
EnviData<-read.csv2("F://Analyse_Microbio/Honeydew/Metadata_Honeydew2.csv", row.names=1)
SAM = sample_data(EnviData)

#Changement des labels dans l'arbre phylo pour qu'elle correspondent avec le reste du jeu de donn?es
#Change.label<-read.csv2("D://Analyse_Microbio/Honeydew/Metadata_Honeydew2.csv")
#Change.label<-as.character(Change.label$Taxa)


OTU_tree<-phyloseq::read_tree("F://Analyse_Microbio/Honeydew/01032023/rootedtree_Honeydew.nwk")
#OTU_tree$tip.label<-Change.label


#OTU_tree <- compute.brlen(OTU_tree, method = "Grafen")#Mani?re de recalculer l'arbre phylo faut checker l'utilit? de cette commande dans la documentation

#Fusion de tout les fichiers pour cr?er l'objet phyloseq selon deux m?thodes
physeq1 = merge_phyloseq(physeq, SAM, OTU_tree)
physeq1 = merge_phyloseq(physeq, SAM)
physeq1
#physeq2 = phyloseq(OTU, TAX, SAM, OTU_tree)
#physeq2

#Verification que les deux objets phyloseq sont identiques, si TRUE ?a veut dire que c'est bon on peut passer ? la suite
#identical(physeq1, physeq2)


#Visualization of the data characteristics

sample_names(physeq1)
rank_names(physeq1)
sample_variables(physeq1)
#######################################################
###Analyse sur les donn?es des ?chantillons du miellat##
#######################################################

plot_bar(physeq1, fill="Phylum")

to_remove <- c("P1A07", "P1B09", "P1C04BLANK","P1C04BLANK_01","P1E01","P2C05WATER","P2F04","P2H12","P2F02")

physeq1_filtered <- prune_samples(!(sample_names(physeq1) %in% to_remove), physeq1)
physeq1_filtered
##Bon ici on voit d?j? que dans certains ECH P1A07, P1B09, P1C04BLANK, P1E01, P2C05WATER, P2F04, P2H12, y a a pas grand chose ==> on les vire, on va prendre un seuil ? >1000 reads, 
#physeq_filter <- phyloseq::genefilter_sample(physeq1, filterfun_sample(function(x) x >= 5))
#physeq1_filtered <- prune_taxa(physeq_filter, physeq1)
#plot_bar(GPfr , fill="Phylum")
plot_bar(physeq1_filtered, fill="Phylum")
#Int?ressant de voir la somme des comptage des identifications des bact?ries ? travers le NGS
# Create table, number of features for each phyla to see distribution and NA
#table(tax_table(GPfr_filtered)[, "Phylum"], exclude = NULL)



################################################################
#analyse en multivari?e des ?chantillons par esp?ces de pucerons
################################################################
#On ne prends que les taxons qui sont au dessus de 1000, v?rifier ici quelle est le seuil de filtrage par taxon
GPr  = transform_sample_counts(physeq1, function(x) x / sum(x) )
GPfr = phyloseq::filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)


##Par esp?ce
# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ord <- ordinate(physeq1, "MDS", "bray")
plot_ordination(physeq1, ord, color = "Species") +
  geom_point(size = 5)



##Par viellissement de miellat
# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ord <- ordinate(physeq1, "MDS", "bray")
plot_ordination(physeq1, ord, color = "Ageing", shape = "Ageing") +
  geom_point(size = 5)



####Sortir les ?chantillons par viellissement de miellat
tabletopaper<-as.data.frame(table(meta(physeq1_filtered)$Species, meta(physeq1_filtered)$Ageing))
write.csv2(tabletopaper,"table_to_article.csv")
physeq1_filtered_pisum <- subset_samples(physeq1_filtered, Species == "A_pisum")
physeq1_filtered_fabae <- subset_samples(physeq1_filtered, Species == "A_fabae")
print(physeq1_filtered_pisum)
table(meta(physeq1_filtered_pisum)$Species, meta(physeq1_filtered_pisum)$Ageing)

GPr  = transform_sample_counts(physeq1_filtered_fabae, function(x) x / sum(x) )
GPfr = phyloseq::filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
write.csv2(GPfr@tax_table,"Preliminary_results/taxo_fabae.csv")



#####################
##Venn Diagram#######
#####################
library(MicrobiotaProcess)
library(VennDiagram)
library(UpSetR)

##Per Aphid species
vennL<-get_vennlist(physeq1_filtered, factorNames="Species")
names(vennL)<-c("Aphis fabae","Acyrthosiphon pisum")
vennp <- venn.diagram(vennL,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#46495c", "#14a360"),
                      cat.col=c("#46495c", "#14a360"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
                    

veen_final<-grid::grid.draw(vennp)

ggsave("Venn_Species.png",veen_final,dpi = 500)
dev.off()


#per ageing
vennL<-get_vennlist(physeq1_filtered, factorNames="Ageing")

vennp1 <- venn.diagram(vennL,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#CBD588", "#5F7FC7", "orange","#DA5724"),
                      cat.col=c("#CBD588", "#5F7FC7", "orange","#DA5724"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")

grid::grid.draw(vennp1)

vennL<-get_vennlist(physeq1_filtered, factorNames="Ageing_Species")

vennp1 <- venn.diagram(vennL,
                       height=8,
                       width=8, 
                       filename=NULL, 
                       fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                              "#AD6F3B", "#673770"),
                       cat.col=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                 "#AD6F3B", "#673770"),
                       alpha = 0.85, 
                       fontfamily = "serif",
                       fontface = "bold",
                       cex = 1.2,
                       cat.cex = 1.3,
                       cat.default.pos = "outer",
                       cat.dist=0.1,
                       margin = 0.1, 
                       lwd = 3,
                       lty ='dotted',
                       imagetype = "svg")

grid::grid.draw(vennp1)

###########################
############MPSE###########
###########################

library(ggplot2)
library(phyloseq)
library(shadowtext)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(MicrobiotaProcess)

#Creation of MPSE object
MPSE_Test_Microbio<-as.MPSE(physeq1_filtered)
#MPSE_Test_Microbio@otutree@phylo$tip.label


# Rarefied species richness
MPSE_Test_Microbio %<>% mp_rrarefy()
# 'chunks' represent the split number of each sample to calculate alpha
# diversity, default is 400. e.g. If a sample has total 40000
# reads, if chunks is 400, it will be split to 100 sub-samples
# (100, 200, 300,..., 40000), then alpha diversity index was
# calculated based on the sub-samples. 
# '.abundance' the column name of abundance, if the '.abundance' is not be 
# rarefied calculate rarecurve, user can specific 'force=TRUE'.
MPSE_Test_Microbio %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 100
  )
# The RareAbundanceRarecurve column will be added the colData slot 
# automatically (default action="add")
MPSE_Test_Microbio %>% print(width=180)

p3 <- MPSE_Test_Microbio %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Species, 
    plot.group = TRUE
  ) + scale_x_discrete(labels=c('Aphis fabae', 'Acyrtosyphon pisum')) +
 scale_color_manual(values=c("#46495c", "#14a360")) +
  scale_fill_manual(values=c("#46495c", "#14a360"))  
p3<-p3 + theme_bw(base_size = 18)
p3 <- p3 +  theme(legend.position = "none") 
p3
#ggsave("RarefactionCurve05122024.png", plot = p3, dpi=1500)


MPSE_Test_Microbio %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
MPSE_Test_Microbio
library(gghalves)

f1 <- MPSE_Test_Microbio %>% 
  mp_plot_alpha(
    .group=Species, 
    .alpha=c(Observe, Chao1, Shannon)
  ) +
  scale_fill_manual(values=c("#46495c", "#14a360")) +
  scale_color_manual(values=c("#46495c", "#14a360"))

f1<-f1 + scale_x_discrete(labels=c('Aphis fabae', 'Acyrthosyphon pisum')) +  theme(legend.position = "none")
#graph


ggsave("AlphaDiversityMetrics05122024.png", plot = f1, dpi=1500)


MPSE_Test_Microbio %<>% 
  mp_cal_alpha(.abundance=Abundance, force = T)
MPSE_Test_Microbio

###Per species and ageing
f1 <- MPSE_Test_Microbio %>% 
  mp_plot_alpha(
    .group=c(Ageing,Species), 
    .alpha=c(Observe, Chao1, Shannon),
    margin_top = 0.5,
    step_increase =0.2,
    tip_length = 0.03,
    size = 0.5,
    textsize = 2.5
  ) 
f1

#graph
ggsave("AlphaDiversityMetrics_Ageing_15122023.png", plot = f1, dpi=500)


###ABundance
MPSE_Test_Microbio %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Species
  )
MPSE_Test_Microbio@assays

library(ggalluvial)
library(ggh4x)
# visualize the relative abundance of top 20 phyla for each sample.
p1 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class = Phylum, 
    topn = 5,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p1<-p1+
  theme(legend.position = "right") + ggtitle("")+
                     labels=c("A. fabae", "A. pisum"))+
 
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p1
p1<-p1 +  scale_fill_manual(values=c("darkolivegreen3","cornsilk2","gold", "#00CDCD","#FF8247", "#4F94CD"),
                          labels=c("Others", "Plantomycetota","Bacteriodota", "Verrucomicrobiota","Firmicutes","Proteobacteria"))

p1 
ggsave("Phylum_Level11122023.png", plot = p1, dpi=500)

#Write the proportion !
write.csv2(ggplot_build(p1)$data[[1]],"data_phylum.csv")

#Family
p2 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class = Family, 
    topn = 6,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p2<-p2+
  theme(legend.position = "right") + ggtitle("")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A. fabae", "A. pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p2
p2<-p2 +  scale_fill_manual(values=c("darkolivegreen3","cornsilk2","gold", "#00CDCD","#FF8247", "#4F94CD","firebrick1"),
                               labels=c("Others", "Erwiniaceae","Verrucomicrobiaceae", "Moraxellaceae","Pseudomonadaceae","Staphylococcaceae", "Morganellaceae"))

ggsave("Family_Level11122023.png", plot = p2, dpi=500)
write.csv2(ggplot_build(p2)$data[[1]],"data_family.csv")


p3 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class =Genus, 
    topn = 6,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p3<-p3+
  theme(legend.position = "right") + ggtitle("")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A. fabae", "A. pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p3<-p3 +  scale_fill_manual(values=c("darkolivegreen3","cornsilk2","gold", "#00CDCD","#FF8247", "#4F94CD","firebrick1"),
                        labels=c("Others", "Erwinia","Prosthecobacter", "Acinetobacter","Pseudomonas","Staphylococcus", "Buchnera"))
p3
ggsave("Genus_Level11122023.png", plot = p3, dpi=500)
#write.csv2(ggplot_build(p3)$data[[1]],"data_genus.csv")



p4<-p4+
  theme(legend.position = "right") + ggtitle("")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A. fabae", "A. pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p3<-p3 +  scale_fill_manual(values=c("darkolivegreen3","cornsilk2","gold", "#00CDCD","#FF8247", "#4F94CD","firebrick1"),
                            labels=c("Others", "Erwinia","Prosthecobacter", "Acinetobacter","Pseudomonas","Staphylococcus", "Buchnera"))
p3


###With ageing
library(ggfittext)
p4 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=c(Ageing,Species), 
    taxa.class =Genus, 
    topn = 6,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p4<-p4 + scale_x_discrete(labels = c("24h" = "Fresh", "48h" = "24h", "72h" = "48h", "96h" = "72h")) 

write.csv2(ggplot_build(p4)$data[[1]],"Ageing_genus.csv")


#p4 + geom_fit_text(aes(label = paste0(round(RareAbundance,1), "%")), position=position_stack(vjust=.5), show.legend=F, color='white')


p4 <- p4 +  scale_fill_manual(values=c("#54278F","#557982","#B76242", "#4D4D4D","#998EC0", "#D6988A","#1F78B4"),
                              labels=c("Other genera", "Erwinia spp.","Prosthecobacter spp.", "Acinetobacter spp.","Pseudomonas spp.","Staphylococcus spp.", "Buchnera spp."))
p4



p5 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=c(Ageing,Species), 
    taxa.class =OTU, 
    topn = 6,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p5
p5 <- p5 +  scale_fill_manual(values=c("#54278F","#557982","#B76242", "#4D4D4D","#998EC0", "#D6988A","#1F78B4"),

                                                            labels=c("Other ASVs", "Acinetobacter soli","Pseudomonas", "Prosthecobacter","Buchnera","Staphylococcus sciuri", "Buchnera")) + scale_x_discrete(labels = c("24h" = "Fresh", "48h" = "24h", "72h" = "48h", "96h" = "72h")) 
p5
p4
p6 <- p4/p5
p6

#ggsave("Ageing_prop05122024.svg", plot = p6, width = 16, height = 24, dpi = 1500, units = "cm")
ggsave("Ageing_prop05122024_Genus.svg", plot = p4, dpi=1500)
ggsave("Ageing_prop05122024_ASV.svg", plot = p5, dpi=1500)

##Heatmap
h1 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Species,
    taxa.class = Family,
    relative = TRUE,
    topn = 10,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h1


h2 <- MPSE_Test_Microbio %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Species,
    taxa.class = Phylum,
    relative = FALSE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
#aplot::plot_list(gglist=list(h1, h2), tag_levels="A")



#####Beta Diversity

MPSE_Test_Microbio %<>% 
  mp_decostand(.abundance=Abundance)
MPSE_Test_Microbio


library(corrr)
# calculate the distance between the samples.
# the distance will be generated a nested tibble and added to the
# colData slot.
MPSE_Test_Microbio %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
MPSE_Test_Microbio

# mp_plot_dist provides there methods to visualize the distance between the samples or groups
# when .group is not provided, the dot heatmap plot will be return
p1 <- MPSE_Test_Microbio %>% mp_plot_dist(.distmethod = bray)
p1


# when .group is provided, the dot heatmap plot with group information will be return.
p2 <-  MPSE_Test_Microbio %>% mp_plot_dist(.distmethod = bray, .group = Species)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("orange", "deepskyblue"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = Species # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )

# when .group is provided and group.test is TRUE, the comparison of different groups will be returned
p3 <-  MPSE_Test_Microbio %>% mp_plot_dist(.distmethod = bray, .group = Species, group.test=TRUE, textsize=2)
p3 

####The PCOA ANalysis#####

MPSE_Test_Microbio %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
MPSE_Test_Microbio

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
MPSE_Test_Microbio %<>%
  mp_adonis(.abundance=hellinger, .formula=~Species, distmethod="bray", permutations=9999, action="add")
MPSE_Test_Microbio %>% mp_extract_internal_attr(name=adonis)
#ado_test$aov.tab


# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
library(ggside)
p2 <- MPSE_Test_Microbio  %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Species, 
    .color = Species, 
    .size = Observe, 
    show.side = F,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#46495c", "#14a360"), labels = c("Aphis fabae", "Acyrthosiphon pisum"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=9))
  ) +
  scale_color_manual(
    values=c("#46495c", "#14a360"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=9))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=9))
  )
p2

ggsave("Betadiversity.png", plot =p2, dpi=500)

####Hierarchical cluster analysis

MPSE_Test_Microbio  %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )
MPSE_Test_Microbio 

# if action = 'add', the result of hierarchical cluster will be added to the MPSE object
# mp_extract_internal_attr can extract it. It is a treedata object, so it can be visualized
# by ggtree.
sample.clust <- MPSE_Test_Microbio  %>% mp_extract_internal_attr(name='SampleClust')
sample.clust


library(ggtree)
p <- ggtree(sample.clust) + 
  geom_tippoint(aes(color=Species)) +
  geom_tiplab(as_ylab = TRUE) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))
p


library(ggtreeExtra)
library(ggplot2)
phyla.tb <- MPSE_Test_Microbio %>% 
  mp_extract_abundance(taxa.class=Phylum, topn=5)
# The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
phyla.tb



p1 <- p + 
  geom_fruit(
    data=phyla.tb,
    geom=geom_col,
    mapping = aes(x = RelRareAbundanceBySample, 
                  y = Sample, 
                  fill = Phyla
    ),
    orientation = "y",
    #offset = 0.4,
    pwidth = 3, 
    axis.params = list(axis = "x", 
                       title = "The relative abundance of phyla (%)",
                       title.size = 4,
                       text.size = 2, 
                       vjust = 1),
    grid.params = list()
  )
p1


#########BIOMARKER Discovery#######

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)
MPSE_Test_Microbio %>% print(width=150)

MPSE_Test_Microbio %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Species,
    first.test.alpha = 0.05
  )
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- MPSE_Test_Microbio %>% 
  mp_extract_tree(type="taxatree")
taxa.tree

# And the result tibble of different analysis can also be extracted 
# with tidytree (>=0.3.5)
taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Species, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))

##Add Tip Change

tips<-read.csv2("TpLAbel_tsfo.csv")
taxa.tree@phylo$tip.label<-tips$TipLabel


##extract the taxonomy information
MPSE_Test_Microbio %>% mp_extract_taxonomy() 

# Since taxa.tree is treedata object, it can be visualized by ggtree and ggtreeExtra
p1 <- ggtree(
  taxa.tree,
  layout = "circular",
  size = 0.3
  ) +
  geom_point(
    data = td_filter(isTip),
    fill="white",
    size=1,
    shape=21
  )
# display the high light of phylum clade.
p1

p2 <- p1 +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  )
p2
# display the relative abundance of features(OTU)

p3 <- p2 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Species, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Species,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#46495c","#14a360"))
p3
# display the tip labels of taxa tree
p4 <- p3 + geom_tiplab(size=2, offset=7.2)
p4
# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_Species,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )
p5

# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_Species)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_Species,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#46495c","#14a360"))

p7<- p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)

p7

MPSE_Test_Microbio@otutree@phylo$tip.label<-tips$TipLabel
MPSE_Test_Microbio@taxatree@phylo$tip.label<-tips$TipLabel

p <- MPSE_Test_Microbio %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.1
  ) +
  scale_fill_manual(values=c("#46495c", "#14a360")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#46495c", "#14a360")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999","#AD6F3B", "#673770","#D14285"
    )
  )
p + 
  geom_tiplab(size=2, offset=7.2)+ 
  theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)



f <- MPSE_Test_Microbio %>%
  mp_plot_diff_cladogram(
    label.size = 2.5,
    hilight.alpha = .3,
    bg.tree.size = .5,
    bg.point.size = 2,
    bg.point.stroke = .25
  ) +
  scale_fill_diff_cladogram( # set the color of different group.
    values = c("#46495c", "#14a360")
  ) +
  scale_size_continuous(range = c(1, 4))
f
ggsave("cladogram.png", plot = f, dpi=500)

f.box <- MPSE_Test_Microbio %>%
  mp_plot_diff_boxplot(
    .group = Species,
  ) %>%
  set_diff_boxplot_color(
    values = c("#46495c", "#14a360"),
    guide = guide_legend(title=NULL)
  )
f.box
f.bar <-  MPSE_Test_Microbio %>%
  mp_plot_diff_boxplot(
    taxa.class = c(Genus, OTU), # select the taxonomy level to display
    group.abun = TRUE, # display the mean abundance of each group
    removeUnknown = TRUE, # whether mask the unknown taxa.
  ) %>%
  set_diff_boxplot_color(
    values = c("#46495c", "#14a360"),
    guide = guide_legend(title=NULL)
  )
f.bar 
class(f.bar)



#ff<-aplot::plot_list(f.box, f.bar)
ggsave("bar_plotto.svg", plot = f.bar, dpi=500)

f.mahattan <- MPSE_Test_Microbio %>%
  mp_plot_diff_manhattan(
    .group = Sign_Species,
    .y = fdr,
    .size = 2.4,
    taxa.class = c('OTU', 'Genus'),
    anno.taxa.class = Phylum
  )
f.mahattan
ggsave("f_mahattan.png", plot = f.mahattan, dpi=500)


########


install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)


########24h VERSUUUUUUS#######

library(dplyr)
#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)
library(microViz)
physeq1_fabae <- microViz::ps_filter(physeq1_filtered, Species == "A_fabae")

MPSE_Test_Microbio_fabae<-as.MPSE(physeq1_fabae)
MPSE_Test_Microbio_fabae

# Rarefied species richness
MPSE_Test_Microbio_fabae %<>% mp_rrarefy()
# 'chunks' represent the split number of each sample to calculate alpha
# diversity, default is 400. e.g. If a sample has total 40000
# reads, if chunks is 400, it will be split to 100 sub-samples
# (100, 200, 300,..., 40000), then alpha diversity index was
# calculated based on the sub-samples. 
# '.abundance' the column name of abundance, if the '.abundance' is not be 
# rarefied calculate rarecurve, user can specific 'force=TRUE'.
MPSE_Test_Microbio_fabae%<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 100
  )
# The RareAbundanceRarecurve column will be added the colData slot 
# automatically (default action="add")
MPSE_Test_Microbio_fabae %>% print(width=180)

p3 <- MPSE_Test_Microbio_fabae %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Ageing, 
    plot.group = TRUE
  ) +
  scale_color_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019")) +
  scale_fill_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019"))  
p3<-p3 + ggtitle ("A. fabae - Ageing")+ theme_bw(base_size = 18)
p3
ggsave("Curve_Acc_fabae10032023.png", plot = p3, dpi=500)

MPSE_Test_Microbio_fabae %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
MPSE_Test_Microbio_fabae
library(gghalves)

f1 <- MPSE_Test_Microbio_fabae %>% 
  mp_plot_alpha(
    .group=Ageing, 
    .alpha=c(Observe, Chao1, Shannon)
  ) + 
  scale_fill_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019")) +
  scale_color_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019"))
 
f1
f1<-f1 + ggtitle("Ageing - A. fabae") +theme_bw(base_size = 18)
f1

ggsave("Alphadiveristy_fabae10032023.png", plot = f1, dpi=500)
###ABundance

MPSE_Test_Microbio_fabae %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Ageing
  )
#MPSE_Test_Microbio_fabae

library(ggalluvial)
library(ggh4x)

# visualize the relative abundance of top 20 phyla for each sample.
p1 <- MPSE_Test_Microbio_fabae %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class = Phylum, 
    topn = 5,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p1<-p1+
  theme(legend.position = "right") + ggtitle("Bacterial phylum proportion into A. fabae Ageing")
p1  


#Family
p2 <- MPSE_Test_Microbio_fabae %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class = Family, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p2<-p2+
  theme(legend.position = "right") + ggtitle("Bacterial family proportion into A. fabae Ageing")
  
p2

p3 <- MPSE_Test_Microbio_fabae%>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class =Genus, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p3<-p3+
  theme(legend.position = "right") + ggtitle("Bacterial genus proportion into A. fabae Ageing")

p3
ggsave("PhylumAphid_fabae10032023.png", plot = p1, dpi=500)
ggsave("FamilyAphid_fabae10032023.png", plot = p2, dpi=500)
ggsave("GenusAphid_fabae10032023.png", plot = p3, dpi=500)




##Heatmap
h1 <- MPSE_Test_Microbio_fabae %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Ageing,
    taxa.class = Family,
    relative = TRUE,
    topn = 10,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h1
ggsave("Heatmap_family_fabae10032023.png", plot = h1, dpi=500)

# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
#aplot::plot_list(gglist=list(h1, h2), tag_levels="A")

#####Beta Diversity

MPSE_Test_Microbio_fabae %<>% 
  mp_decostand(.abundance=Abundance)
MPSE_Test_Microbio_fabae


library(corrr)
# calculate the distance between the samples.
# the distance will be generated a nested tibble and added to the
# colData slot.
MPSE_Test_Microbio_fabae %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
MPSE_Test_Microbio_fabae
# mp_plot_dist provides there methods to visualize the distance between the samples or groups
# when .group is not provided, the dot heatmap plot will be return
p1 <-  MPSE_Test_Microbio_fabae%>% mp_plot_dist(.distmethod = bray)
p1


# when .group is provided, the dot heatmap plot with group information will be return.
p2 <-   MPSE_Test_Microbio_fabae %>% mp_plot_dist(.distmethod = bray, .group = Ageing)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("#19a000", "#00a087FF", "#8700a0","#a00019"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = Species # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )

# when .group is provided and group.test is TRUE, the comparison of different groups will be returned
p3 <-  MPSE_Test_Microbio_fabae%>% mp_plot_dist(.distmethod = bray, .group = Ageing, group.test=TRUE, textsize=2)
p3 

####The PCOA ANalysis#####

MPSE_Test_Microbio_fabae%<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
MPSE_Test_Microbio_fabae

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
MPSE_Test_Microbio_fabae %<>%
  mp_adonis(.abundance=hellinger, .formula=~Ageing, distmethod="bray", permutations=9999, action="add")
MPSE_Test_Microbio_fabae %>% mp_extract_internal_attr(name=adonis)

# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
library(ggside)
library(ggpp)
p2 <- MPSE_Test_Microbio_fabae  %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Ageing,
    show.side = F,
    show.adonis = F,
    .color = Ageing, 
    .size = Observe, 
    ellipse = T,
    show.legend = F # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#36C144", "#0ABDCA", "#8700a0","#CA0A31"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#36C144", "#0ABDCA", "#8700a0","#CA0A31"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
p2
ggsave("PCoA_fabae21122023.png", plot = p2, dpi=500)


######Time A pisum##########
#library(dplyr)

physeq1_pisum <- microViz::ps_filter(physeq1, Species == "A_pisum")
physeq1_pisum


MPSE_Test_pisum<-as.MPSE(physeq1_pisum)
MPSE_Test_pisum

# Rarefied species richness
MPSE_Test_pisum %<>% mp_rrarefy()
# 'chunks' represent the split number of each sample to calculate alpha
# diversity, default is 400. e.g. If a sample has total 40000
# reads, if chunks is 400, it will be split to 100 sub-samples
# (100, 200, 300,..., 40000), then alpha diversity index was
# calculated based on the sub-samples. 
# '.abundance' the column name of abundance, if the '.abundance' is not be 
# rarefied calculate rarecurve, user can specific 'force=TRUE'.
MPSE_Test_pisum %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 100
  )
# The RareAbundanceRarecurve column will be added the colData slot 
# automatically (default action="add")
p3 <- MPSE_Test_pisum%>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Ageing, 
    plot.group = TRUE
  ) +
  scale_color_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019")) +
  scale_fill_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019"))  
p3<-p3 + ggtitle ("A. pisum - Ageing")+ theme_bw(base_size = 18)
p3
ggsave("Curve_Acc_pisum10032023.png", plot = p3, dpi=500)


MPSE_Test_pisum %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
MPSE_Test_pisum
library(gghalves)

f1 <- MPSE_Test_pisum %>% 
  mp_plot_alpha(
    .group=Ageing, 
    .alpha=c(Observe, Chao1, Shannon)
  ) + 
  scale_fill_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019")) +
  scale_color_manual(values=c("#19a000", "#00a087FF", "#8700a0","#a00019"))

f1
f1<-f1 + ggtitle("Ageing - A. pisum") +theme_bw(base_size = 18)
f1

ggsave("Alphadiveristy_pisum10032023.png", plot = f1, dpi=500)
###ABundance

MPSE_Test_pisum %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Ageing
  )
MPSE_Test_pisum

# visualize the relative abundance of top 20 phyla for each sample.
p1 <- MPSE_Test_pisum %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class = Phylum, 
    topn = 5,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p1<-p1+
  theme(legend.position = "right") + ggtitle("Bacterial phylum proportion into A. pisum Ageing")
p1  


#Family
p2 <- MPSE_Test_pisum%>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class = Family, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p2<-p2+
  theme(legend.position = "right") + ggtitle("Bacterial family proportion into A. pisum Ageing")

p2

p3 <- MPSE_Test_pisum%>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Ageing, 
    taxa.class =Genus, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p3<-p3+
  theme(legend.position = "right") + ggtitle("Bacterial genus proportion into A. pisum Ageing")

p3
ggsave("PhylumAphid_pisum10032023.png", plot = p1, dpi=500)
ggsave("FamilyAphid_pisum10032023.png", plot = p2, dpi=500)
ggsave("GenusAphid_pisum10032023.png", plot = p3, dpi=500)

##Heatmap
h1 <- MPSE_Test_pisum %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Ageing,
    taxa.class = Family,
    relative = TRUE,
    topn = 10,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h1

ggsave("Heatmap_familyAphid_pisum10032023.png", plot = h1, dpi=500)

# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
#aplot::plot_list(gglist=list(h1, h2), tag_levels="A")



#####Beta Diversity

MPSE_Test_pisum %<>% 
  mp_decostand(.abundance=Abundance)
MPSE_Test_pisum


library(corrr)
# calculate the distance between the samples.
# the distance will be generated a nested tibble and added to the
# colData slot.
MPSE_Test_pisum %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
MPSE_Test_pisum

# mp_plot_dist provides there methods to visualize the distance between the samples or groups
# when .group is not provided, the dot heatmap plot will be return
p1 <- MPSE_Test_pisum %>% mp_plot_dist(.distmethod = bray)
p1


# when .group is provided, the dot heatmap plot with group information will be return.
p2 <-  MPSE_Test_pisum %>% mp_plot_dist(.distmethod = bray, .group = Ageing)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("orange", "deepskyblue","green","red"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = Species # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )

# when .group is provided and group.test is TRUE, the comparison of different groups will be returned
p3 <-  MPSE_Test_pisum %>% mp_plot_dist(.distmethod = bray, .group = Ageing, group.test=TRUE, textsize=2)
p3 

####The PCOA ANalysis#####

MPSE_Test_pisum%<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
MPSE_Test_pisum

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
MPSE_Test_pisum %<>%
  mp_adonis(.abundance=hellinger, .formula=~Ageing, distmethod="bray", permutations=9999, action="add")
MPSE_Test_pisum %>% mp_extract_internal_attr(name=adonis)


# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
library(ggside)
p2 <- MPSE_Test_pisum  %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Ageing, 
    .color = Ageing, 
    .size = Observe, 
    show.side = F,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#36C144", "#0ABDCA", "#8700a0","#CA0A31"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#36C144", "#0ABDCA", "#8700a0","#CA0A31"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
p2

ggsave("PCoA_pisum21122023.png", plot = p2, dpi=500)
####Hierarchical cluster analysis

MPSE_Test_pisum  %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )
MPSE_Test_pisum 

# if action = 'add', the result of hierarchical cluster will be added to the MPSE object
# mp_extract_internal_attr can extract it. It is a treedata object, so it can be visualized
# by ggtree.
sample.clust <- MPSE_Test_pisum  %>% mp_extract_internal_attr(name='SampleClust')
sample.clust


library(ggtree)
p <- ggtree(sample.clust) + 
  geom_tippoint(aes(color=Ageing)) +
  geom_tiplab(as_ylab = TRUE) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))
p


library(ggtreeExtra)
library(ggplot2)
phyla.tb <- MPSE_Test_pisum %>% 
  mp_extract_abundance(taxa.class=Genus, topn=5)
# The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
phyla.tb



p1 <- p + 
  geom_fruit(
    data=phyla.tb,
    geom=geom_col,
    mapping = aes(x = RelRareAbundanceBySample, 
                  y = Sample, 
                  fill = Phyla
    ),
    orientation = "y",
    #offset = 0.4,
    pwidth = 3, 
    axis.params = list(axis = "x", 
                       title = "The relative abundance of phyla (%)",
                       title.size = 4,
                       text.size = 2, 
                       vjust = 1),
    grid.params = list()
  )
p1










#####################################
#####Pisum vs Fabae at 24h###########
#####################################
library(microViz)
physeq1_24h <- ps_filter(physeq1, Ageing == "24h")

MPSE_Test_Microbio_24h<-as.MPSE(physeq1_24h)
MPSE_Test_Microbio_24h

# Rarefied species richness
MPSE_Test_Microbio_24h %<>% mp_rrarefy()
# 'chunks' represent the split number of each sample to calculate alpha
# diversity, default is 400. e.g. If a sample has total 40000
# reads, if chunks is 400, it will be split to 100 sub-samples
# (100, 200, 300,..., 40000), then alpha diversity index was
# calculated based on the sub-samples. 
# '.abundance' the column name of abundance, if the '.abundance' is not be 
# rarefied calculate rarecurve, user can specific 'force=TRUE'.
MPSE_Test_Microbio_24h %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 100
  )
# The RareAbundanceRarecurve column will be added the colData slot 
# automatically (default action="add")
MPSE_Test_Microbio_24h %>% print(width=180)

p3 <- MPSE_Test_Microbio_24h %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = "Observe", 
    .group = Species, 
    plot.group = TRUE
  ) +
  scale_color_manual(values=c("#46495c", "#14a360")) +
  scale_fill_manual(values=c("#46495c", "#14a360"))  
p3<-p3 + theme_bw(base_size = 18)
p3
ggsave("RarefactionCurve_24h_06032023.png", plot = p3, dpi=500)


MPSE_Test_Microbio_24h%<>% 
  mp_cal_alpha(.abundance=RareAbundance)
MPSE_Test_Microbio_24h
library(gghalves)

f1 <- MPSE_Test_Microbio_24h %>% 
  mp_plot_alpha(
    .group=Species, 
    .alpha=c(Observe, Chao1, Shannon)
  ) +
  scale_fill_manual(values=c("#46495c", "#14a360")) +
  scale_color_manual(values=c("#46495c", "#14a360"))

f1
ggsave("AlphaDiversityMetrics_24h_06032023.png", plot = f1, dpi=500)

###ABundance
MPSE_Test_Microbio_24h %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Species
  )
MPSE_Test_Microbio_24h

# visualize the relative abundance of top 20 phyla for each sample.
p1 <- MPSE_Test_Microbio_24h %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class = Phylum, 
    topn = 5,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p1<-p1+
  theme(legend.position = "right") + ggtitle("Bacterial phylum proportion between A. fabae and A.pisum at 24h")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A.fabae", "A.pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p1
#Family
p2 <- MPSE_Test_Microbio_24h%>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class = Family, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p2<-p2+
  theme(legend.position = "right") + ggtitle("Bacterial family proportion between A. fabae and A.pisum at 24h")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A.fabae", "A.pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p2
#p1+p2

p3 <- MPSE_Test_Microbio_24h %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Species, 
    taxa.class =Genus, 
    topn = 10,
    geom="flowbar",
    relative = T,
    plot.group = T,
    force=T
  )
p3<-p3+
  theme(legend.position = "right") + ggtitle("Bacterial genus proportion between A. fabae and A.pisum at 24h")+
  scale_x_discrete(breaks=c("A_fabae","A_pisum"),
                   labels=c("A.fabae", "A.pisum"))+
  
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(face="italic"   ))
p3
ggsave("PhylumAphid_24h_06032023.png", plot = p1, dpi=500)
ggsave("FamilyAphid_24h_06032023.png", plot = p2, dpi=500)
ggsave("GenusAphid_24h_06032023.png", plot = p3, dpi=500)

##Heatmap
h1 <- MPSE_Test_Microbio_24h %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Species,
    taxa.class = Family,
    relative = TRUE,
    topn = 10,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h1
ggsave("Heatmap_24h_06032023.png", plot = h1, dpi=500)

# the character (scale or theme) of figure can be adjusted by set_scale_theme
# refer to the mp_plot_dist
#aplot::plot_list(gglist=list(h1, h2), tag_levels="A")



#####Beta Diversity

MPSE_Test_Microbio_24h %<>% 
  mp_decostand(.abundance=Abundance)
MPSE_Test_Microbio_24h


library(corrr)
# calculate the distance between the samples.
# the distance will be generated a nested tibble and added to the
# colData slot.
MPSE_Test_Microbio_24h%<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
MPSE_Test_Microbio_24h

# mp_plot_dist provides there methods to visualize the distance between the samples or groups
# when .group is not provided, the dot heatmap plot will be return
p1 <- MPSE_Test_Microbio_24h %>% mp_plot_dist(.distmethod = bray)
p1


# when .group is provided, the dot heatmap plot with group information will be return.
p2 <-  MPSE_Test_Microbio_24h %>% mp_plot_dist(.distmethod = bray, .group = Species)
# The scale or theme of dot heatmap plot can be adjusted using set_scale_theme function.
p2 %>% set_scale_theme(
  x = scale_fill_manual(
    values=c("orange", "deepskyblue"), 
    guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
    )
  ), 
  aes_var = Species # specific the name of variable 
) %>%
  set_scale_theme(
    x = scale_color_gradient(
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  ) %>%
  set_scale_theme(
    x = scale_size_continuous(
      range = c(0.1, 3),
      guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
    ),
    aes_var = bray
  )

# when .group is provided and group.test is TRUE, the comparison of different groups will be returned
p3 <-  MPSE_Test_Microbio_24h%>% mp_plot_dist(.distmethod = bray, .group = Species, group.test=TRUE, textsize=2)
p3 

####The PCOA ANalysis#####

MPSE_Test_Microbio_24h %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
MPSE_Test_Microbio_24h

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
MPSE_Test_Microbio_24h %<>%
  mp_adonis(.abundance=hellinger, .formula=~Species, distmethod="bray", permutations=9999, action="add")
MPSE_Test_Microbio_24h %>% mp_extract_internal_attr(name=adonis)
#ado_test$aov.tab


# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
library(ggside)
p2 <- MPSE_Test_Microbio_24h  %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Species, 
    .color = Species, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#46495c", "#14a360"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#46495c", "#14a360"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
p2<-p2+ggtitle("PCoA at 24h")
p2
ggsave("PCOA_24h_06032023.png", plot = p3, dpi=500)

####Hierarchical cluster analysis

MPSE_Test_Microbio_24h  %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average", # (UPGAE)
    action = "add" # action is used to control which result will be returned
  )
MPSE_Test_Microbio_24h

# if action = 'add', the result of hierarchical cluster will be added to the MPSE object
# mp_extract_internal_attr can extract it. It is a treedata object, so it can be visualized
# by ggtree.
sample.clust <- MPSE_Test_Microbio_24h  %>% mp_extract_internal_attr(name='SampleClust')
sample.clust


library(ggtree)
p <- ggtree(sample.clust) + 
  geom_tippoint(aes(color=Species)) +
  geom_tiplab(as_ylab = TRUE) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))
p


library(ggtreeExtra)
library(ggplot2)
phyla.tb <- MPSE_Test_Microbio_24h%>% 
  mp_extract_abundance(taxa.class=Phylum, topn=5)
# The abundance of each samples is nested, it can be flatted using the unnest of tidyr.
phyla.tb %<>% tidyr::unnest(cols=RareAbundanceBySample) %>% dplyr::rename(Phyla="label")
phyla.tb



p1 <- p + 
  geom_fruit(
    data=phyla.tb,
    geom=geom_col,
    mapping = aes(x = RelRareAbundanceBySample, 
                  y = Sample, 
                  fill = Phyla
    ),
    orientation = "y",
    #offset = 0.4,
    pwidth = 3, 
    axis.params = list(axis = "x", 
                       title = "The relative abundance of phyla (%)",
                       title.size = 4,
                       text.size = 2, 
                       vjust = 1),
    grid.params = list()
  )
p1


#########BIOMARKER Discovery#######

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)
MPSE_Test_Microbio_24h %>% print(width=150)

MPSE_Test_Microbio_24h%<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Species,
    first.test.alpha = 0.05
  )
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- MPSE_Test_Microbio_24h %>% 
  mp_extract_tree(type="taxatree")
taxa.tree

# And the result tibble of different analysis can also be extracted 
# with tidytree (>=0.3.5)
taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Species, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))

##Add Tip Change

#tips<-read.csv2("TpLAbel_tsfo.csv")
#taxa.tree@phylo$tip.label<-tips$TipLabel

# Since taxa.tree is treedata object, it can be visualized by ggtree and ggtreeExtra
p1 <- ggtree(
  taxa.tree,
  layout = "circular",
  size = 0.3
) +
  geom_point(
    data = td_filter(isTip),
    fill="white",
    size=1,
    shape=21
  )
# display the high light of phylum clade.
p1

p2 <- p1 +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  )
p2
# display the relative abundance of features(OTU)

p3 <- p2 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Species, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Species,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))
p3
# display the tip labels of taxa tree
p4 <- p3 + geom_tiplab(size=2, offset=7.2)
p4
# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_Species,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )
p5

# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_Species)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_Species,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))

p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)



#MPSE_Test_Microbio_24h@otutree@phylo$tip.label<-tips$TipLabel
#MPSE_Test_Microbio_24h@taxatree@phylo$tip.label<-tips$TipLabel

p <- MPSE_Test_Microbio_24h %>%
  mp_plot_diff_res(
    group.abun = TRUE,
    pwidth.abun=0.1
  ) +
  scale_fill_manual(values=c("#46495c", "#14a360")) +
  scale_fill_manual(
    aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
    values = c("#46495c", "#14a360")
  ) +
  scale_fill_manual(
    aesthetics = "fill_new_new", # The fill aes for hight light layer of tree was renamed to 'fill_new_new'
    values = c("#E41A1C", "#377EB8", "#4DAF4A",
               "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999","#AD6F3B", "#673770","#D14285"
    )
  )
p + 
  geom_tiplab(size=2, offset=7.2)+ 
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
  )



f <- MPSE_Test_Microbio_24h %>%
  mp_plot_diff_cladogram(
    label.size = 2.5,
    hilight.alpha = .3,
    bg.tree.size = .5,
    bg.point.size = 2,
    bg.point.stroke = .25
  ) +
  scale_fill_diff_cladogram( # set the color of different group.
    values = c("#46495c", "#14a360")
  ) +
  scale_size_continuous(range = c(1, 4))
f
ggsave("cladogram24h.png", plot = f, dpi=500)

f.box <- MPSE_Test_Microbio_24h %>%
  mp_plot_diff_boxplot(
    .group = Species,
  ) %>%
  set_diff_boxplot_color(
    values = c("#46495c", "#14a360"),
    guide = guide_legend(title=NULL)
  )
f.box
f.bar <-  MPSE_Test_Microbio_24h %>%
  mp_plot_diff_boxplot(
    taxa.class = c(Genus, OTU), # select the taxonomy level to display
    group.abun = TRUE, # display the mean abundance of each group
    removeUnknown = TRUE, # whether mask the unknown taxa.
  ) %>%
  set_diff_boxplot_color(
    values = c("#46495c", "#14a360"),
    guide = guide_legend(title=NULL)
  )
f.bar
ff<-aplot::plot_list(f.box, f.bar)
ggsave("bar_plotto_24h.png", plot = ff, dpi=500)

f.mahattan <- MPSE_Test_Microbio_24h %>%
  mp_plot_diff_manhattan(
    .group = Sign_Species,
    .y = fdr,
    .size = 2.4,
    taxa.class = c('OTU', 'Genus'),
    anno.taxa.class = Phylum
  )
f.mahattan
ggsave("f_mahattan.png", plot = f.mahattan, dpi=500)


########
####Analysis of the dynamic of ASVs of the bacteria
#phylum abundance across all samples

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(metagMisc)
library(reshape2)
library(scales)
library(ggpubr)
library(ggradar)
library(cowplot)
#tax_table(DATA)[tax_table(DATA)[,'phylum'] == ""] <- NA 
tax_table(physeq1_filtered)[tax_table(physeq1_filtered)[,'Genus'] == ""] <- NA 

DATA.phyl <- tax_glom(physeq1_filtered, 'Genus', NArm = FALSE)
smp <- rep('sample', 77)
DATA.smp <- merge_samples(DATA.phyl, smp)
DATA.smp  = transform_sample_counts(DATA.smp, function(x) x / sum(x) *100 )



Abund.phyl <- tibble(abund = t(data.frame(otu_table(DATA.smp))),
                     name = get_taxa(DATA.smp, 'Genus')) %>%
  mutate(abund = round(abund, digit = 2),
         name = gsub('g:', '', name)) %>%
  arrange(desc(abund)) 


ggplot(Abund.phyl, aes(x = reorder(name,abund), y = abund)) + geom_point()+
  geom_segment(aes(xend=name, yend=0), color="black") +
  geom_text(aes(name, abund + 5, label = abund), position = position_dodge(width = 1))+
  xlab('Phylum')+
  ylab('Percentage')+
  coord_flip()


#For all the data
DATA.spec <- subset_taxa(physeq1_filtered, Genus != "")
DATA.spec <- tax_glom(DATA.spec, 'Genus')
GenusTable <- data.frame(t(otu_table(DATA.spec)))


#PresAbs <- data.frame(otu_table(DATA.spec))
colnames(GenusTable) <- get_taxa_unique(DATA.spec, 'Genus')
colnames(GenusTable) <- gsub(" ", "", colnames(GenusTable))
DataGenus <- merge(as.data.frame(physeq1_filtered@sam_data),GenusTable,by.x = "row.names", by.y = "row.names")
# Assuming your dataset is named 'my_data' and you want to keep columns 'column1' and 'column2'
DataGenus_filtered <- DataGenus %>%
  select(sample,Ageing,Species,Buchnera,Erwinia,Pseudomonas,Prosthecobacter,Acinetobacter,Staphylococcus,Serratia,Massilia,Carnimonas)

# Assuming your data frame is named 'my_data'
# Remove the letter 'a' from all elements in the column 'my_column'
DataGenus_filtered$Ageing <- gsub("h", "", DataGenus_filtered$Ageing)
DataGenus_filtered$Ageing<- as.numeric(DataGenus_filtered$Ageing)
DataGenus_filtered$Species<- as.factor(DataGenus_filtered$Species)
DataGenus_filtered$sample<- as.factor(DataGenus_filtered$sample)
str(DataGenus_filtered)

###Model GLMM###
library(lme4)
library(multcomp)
library(glmmTMB)
library(DHARMa)
#Useful script to source
#les scripts de Gilles San Martin
source("C://Users/grego/OneDrive/Assistanat/2018-2019/MesureBiodiversite/Analyses spatiales 2/mytoolbox.R")
source("C://Users/grego/OneDrive/Assistanat/2018-2019/MesureBiodiversite/Analyses spatiales 2/model.select_0.4.R")
####
###GLMM Formulae
glmm_Buchnera<-glmer(Buchnera ~ Ageing + Species + (1|sample) , data=DataGenus_filtered, family = poisson)
summary(glmm_Buchnera)

ggplot(data=DataGenus_filtered, aes(y = Buchnera, x = Ageing, fill = Species)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = T, method.args = list(family = "poisson"))+
  theme_bw()

###By APhid
DataGenus_filtered_fabae <- subset(DataGenus_filtered, Species == "A_fabae")
DataGenus_filtered_pisum <- subset(DataGenus_filtered, Species == "A_pisum")

###GLMM Formulae
#Buchnera_fabae #Cela ne change pas au cours des heures
glmm_Buchnera<-glmer(Buchnera ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Buchnera)

ggplot(data=DataGenus_filtered_fabae, aes(y = Buchnera, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Erwinia_fabae ##Il y a un pic d'Erwinia à 48h
glmm_Erwinia<-glmer(Erwinia ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Erwinia)

ggplot(data=DataGenus_filtered_fabae, aes(y = Erwinia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Pseudomonas_fabae #Il n'y a pas de changement au cours du temps
glmm_Pseudomonas<-glmer(Pseudomonas ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Pseudomonas)

ggplot(data=DataGenus_filtered_fabae, aes(y = Pseudomonas, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Prosthecobacter_fabae #Il y a une augmentation significative des Prosthecobacter au cours du temps chez fabae
glmm_Prosthecobacter<-glmer(Prosthecobacter ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Prosthecobacter)

ggplot(data=DataGenus_filtered_fabae, aes(y = Prosthecobacter, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Acinetobacter_fabae #Il y a une augmentation significative des Acinetobacter au cours du temps chez fabae
glmm_Acinetobacter<-glmer(Acinetobacter ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Acinetobacter)

ggplot(data=DataGenus_filtered_fabae, aes(y = Acinetobacter, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Staphylococcus_fabae #Il y a une augmentation significative des Staphylococcus au cours du temps chez fabae
glmm_Staphylococcus<-glmer(Staphylococcus ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Staphylococcus)

ggplot(data=DataGenus_filtered_fabae, aes(y = Staphylococcus, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Serratia_fabae #Pas de changement
glmm_Serratia<-glmer(Serratia ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Serratia)

ggplot(data=DataGenus_filtered_fabae, aes(y = Serratia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Massilia_fabae #Pas de changement
glmm_Massilia<-glmer(Massilia ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Massilia)

ggplot(data=DataGenus_filtered_fabae, aes(y = Massilia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Carnimonas_fabae #Pas de changement
glmm_Carnimonas<-glmer(Carnimonas ~ Ageing  + (1|sample) , data=DataGenus_filtered_fabae, family = poisson)
summary(glmm_Carnimonas)

ggplot(data=DataGenus_filtered_fabae, aes(y = Carnimonas, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()




##########pisum






#Buchnera_pisum #Cela ne change pas au cours des heures
glmm_Buchnera<-glmer(Buchnera ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Buchnera)

ggplot(data=DataGenus_filtered_pisum, aes(y = Buchnera, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Erwinia_pisum ##Cela ne change pas au cours des heures
glmm_Erwinia<-glmer(Erwinia ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Erwinia)

ggplot(data=DataGenus_filtered_pisum, aes(y = Erwinia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Pseudomonas_pisum #Il y a une augmentation au cours du temps
glmm_Pseudomonas<-glmer(Pseudomonas ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Pseudomonas)

ggplot(data=DataGenus_filtered_pisum, aes(y = Pseudomonas, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Prosthecobacter_pisum #Cela ne change pas au cours des heures
glmm_Prosthecobacter<-glmer(Prosthecobacter ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Prosthecobacter)

ggplot(data=DataGenus_filtered_pisum, aes(y = Prosthecobacter, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Acinetobacter_pisum #Cela ne change pas au cours des heures
glmm_Acinetobacter<-glmer(Acinetobacter ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Acinetobacter)

ggplot(data=DataGenus_filtered_pisum, aes(y = Acinetobacter, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Staphylococcus_pisum #Presque significatif
glmm_Staphylococcus<-glmer(Staphylococcus ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Staphylococcus)

ggplot(data=DataGenus_filtered_pisum, aes(y = Staphylococcus, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Serratia_pisum #Pas de changement
glmm_Serratia<-glmer(Serratia ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Serratia)

ggplot(data=DataGenus_filtered_pisum, aes(y = Serratia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Massilia_pisum #Pas de changement
glmm_Massilia<-glmer(Massilia ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Massilia)

ggplot(data=DataGenus_filtered_pisum, aes(y = Massilia, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()

#Carnimonas_pisum #NA trop peu de données
glmm_Carnimonas<-glmer(Carnimonas ~ Ageing  + (1|sample) , data=DataGenus_filtered_pisum, family = poisson)
summary(glmm_Carnimonas)

ggplot(data=DataGenus_filtered_pisum, aes(y = Carnimonas, x = Ageing)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm",se = F, method.args = list(family = "poisson"))+
  theme_bw()
