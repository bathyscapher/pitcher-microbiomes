### Read taxa (ASV), add metadata and much more

library("ade4")
library("cocorresp")
library("phyloseq")
library("ggplot2")
library("GGally")
theme_set(theme_bw(base_size = 20) +
            theme(rect = element_rect(fill = "transparent")))
library("reshape2")


rm(list = ls()); gc()
set.seed(34706)


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")


################################################################################
### Moss taxa
## Merged
moss.prok.a <- readRDS("rds/SMP_moss20.prok.a.RDS")
moss.euk.a <- readRDS("rds/SMP_moss20.euk.a.RDS")


################################################################################
### nMDS pitchers
## Choose pro- or eukaryotes
primer <- "16S"
# primer <- "18S"


if(primer == "16S") {
  smp <- prok.a
} else if(primer == "18S") {
  smp <- euk.a
}


## log-transform
smp.log <- transform_sample_counts(smp, function(otu) {log1p(otu)})


## nMDS
smp.nmds <- ordinate(smp.log, method = "NMDS", distance = "bray", k = 4,
                     autotransform = FALSE, trymax = 100)
smp.nmds
stressplot(smp.nmds)


## Plot nMDS
plot_ordination(smp, smp.nmds, shape = "Succession", color = "Site",
                title = NULL, axes = 1:2) +
  stat_ellipse(aes(group = Site), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("img/SMP_16S_nMDS_ASV.pdf", width = 8.27, height = 8.27)
# ggsave("img/SMP_18S_nMDS_ASV.pdf", width = 8.27, height = 8.27)


### nMDS pitchers and mosses
if(primer == "16S") {
  smp.moss <- merge_phyloseq(prok.a, moss.prok.a)
} else if(primer == "18S") {
  smp.moss <- merge_phyloseq(euk.a, moss.euk.a)
}


## log-transform
smp.moss.log <- transform_sample_counts(smp.moss, function(otu) {log1p(otu)})


## nMDS
smp.moss.nmds <- ordinate(smp.moss.log, method = "NMDS", distance = "bray",
                          k = 4, autotransform = FALSE, trymax = 100,
                          weakties = FALSE)
smp.moss.nmds
stressplot(smp.moss.nmds)


## Plot nMDS
plot_ordination(smp.moss, smp.moss.nmds,
                shape = "Site", color = "Succession",
                title = NULL, axes = 1:2) +
  stat_ellipse(aes(group = Succession), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = gradCol[c(2, 1, 4)]) +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("img/SMP_16S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)
# ggsave("img/SMP_18S_nMDS_withMoss.pdf", width = 8.27, height = 8.27)





################################################################################
### Symmetric CoCA #############################################################
## Get ASVs
otu.p <- data.frame(otu_table(prok.a))
otu.p <- otu.p[!(row.names(otu.p) %in% "LT4a06E"), ]
otu.e <- data.frame(otu_table(euk.a))
# setdiff(rownames(otu.p), rownames(otu.e))


## sCoCA
source("~/Seafile/SarraceniaMicrobiomeProject/analysis/CoCa_Alric/sCoCA-permutation.R")


sCoCA.smp <- coca(otu.p ~ ., otu.e, method = "symmetric", symmetric = TRUE)
summary(sCoCA.smp)


## Permutation test to evaluate the significance of the association between two
## tables
eig.test <- randtest.coca(X = otu.p, Y = otu.e, nrepet = 999)


## Amount of covariance explained by each axis of the sCoCA. Keep 3 axes it says
# pdf("SMP_sCoCA_eigen.pdf", height = 4, width = 11.69)
screeplot(sCoCA.smp, xlab = "sCoCA axes")
# dev.off()


## Common variance. sum(sCoCA.smp$lambda) = "trace" = "total inertia"
100 * (sCoCA.smp$lambda[1:3]) / sum(sCoCA.smp$lambda)
sum(100 * (sCoCA.smp$lambda[1:3]) / sum(sCoCA.smp$lambda))


## Correlation between CoCA axes of X and Y
head(corAxis(sCoCA.smp))


## Refit
sCoCA.smp4 <- coca(otu.p ~ ., otu.e, method = "symmetric", symmetric = TRUE,
                   n.axes = 3)


## Variance explained for each table
eig.p <- 100 * (sCoCA.smp4$inertia$total$Y - sCoCA.smp4$inertia$residual$Y) /
  sCoCA.smp4$inertia$total$Y
eig.e <- 100 * (sCoCA.smp4$inertia$total$X - sCoCA.smp4$inertia$residual$X) /
  sCoCA.smp4$inertia$total$X
eig.e; eig.p






################################################################################
### Plot taxa
## Load pro- and eukaryotic taxa together
smp <- readRDS("rds/SMP_smp.RDS")


### Dataframe with taxa and abundances
## Pitchers
tax.l <- as.data.frame(smp@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(smp)))

tax.otu.l <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.l) <- tax.otu.l$Row.names
colnames(tax.otu.l)[1] <- "Taxa"

rm(tax.l, otu.l)


### log-transform abundances
tax.otu.l[, -c(1:7)] <- log1p(tax.otu.l[, -c(1:7)])


tax.otu.lm <- reshape2::melt(tax.otu.l,
                             id.vars = c("Taxa", "Domain", "Phylum", "Class",
                                         "Order", "Family", "Genus"),
                             variable.name = "FullID", value.name = "Count")
tax.otu.lm$Site <- substr(tax.otu.lm$FullID, 1, 2)
tax.otu.lm$Site <- ordered(tax.otu.lm$Site,
                           levels = c("CB", "LT", "LE", "LV", "LM"))
tax.otu.lm$Succession <- as.factor(substr(tax.otu.lm$FullID, 7, 7))
levels(tax.otu.lm$Succession) <- c("Early", "Late")
tax.otu.lm$Domain[tax.otu.lm$Domain == "Archaea"] <- "A."


## Reorder by taxonomy
tax.otu.lm[, c(2:7)] <- lapply(tax.otu.lm[, c(2:7)], as.factor)


tax <- c("Thaumarchaeota", ## Archaea
         ## Proteobacteria
         "Proteobacteria", "Epsilonbacteraeota",
         "Chlamydiae", "Planctomycetes", "Verrucomicrobia", ## PVC group
         "Gemmatimonadetes", ## FCB group
         ## Singles
         "Acidobacteria", "Bacteroidetes", "Chloroflexi", "Dependentiae",
         ## Terrabacteria
         "Actinobacteria", "Armatimonadetes", "Cyanobacteria", "Firmicutes",
         "WPS-2", ## ?
         "Schizoplasmodiida", "Discosea", ## Amoeba
         ## SAR:Heterokonta
         "Bicosoecida", "Ochrophyta", "Peronosporomycetes",
         "Cercozoa", ## SAR: Rhizaria
         ## Harosa:Alveolata
         "Ciliophora", "Dinoflagellata", "Protalveolata",
         "Euglenozoa", "Heterolobosea", ## Excavata
         ## Algae
         "Chlorophyta_ph", "Phragmoplastophyta",
         "Cryptomonadales", "Florideophycidae",
         "Nucleariidae_and_Fonticula_group", ## Holomycota
         ## Fungi
         "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota",
         "Opisthokonta_ph", "Cryptomycota", "Zoopagomycota",
         "Choanoflagellida", ## Choanozoa
         ## Metazoa
         "Annelida", "Platyhelminthes", "Nematoda", "Rotifera", "Tardigrada",
         "Arthropoda")

## Pronounce small values (otherwise they will be barely visible)
summary(tax.otu.lm$Count[tax.otu.lm$Count < 0.01])
tax.otu.lm$Count[tax.otu.lm$Count < 0.01 & tax.otu.lm$Count > 0] <- 0.01
summary(tax.otu.lm$Count)




## Plot
ggplot(tax.otu.lm, aes(#x = factor(Phylum),
                       x = factor(Phylum, level = rev(tax)),
                       y = Count, fill = Site)) +
  geom_bar(stat = "identity") +
  facet_grid(Domain ~ Succession, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(paste(ln[e], "(Relative abundance + 1)"))) +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 15)
# ggsave("img/SMP_PhylaOverview.pdf", width = 8.27, height = 11.69)


aggregate(Count ~ Class, data = tax.otu.lm, FUN = sum)


### Inspect taxa and frequencies
unique(tax.otu.l$Domain)

dim(tax.otu.l[tax.otu.l$Domain == "Bacteria", ])
dim(tax.otu.l[tax.otu.l$Domain == "Eukaryota", ])
unique(tax.otu.l$Phylum[tax.otu.l$Domain == "Bacteria"])
unique(tax.otu.l$Phylum[tax.otu.l$Domain == "Eukaryota"])

df <- tax.otu.l[tax.otu.l$Phylum == "Basidiomycota", ]
summary(colSums(df[, -c(1:7)]) == 0)
plot(colSums(df[, -c(1:7)]))

100 / 160 * 5


### Plot incidence
## Presence/absence
tax.otu.pa <- tax.otu.l
tax.otu.pa[, -c(1:7)][tax.otu.pa[, -c(1:7)] > 0] <- 1



## Sum by Phylum and replace with 1
inc <- aggregate(tax.otu.pa[, -c(1:7)],
                 by = list(tax.otu.pa$Domain, tax.otu.pa$Phylum), FUN = sum)
names(inc)[1:2] <- c("Domain", "Phylum")

inc[, -c(1:2)][inc[, -c(1:2)] > 0] <- 1
inc$Incidence <- rowSums(inc[, -c(1:2)])


## Sort
inc <- inc[order(inc$Incidence), ]
inc$Phylum <- factor(inc$Phylum, levels = unique(as.character(inc$Phylum)))
inc$IncidencePercent <- 100 / 160 * inc$Incidence

inc$Incidence[inc$Domain == "Archaea"]
inc$Incidence[inc$Phylum == "Actinobacteria"]
inc$Incidence[inc$Phylum == "Basidiomycota"]


tax.otu.l$Genus[tax.otu.l$Phylum == "Basidiomycota"]


# euks <- tax.otu.l[tax.otu.l$Domain == "Eukaryota", ]
# sort(unique(euks$Phylum))
#
# algae <- inc[inc$Phylum %in% c("Chlorophyta_ph", "Phragmoplastophyta", "Cryptomonadales", "Florideophycidae"), ]
# algae <- euks[euks$Phylum %in% c("Chlorophyta_ph", "Phragmoplastophyta", "Cryptomonadales", "Florideophycidae"), ]
#
# Chlorophyta_ph
# Florideophycidae

## Abbreviate archaea
inc$Domain[inc$Domain == "Archaea"] <- "A."


## Plot
ggplot(inc, aes(x = Phylum, y = IncidencePercent)) +
  geom_point() +
  facet_grid(Domain ~ ., scales = "free", space = "free_y") +
  theme(legend.position = "none") +
  ylab("Incidence [%]") +
  xlab("") +
  coord_flip() +
  theme_bw(base_size = 10)
# ggsave("img/SMP_Incidence.pdf", width = 8.27, height = 11.69 / 2)
