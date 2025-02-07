### SMP Pitcher microbiome and prey composition


library("phyloseq")
library("ggplot2")
theme_set(theme_bw(base_size = 15) +
            theme(rect = element_rect(fill = "transparent")))
library("vegan")
library("cocorresp")


rm(list = ls()); gc()
set.seed(34706)


################################################################################
## Colorblind temperature scale for the 5 sites
gradCol <- c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9")


################################################################################
### Pitcher data
prok.a <- readRDS("rds/SMP_prok.a.RDS")
euk.a <- readRDS("rds/SMP_euk.a.RDS")
smp <- readRDS("rds/SMP_smp.RDS")


smpMeta.df <- data.frame(sample_data(smp))


################################################################################
## Samples without prey

empty <- smpMeta.df[smpMeta.df$prey.items == 0, ]
rownames(empty)


smp.e <- prune_samples(smp@sam_data$FullID %in% rownames(empty), smp)


## Remove empty taxa
any(taxa_sums(smp.e) == 0)
sum(taxa_sums(smp.e) == 0)
smp.e <- prune_taxa(taxa_sums(smp.e) > 0, smp.e) # remove unobserved
smp.e
# tax_table(smp.e)


################################################################################
### Plot taxa
tax.l <- as.data.frame(smp.e@tax_table@.Data)
otu.l <- t(as.data.frame(otu_table(smp.e)))

tax.otu.l <- merge(tax.l, otu.l, by = 0, all = TRUE)
rownames(tax.otu.l) <- tax.otu.l$Row.names
colnames(tax.otu.l)[1] <- "Taxa"

rm(tax.l, otu.l)


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
tax.otu.lm$FullID <- ordered(tax.otu.lm$FullID,
                             levels = c("LT4a06E", "LT4a07E", "LE3a14L",
                                        "LE3a15E", "LV2b05E", "LV2b09E",
                                        "LV2b10E"))


## Reorder taxonomically
tax.otu.lm[, c(2:7)] <- lapply(tax.otu.lm[, c(2:7)], as.factor)

# setdiff(tax, sort(unique(tax.otu.lm$Phylum)))

tax <- c("Thaumarchaeota", ## Archaea
         ## Proteobacteria
         "Proteobacteria",
         "Chlamydiae", "Planctomycetes", "Verrucomicrobia", ## PVC group
         "Gemmatimonadetes", ## FCB group
         ## Singles
         "Acidobacteria", "Bacteroidetes", "Chloroflexi",
         ## Terrabacteria
         "Actinobacteria", "Armatimonadetes", "Cyanobacteria", "Firmicutes",
         ## SAR:Heterokonta
         "Ochrophyta", "Peronosporomycetes",
         "Cercozoa", ## SAR: Rhizaria
         ## Harosa:Alveolata
         "Ciliophora",
         "Euglenozoa", ## Excavata
         ## Algae
         "Chlorophyta_ph",
         ## Fungi
         "Basidiomycota",  "Chytridiomycota",
         "Zoopagomycota",
         "Choanoflagellida", ## Choanozoa
         ## Metazoa
         "Arthropoda")


## Pronounce small values even more
tax.otu.lm$Count[tax.otu.lm$Count < 0.1] <- 0.01
summary(tax.otu.lm$Count)


## Plot
ggplot(tax.otu.lm, aes(x = factor(Phylum, level = rev(tax)), y = Count,
                       fill = Site)) +
  geom_bar(stat = "identity") +
  facet_grid(Domain ~ FullID, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = gradCol) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45,
                                                             hjust = 1)) +
  ylab(expression(paste(ln[e], "(Relative abundance + 1)"))) +
  xlab("") +
  coord_flip()
# ggsave("img/SMP_Preyfree.pdf", width = 8.27, height = 6)


################################################################################
### nMDS
## Select domain
# ps <- prok.a
ps <- euk.a


## Group insect counts
sample_data(ps)$prey.items.grp <- sample_data(ps)$prey.items
sample_data(ps)$prey.items.grp[sample_data(ps)$prey.items > 0 &
                                  sample_data(ps)$prey.items < 6] <- "1 - 5"
sample_data(ps)$prey.items.grp[sample_data(ps)$prey.items > 5 &
                                  sample_data(ps)$prey.items < 21] <- "6 - 20"
sample_data(ps)$prey.items.grp[sample_data(ps)$prey.items > 20 &
                                  sample_data(ps)$prey.items < 51] <- "21 - 50"
sample_data(ps)$prey.items.grp[sample_data(ps)$prey.items > 50 &
                                  sample_data(ps)$prey.items < 71] <- "51 - 70"
sample_data(ps)$prey.items.grp[sample_data(ps)$prey.items > 70] <- "> 70"
table(sample_data(ps)$prey.items.grp)


sample_data(ps)$prey.items.grp <- as.factor(sample_data(ps)$prey.items.grp)
# levels(sample_data(ps)$prey.items.grp) <- c("0", "1 - 5", "6 - 20", "21 - 50",
#                                         "51 - 70", "> 70")


smp.log <- transform_sample_counts(ps, function(otu) {log1p(otu)})


## Ordinate
smp.nmds <- ordinate(smp.log, method = "NMDS", distance = "euclidean", k = 4,
                     autotransform = FALSE, trymax = 100)
smp.nmds
stressplot(smp.nmds)


plot_ordination(ps, smp.nmds, color = "prey.items", #shape = "Site",
                title = NULL, axes = 1:2) +
  stat_ellipse(aes(group = prey.items.grp), type = "t", linetype = 2, size = 0.2) +
  geom_point(size = 3) +
  scale_color_continuous(type = "viridis") +
  coord_fixed(ratio = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "vertical") +
  xlab("nMDS1") +
  ylab("nMDS2")
# ggsave("img/SMP_16S_nMDS_Prey.pdf", width = 8.27, height = 8.27)
# ggsave("img/SMP_18S_nMDS_Prey.pdf", width = 8.27, height = 8.27)


################################################################################
################################################################################
### Predictive co-correspondence analysis
set.seed(14262)


prey <- read.table("csv/SMP_Prey_clean.csv",
                   sep = "\t", header = TRUE)


rownames(prey) <- prey$FullID
prey$FullID <- NULL
prey$Site <- NULL


## Get non-empty samples
not.empty <- rownames(prey[rowSums(prey) > 0, ])


## Remove empty columns and log transform
prey <- prey[, colSums(prey) > 0]
prey <- prey[rowSums(prey) > 0, ]
prey <- log1p(prey)


## Pro- or eukaryotes
# smp.a <- prok.a
smp.a <- euk.a


## Subset microbiome to samples with prey only
smp.ne <- prune_samples(smp.a@sam_data$FullID %in% not.empty, smp.a)
smp.otu <- data.frame(otu_table(smp.ne))


smp.otu <- smp.otu[, colSums(smp.otu) > 0]
smp.otu <- smp.otu[rowSums(smp.otu) > 0, ]


## log transform
smp.otu <- log1p(smp.otu)



################################################################################
# predictive coca
pCoCA.smp <- coca(smp.otu ~ ., data = prey, method = "predictive",
                  reg.method = "simpls")

# Explained variance
summary(pCoCA.smp)


cross.smp <- crossval(smp.otu, prey)
which.max(cross.smp$CVfit)
cross.smp


pCoCA.smp.perm <- permutest(pCoCA.smp, permutations = 999)
pCoCA.smp.perm


pCoCA.smp.1 <- coca(smp.otu ~ ., data = prey, method = "predictive",
                    reg.method = "simpls", n.axes = 1)
summary(pCoCA.smp.1)
