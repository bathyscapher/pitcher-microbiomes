### Structural equation modeling

library("gplots")
library("lavaan")
library("semPlot")


rm(list = ls()); gc()
set.seed(34706)


smpMeta.df <- read.table("csv/SMP_Metadata_SEM_scaled.csv", sep = "\t")


dim(smpMeta.df) # should be a data.frame with 160 rows and 29 columns
colnames(smpMeta.df)


################################################################################
### Model with composite
### Composite "microhabitat"
micro <-
'mb.leaf.prok ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items
mb.leaf.euk ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items
prey ~ SampleVolume_mL + Age_d + SampleColorDummy + pH + prey.items'


fit.micro <- sem(micro, data = smpMeta.df)
summary(fit.micro, rsq = TRUE)


## Extract coefficients
SampleVolume.lp <- coef(fit.micro)[[1]]
Age_d.lp <- coef(fit.micro)[[2]]
SampleColorDummy.lp <- coef(fit.micro)[[3]]
pH.lp <- coef(fit.micro)[[4]]
prey.items.lp <- coef(fit.micro)[[5]]


coef(fit.micro)[1:5]; SampleVolume.lp; Age_d.lp; SampleColorDummy.lp; pH.lp; prey.items.lp


SampleVolume.le <- coef(fit.micro)[[6]]
Age_d.le <- coef(fit.micro)[[7]]
SampleColorDummy.le <- coef(fit.micro)[[8]]
pH.le <- coef(fit.micro)[[9]]
prey.items.le <- coef(fit.micro)[[10]]

coef(fit.micro)[6:10]; SampleVolume.le; Age_d.le; SampleColorDummy.le; pH.le; prey.items.le


SampleVolume.p <- coef(fit.micro)[[11]]
Age_d.p <- coef(fit.micro)[[12]]
SampleColorDummy.p <- coef(fit.micro)[[13]]
pH.p <- coef(fit.micro)[[14]]
prey.items.p <- coef(fit.micro)[[15]]

coef(fit.micro)[11:15]; SampleVolume.p; Age_d.p; SampleColorDummy.p; pH.p; prey.items.p


## Compute composite "microhabitat"
smpMeta.df$comp.micro.prok <- SampleVolume.lp * smpMeta.df$SampleVolume_mL +
  Age_d.lp * smpMeta.df$Age_d +
  prey.items.lp * smpMeta.df$prey.items +
  SampleColorDummy.lp * smpMeta.df$SampleColorDummy +
  pH.lp * smpMeta.df$pH

smpMeta.df$comp.micro.euk <- SampleVolume.le * smpMeta.df$SampleVolume_mL +
  Age_d.le * smpMeta.df$Age_d +
  prey.items.le * smpMeta.df$prey.items +
  SampleColorDummy.le * smpMeta.df$SampleColorDummy +
  pH.le * smpMeta.df$pH

smpMeta.df$comp.micro.prey <- SampleVolume.p * smpMeta.df$SampleVolume_mL +
  Age_d.p * smpMeta.df$Age_d +
  prey.items.p * smpMeta.df$prey.items +
  SampleColorDummy.p * smpMeta.df$SampleColorDummy +
  pH.p * smpMeta.df$pH


rm(SampleVolume.le, Age_d.le, prey.items.le, SampleColorDummy.le, pH.le,
   SampleVolume.lp, Age_d.lp, prey.items.lp, SampleColorDummy.lp, pH.lp,
   SampleVolume.p, Age_d.p, prey.items.p, SampleColorDummy.p, pH.p)


### Run SEM with the manually computed composite
comp.micro <-
'mb.leaf.prok ~ comp.micro.prok
mb.leaf.euk ~ comp.micro.euk
prey ~ comp.micro.prey'


fit.micro.comp <- sem(comp.micro, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.micro.comp, rsq = TRUE, standardized = TRUE)
summary(fit.micro, rsq = TRUE, standardized = TRUE)


rm(fit.micro.comp, micro, comp.micro)


### Composite "macrohabitat": estimate separate composites for the response
### variables
macro <-
  'mb.leaf.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
mb.leaf.euk ~ Distance + MeanTemperature + Altitude + CanopyCover
prey ~ Distance + MeanTemperature + Altitude + CanopyCover
#mb.moss.prok ~ Distance + MeanTemperature + Altitude + CanopyCover
#mb.moss.euk ~ Distance + MeanTemperature + Altitude + CanopyCover
'


fit.macro <- sem(macro, data = smpMeta.df)
summary(fit.macro, rsq = TRUE, fit.measures = TRUE)


## Extract coefficients
## For the pitchers/leaves
Distance.lp <- coef(fit.macro)[[1]]
meanTemperature.lp <- coef(fit.macro)[[2]]
Altitude.lp <- coef(fit.macro)[[3]]
CanopyCover.lp <- coef(fit.macro)[[4]]

coef(fit.macro)[1:4]; Distance.lp; meanTemperature.lp; Altitude.lp; CanopyCover.lp


Distance.le <- coef(fit.macro)[[5]]
meanTemperature.le <- coef(fit.macro)[[6]]
Altitude.le <- coef(fit.macro)[[7]]
CanopyCover.le <- coef(fit.macro)[[8]]

coef(fit.macro)[5:8]; Distance.le; meanTemperature.le; Altitude.le; CanopyCover.le


## For the prey
Distance.p <- coef(fit.macro)[[9]]
meanTemperature.p <- coef(fit.macro)[[10]]
Altitude.p <- coef(fit.macro)[[11]]
CanopyCover.p <- coef(fit.macro)[[12]]

coef(fit.macro)[9:12]; Distance.p; meanTemperature.p; Altitude.p; CanopyCover.p


# ## For the mosses
# Distance.mp <- coef(fit.macro)[[13]]
# meanTemperature.mp  <- coef(fit.macro)[[14]]
# Altitude.mp <- coef(fit.macro)[[15]]
# CanopyCover.mp <- coef(fit.macro)[[16]]
#
# coef(fit.macro)[13:16]; Distance.mp; meanTemperature.mp; Altitude.mp; CanopyCover.mp
#
#
# Distance.me <- coef(fit.macro)[[17]]
# meanTemperature.me  <- coef(fit.macro)[[18]]
# Altitude.me <- coef(fit.macro)[[19]]
# CanopyCover.me <- coef(fit.macro)[[20]]
#
# coef(fit.macro)[17:20]; Distance.me; meanTemperature.me; Altitude.me; CanopyCover.me


### Compute composites "macrohabitat" for the leaves, the prey and the moss
## For the leaves
smpMeta.df$comp.macro.lp <- Distance.lp * smpMeta.df$Distance +
  meanTemperature.lp * smpMeta.df$MeanTemperature +
  Altitude.lp * smpMeta.df$Altitude +
  CanopyCover.lp * smpMeta.df$CanopyCover

smpMeta.df$comp.macro.le <- Distance.le * smpMeta.df$Distance +
  meanTemperature.le * smpMeta.df$MeanTemperature +
  Altitude.le * smpMeta.df$Altitude +
  CanopyCover.le * smpMeta.df$CanopyCover

## For the prey
smpMeta.df$comp.macro.prey <- Distance.p * smpMeta.df$Distance +
  meanTemperature.p * smpMeta.df$MeanTemperature +
  Altitude.p * smpMeta.df$Altitude +
  CanopyCover.p * smpMeta.df$CanopyCover

## For the moss
# smpMeta.df$comp.macro.mp <- Distance.mp * smpMeta.df$Distance +
#   meanTemperature.mp * smpMeta.df$MeanTemperature +
#   Altitude.mp * smpMeta.df$Altitude +
#   CanopyCover.mp * smpMeta.df$CanopyCover
#
# smpMeta.df$comp.macro.me <- Distance.me * smpMeta.df$Distance +
#   meanTemperature.me * smpMeta.df$MeanTemperature +
#   Altitude.me * smpMeta.df$Altitude +
#   CanopyCover.me * smpMeta.df$CanopyCover


rm(Distance.le, meanTemperature.le, Altitude.le, CanopyCover.le, Distance.lp,
   meanTemperature.lp, Altitude.lp, CanopyCover.lp, Distance.p,
   meanTemperature.p, Altitude.p,
   # Distance.mp, meanTemperature.mp, Altitude.mp,
   # CanopyCover.mp, Distance.me, meanTemperature.me, Altitude.me,
   # CanopyCover.me,
   CanopyCover.p)


### SEM with composite macrohabitat
comp.macro1 <-
  'mb.leaf.prok ~ comp.macro.lp
mb.leaf.euk ~ comp.macro.le
prey ~ comp.macro.prey'

comp.macro2 <-
  'mb.moss.prok ~ comp.macro.mp
mb.moss.euk ~ comp.macro.me'


fit.macro.comp1 <- sem(comp.macro1, data = smpMeta.df)
# fit.macro.comp2 <- sem(comp.macro2, data = smpMeta.df)


## Control if the standard coefficient and z-value (standardised coefficient)
## are similar
summary(fit.macro.comp1, rsq = TRUE, fit.measures = TRUE)
# summary(fit.macro.comp2, rsq = TRUE, fit.measures = TRUE)
summary(fit.macro, rsq = TRUE, fit.measures = TRUE)


rm(macro, comp.macro1)#, comp.macro2)


################################################################################
### Correlation and covariance
smpMeta.sub <- smpMeta.df[c("Site",
                            "mb.leaf.prok", "comp.macro.lp", "comp.micro.prok",
                            "mb.leaf.euk", "comp.macro.le", "comp.micro.euk",
                            "mb.moss.prok", #"comp.macro.mp",
                            "mb.moss.euk", #"comp.macro.me",
                            "prey.items", "prey",
                            "comp.macro.prey", "comp.micro.prey")]


cor_pal <- colorRampPalette(c("red", "white",  "blue"))(n = 19)
smpMeta.sub <- subset(smpMeta.sub, select = -c(Site))


## Sort correlation data.frame
cors <- data.frame(cor(smpMeta.sub, use = "complete"))
cors <- cors[c("mb.leaf.prok", "mb.leaf.euk", "mb.moss.prok", "mb.moss.euk",
               "comp.macro.lp", "comp.macro.le", "comp.macro.prey",
               # "comp.macro.mp", "comp.macro.me",
               "comp.micro.prok", "comp.micro.euk", "comp.micro.prey",
               "prey", "prey.items")]
cors <- cors[colnames(cors), ]


pdf("img/SMP_Metadata_Correlation.pdf", height = 8.27, width = 8.27)
heatmap.2(as.matrix(cors), dendrogram = "none", key = TRUE, key.title = "",
          cellnote = round(cors, digits = 2), notecol = "black", notecex = 0.55,
          trace = "none", distfun = function(x) {x}, cexCol = 1.4, cexRow = 1.4,
          symm = TRUE, margins = c(16, 16), col = cor_pal,
          tracecol = "black", Rowv = FALSE, Colv = FALSE)
dev.off()


## Sort covariance data.frame
covs <- data.frame(cov(smpMeta.sub, use = "complete"))
covs <- covs[c("mb.leaf.prok", "mb.leaf.euk", "mb.moss.prok", "mb.moss.euk",
               "comp.macro.lp", "comp.macro.le", "comp.macro.prey",
               # "comp.macro.mp", "comp.macro.me",
               "comp.micro.prok", "comp.micro.euk", "comp.micro.prey",
               "prey", "prey.items")]
covs <- covs[colnames(covs), ]


pdf("img/SMP_Metadata_Covariance.pdf", height = 8.27, width = 8.27)
heatmap.2(as.matrix(covs), dendrogram = "none", key = TRUE, key.title = "",
          cellnote = round(covs, digits = 2), notecol = "black", notecex = 0.65,
          trace = "none", distfun = function(x) {x}, cexCol = 1.4, cexRow = 1.4,
          symm = TRUE, margins = c(16, 16), col = cor_pal,
          tracecol = "black", Rowv = FALSE, Colv = FALSE)
dev.off()


################################################################################
### Run SEM with the two composites
## With linking prey to macrohabitat
comp.w <-
'prey ~ comp.macro.prey + comp.micro.prey
# mb.moss.prok ~ comp.macro.mp
# mb.moss.euk ~ comp.macro.me

mb.leaf.prok ~ comp.macro.lp + comp.micro.prok + mb.moss.prok + prey + mb.leaf.euk
mb.leaf.euk ~  comp.macro.le + comp.micro.euk  +  mb.moss.euk + prey + mb.leaf.prok'


fit.comp.w <- sem(comp.w, data = smpMeta.df, estimator = "MLM")
summary(fit.comp.w, rsq = TRUE, fit.measures = TRUE)


resid(fit.comp.w, "cor") # misfit of the bivariate associations
modindices(fit.comp.w, minimum.value = 3, op = "~")


# ## Without linking prey to macrohabitat
# comp.wo <-
# '#prey ~ comp.macro.prey + comp.micro.prey
# mb.leaf.prok ~ comp.macro.lp + comp.micro.prok + mb.moss.prok + prey
# mb.leaf.euk ~  comp.macro.le + comp.micro.euk  +  mb.moss.euk + prey
#
# mb.leaf.prok ~~ mb.leaf.euk
# '
#
#
# fit.comp.wo <- sem(comp.wo, data = smpMeta.df, estimator = "MLM")
# summary(fit.comp.wo, rsq = TRUE, fit.measures = TRUE)
#
#
# resid(fit.comp.wo, "cor") # misfit of the bivariate associations
# modindices(fit.comp.wo, minimum.value = 3, op = "~")


################################################################################
### Plot model
labs_mdl <- c("Prey", "Procaryotes", "Eucaryotes",
              c("Macrohabitat", "Microhabitat"),
              rep(c("Macrohabitat", "Microhabitat", "Moss composition"), 2))



pdf(file = "img/SMP_SEM.pdf", height = 8.27, width = 6)
semPaths(fit.comp.w,
         what = "est", whatLabels = "est",
         residuals = FALSE, intercepts = FALSE,
         sizeMan = 5, sizeMan2 = 3,
         fade = FALSE, layout = "tree", style = "mx",
         nCharNodes = 0, nodeLabels = labs_mdl,
         label.cex = 1, edge.label.cex = 0.35,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",

         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


### Plot composites
labs_micro <- c("Microhabitat prok", "Microhabitat euk", "Microhabitat prey",
                "Habitat size", "Habitat age", "Sample color", "pH", "Prey items")

pdf(file = "img/SMP_composite_micro.pdf", height = 8.27, width = 6)
semPaths(fit.micro,
         what = "est", whatLabels = "est",
         residuals = FALSE, intercepts = FALSE,
         sizeMan = 5, sizeMan2 = 3,
         fade = FALSE, layout = "tree", style = "mx",
         nCharNodes = 0, nodeLabels = labs_micro,
         label.cex = 1, edge.label.cex = 0.35,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()


labs_macro <- c("Macrohabitat prok", "Macrohabitat euk", "Macrohabitat prey",
                "Distance", "Mean temperature", "Altitude", "Canopy cover")

pdf(file = "img/SMP_composite_macro.pdf", height = 8.27, width = 6)
semPaths(fit.macro,
         what = "est", whatLabels = "est",
         residuals = FALSE, intercepts = FALSE,
         sizeMan = 5, sizeMan2 = 3,
         fade = FALSE, layout = "tree", style = "mx",
         nCharNodes = 0,
         nodeLabels = labs_macro,
         label.cex = 1, edge.label.cex = 0.35,
         posCol = "#009e73ff", negCol = "#d55e00ff", edge.label.color = "black",
         layoutSplit = TRUE, curve = 1, curvature = 1, fixedStyle = 1,
         exoCov = FALSE, rotation = 1)
dev.off()
