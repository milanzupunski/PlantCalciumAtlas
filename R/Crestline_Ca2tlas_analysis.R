################################################################################
################################################################################
##
##  Ca2+tlas / Crestline: full statistical and multivariate analysis
##
##  Manuscript: "A calcium-fingerprint atlas for inferring the signaling modes
##               of complex elicitors, applied to a humic biostimulant"
##               Zupunski et al. 2026
##
##  Input:  ONE file, the Supplemental Datasets workbook (sheets "dataset 1"
##          to "dataset 10"). Nothing else is read.
##  Output: every table, statistic and plot is written to the folder set in
##          `out_dir` (Block 0). Every object also stays in your workspace.
##
################################################################################
##
##  WHAT THIS SCRIPT DOES, IN ORDER
##  -------------------------------
##  The Ca2+ response of each root is summarized by six descriptors measured by
##  hand on its crestplot (Fig. 4):
##      sot  = signal onset time (s)
##      dfrt = distance from the root tip at first detection (um)
##      sss  = signal spatial spread (um)
##      sd   = signal duration (s)
##      vt   = tipward velocity (um/s)
##      vs   = shootward velocity (um/s)
##
##  1. Univariate comparison of the 19 reference elicitors (Fig. S4):
##     Kruskal-Wallis per descriptor, Dunn's post hoc test with
##     Benjamini-Hochberg (BH) correction, compact letter display.
##  2. The Ca2+tlas signaling space (Fig. 5a, 6b): PCA of the six descriptors of
##     the 561 reference roots (centered and scaled). The reference roots are
##     "active" individuals and define the axes; everything else (BC4, mutants,
##     GCaMP3, Olympus roots) is projected in afterwards as "supplementary"
##     individuals and never changes the axes.
##  3. Within-treatment descriptor coupling (Fig. 5b): Pearson correlations.
##  4. Class models (Fig. 6c, S6, Tables S2-S4): each of the six a priori
##     classes (MAMP, DAMP, signaling, nutrients/salts, osmotic stress,
##     temperature stress) is modelled as a multivariate normal (MVN)
##     distribution in the first K = 4 principal components, with a shrinkage
##     covariance (corpcor::cov.shrink). For every BC4 root we compute:
##       - Euclidean and squared Mahalanobis distance to every class centroid
##       - posterior class membership (Bayes' rule) under two priors:
##           frequency priors (proportional to class size) and uniform priors
##       - typicality percentile (how typical the root would be as a member)
##       - kNN share (which classes and elicitors its 23 nearest roots belong to)
##     The classifier is validated by stratified 5-fold cross-validation
##     (Table S2, Fig. S6d) and by leave-one-elicitor-out (Table S4).
##  5. Robustness (Fig. S5): GCaMP3 vs R-GECO1, Olympus vs Nikon.
##  6. Mutants (Fig. 7, S8, S9, S11, Table S5): reporter-line validation with
##     FLG22 and NaCl, BC4 descriptors by genotype, PCA of BC4 roots, and
##     projection of mutant BC4 roots into the Ca2+tlas space.
##  7. Root growth (Fig. 7e, S10b): linear model with estimated marginal means.
##  8. Intensity traces (Fig. 1, S1, S7b,d, S9c), plotted from the datasets.
##  9. Exploratory tests not reported in the manuscript (clearly labelled).
## 10. A reproduction check: every number quoted in the manuscript next to the
##     number this script produces.
##
##  HOW TO RUN
##  ----------
##  Run Block 0 and Block 1 first. After that each block runs on its own,
##  except where its header says "Needs Block X". The blocks are plain top-level
##  code: no local(), no wrapper functions. Change anything you like.
##
##  To re-run the whole classifier under the alternative class scheme
##  (PEP1 as DAMP; D-sorbitol and cold merged), set class_scheme <- "relabel"
##  in Block 1 and run Blocks 4 to 8 again.
##
################################################################################
################################################################################



################################################################################
### BLOCK 0. Packages, file locations, output folder
################################################################################

## Install any missing package first (runs only for packages you do not have)
needed  <- c("readxl", "dplyr", "tidyr", "ggplot2", "FactoMineR", "factoextra", "corpcor",
             "mvtnorm", "FSA", "multcompView", "Hmisc", "emmeans", "car")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(readxl)        # read the dataset workbook
library(dplyr)         # data handling
library(tidyr)         # long/wide reshaping
library(ggplot2)       # all plots
library(FactoMineR)    # PCA with supplementary individuals
library(factoextra)    # PCA plots
library(corpcor)       # shrinkage covariance (cov.shrink)
library(mvtnorm)       # multivariate normal density (dmvnorm)
library(FSA)           # Dunn's post hoc test (dunnTest)
library(multcompView)  # compact letter displays
library(Hmisc)         # rcorr (Fig. 5b)
library(emmeans)       # estimated marginal means (Fig. 7e, S10b)
library(car)           # Type II ANOVA (Fig. 7e)

## ---- EDIT THESE TWO LINES ----
setwd("path/to/your/folder")                     # folder that holds the workbook
data_file <- "Supplemental_Datasets.xlsx"         # the dataset workbook (rename yours, or change this line)

out_dir <- "Crestline_outputs"                    # all results go here
dir.create(out_dir, showWarnings = FALSE)

## Descriptor columns, in the order used everywhere below
desc <- c("sot", "dfrt", "sss", "sd", "vt", "vs")
desc_label <- c(sot  = "Signal onset time (s)",
                dfrt = "Distance from the root tip (um)",
                sss  = "Signal spatial spread (um)",
                sd   = "Signal duration (s)",
                vt   = "Tipward velocity (um/s)",
                vs   = "Shootward velocity (um/s)")

## Published values, collected here and compared in Block 18
published <- list()



################################################################################
### BLOCK 1. Load all ten datasets and harmonize names
################################################################################
## Each sheet is read as it is deposited. Column names are shortened so they can
## be typed; the values themselves are not changed. "NA" cells are read as
## missing values (non-responding roots in Dataset 4).

## ---- Dataset 1: Ca2+tlas, 561 reference roots, 19 elicitors (Fig. 5, 6, S4, S5d,e, S6, S11; Tables S2-S4)
d1 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 1", na = c("", "NA")))
names(d1) <- c("genotype", "biosensor", "treatment", "stress_class", "class", "receptor", desc)

## ---- Dataset 2: eATP and L-Glu with R-GECO1 (same roots as Dataset 1) and GCaMP3 (Fig. S5a-c)
d2 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 2", na = c("", "NA")))
names(d2) <- c("genotype", "treatment", "biosensor", "stress_class", "receptor", desc)

## ---- Dataset 3: BC4 in wild type, 20 roots (Fig. 6, S6; same roots as Col-0 BC4 in Dataset 4)
d3 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 3", na = c("", "NA")))
names(d3) <- c("genotype", "biosensor", "treatment", "stress_class", "receptor", desc)

## ---- Dataset 4: Col-0, moca1 and cerk1-2 reporter lines; FLG22, NaCl, BC4 (Fig. 7, S5d,e, S8, S9, S11; Table S5)
d4 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 4", na = c("", "NA")))
names(d4) <- c("genotype", "line", "biosensor", "treatment", "stress_class", "receptor", desc)
d4$line[d4$genotype == "Col-0"] <- "Col-0"                 # wild type has no line number
d4$responder <- !is.na(d4$sot)                              # NA descriptors = no parameterizable response
d4$group <- ifelse(d4$genotype == "Col-0", "Col-0", paste(d4$genotype, d4$line))
d4$group <- factor(d4$group, levels = c("Col-0", "moca1 #1-1", "moca1 #2-1",
                                        "cerk1-2 #1-2", "cerk1-2 #1-3"))

## ---- Dataset 5: raw R-GECO1 intensities, four ROIs, eATP and chitin (Fig. 1)
d5 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 5"))

## ---- Dataset 6: raw R-GECO1 intensities, centerline vs edge ROI (Fig. S1)
d6 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 6"))

## ---- Datasets 7-9: mean +/- SD of dF/F0 traces (Fig. S7b NaCl, S7d FLG22, S9c BC4/mock)
d7 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 7", na = c("", "NA")))
d8 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 8", na = c("", "NA")))
d9 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 9", na = c("", "NA")))
d8$line_group <- gsub("cerk1-2-2", "cerk1-2", d8$line_group)   # typo in 2 rows of the deposited sheet
d9$line_group <- gsub("cerk1-2-2", "cerk1-2", d9$line_group)   # typo in 1 row of the deposited sheet

## ---- Dataset 10: primary root length (cm) at 10 days (Fig. 7e, S10b)
d10 <- as.data.frame(read_xlsx(data_file, sheet = "dataset 10"))
names(d10) <- c("genotype", "line", "treatment", "batch", "length")
d10$batch     <- factor(d10$batch)
d10$treatment <- factor(d10$treatment, levels = c("Control", "Mock", "BC4"))
d10$genotype  <- factor(d10$genotype, levels = c("Col-0", "cerk1-2", "moca1"))

## ---- Elicitor order (by class, as in Fig. S4) ----
elicitor_order <- c("NLP20", "FLG22", "C8", "PG3",                     # MAMP
                    "D-Cellobiose", "eATP", "L-Glu",                   # DAMP
                    "PEP1", "RALF23", "IAA",                           # signaling
                    "NaCl", "KCl", "CaNO32", "KH2PO4", "KNO3",
                    "NH4NO3", "NH42HPO4",                              # nutrients/salts
                    "D-Sorbitol",                                      # osmotic stress
                    "Cold")                                            # temperature stress
d1$treatment <- factor(d1$treatment, levels = elicitor_order)
stopifnot(!anyNA(d1$treatment))                                        # stops if a name does not match

## ---- Class scheme used by the classifier (Blocks 4-8, 14) ----
## "published": the a priori scheme of the manuscript (Methods).
## "relabel"  : sensitivity check, PEP1 moved to DAMP and D-sorbitol + cold
##              merged into one "Osmotic/cold stress" class.
class_scheme <- "published"

d1$class_model <- as.character(d1$class)
if (class_scheme == "relabel") {
  d1$class_model[d1$treatment == "PEP1"] <- "DAMP"
  d1$class_model[d1$treatment %in% c("D-Sorbitol", "Cold")] <- "Osmotic/cold stress"
}
d1$class_model <- factor(d1$class_model)

## ---- Check the input ----
cat("\nDataset 1, roots per elicitor and class:\n")
print(table(d1$treatment, d1$class_model))
cat("\nDataset 4, imaged roots per group and treatment (all roots):\n")
print(with(d4, table(group, treatment)))
cat("\nDataset 4, responding roots only:\n")
print(with(d4[d4$responder, ], table(group, treatment)))



################################################################################
### BLOCK 2. Fig. S4: descriptors of the 19 reference elicitors
################################################################################
## For each descriptor: Kruskal-Wallis test across elicitors, then Dunn's test
## (all pairs, BH-adjusted). Letters: groups sharing a letter do not differ
## (adjusted p >= 0.05); 'a' goes to the group with the highest median.
## Velocities are plotted on a log10 axis for display only; tests use raw values.

published$S4_H <- c(sot = 481.0, dfrt = 349.9, sd = 384.0, sss = 333.8, vt = 461.5, vs = 466.2)

s4_kw <- data.frame()
s4_letters <- list()

for (d in desc) {
  ## Kruskal-Wallis
  kw <- kruskal.test(d1[[d]] ~ d1$treatment)
  s4_kw <- rbind(s4_kw, data.frame(descriptor = d, H = round(unname(kw$statistic), 1),
                                   df = unname(kw$parameter), p = signif(kw$p.value, 2)))

  ## Dunn's test, BH-adjusted
  dunn_res <- dunnTest(as.formula(paste(d, "~ treatment")), data = d1, method = "bh")$res

  ## Symmetric matrix of adjusted p-values, groups ordered by decreasing median
  grp_order <- names(sort(tapply(d1[[d]], d1$treatment, median), decreasing = TRUE))
  p_mat <- matrix(1, length(grp_order), length(grp_order), dimnames = list(grp_order, grp_order))
  pairs <- strsplit(as.character(dunn_res$Comparison), " - ")
  for (i in seq_along(pairs)) {
    p_mat[pairs[[i]][1], pairs[[i]][2]] <- dunn_res$P.adj[i]
    p_mat[pairs[[i]][2], pairs[[i]][1]] <- dunn_res$P.adj[i]
  }
  s4_letters[[d]] <- multcompLetters(p_mat, compare = "<", threshold = 0.05)$Letters

  ## Plot
  lab <- data.frame(treatment = names(s4_letters[[d]]), letter = s4_letters[[d]],
                    y = max(d1[[d]]) * ifelse(d %in% c("vt", "vs"), 2, 1.08))
  p <- ggplot(d1, aes(treatment, .data[[d]], colour = treatment)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(alpha = 0.4, position = position_jitter(width = 0.15, seed = 1)) +
    geom_text(data = lab, aes(treatment, y, label = letter), colour = "black", size = 3) +
    scale_colour_viridis_d(option = "D") +
    labs(x = NULL, y = desc_label[d]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
          legend.position = "none")
  if (d %in% c("vt", "vs")) p <- p + scale_y_log10()
  print(p)
  ggsave(file.path(out_dir, paste0("FigS4_", d, ".png")), p, width = 8, height = 4.5, dpi = 300, bg = "white")
}

s4_kw
s4_letters
write.csv(s4_kw, file.path(out_dir, "FigS4_KruskalWallis.csv"), row.names = FALSE)
write.csv(as.data.frame(s4_letters), file.path(out_dir, "FigS4_letters.csv"))



################################################################################
### BLOCK 3. Fig. 5a,b: the Ca2+tlas PCA and within-treatment correlations
################################################################################
## PCA on the six untransformed descriptors, centered and scaled (each
## descriptor gets unit variance, so onset time in s and spread in um weigh
## equally). Only the 561 reference roots.

published$PC_var <- c(PC1 = 46.0, PC2 = 17.2, cum_PC4 = 87.4)

pca_atlas <- PCA(d1[, desc], scale.unit = TRUE, ncp = 6, graph = FALSE)
pca_atlas$eig                  # variance explained per component
pca_atlas$var$coord            # loadings (correlation of each descriptor with each PC)
pca_atlas$var$cos2             # quality of representation of each descriptor

p5a <- fviz_pca(pca_atlas, habillage = d1$treatment, label = "var", col.var = "black",
                alpha.var = "cos2", alpha.ind = 0, addEllipses = TRUE,
                ellipse.type = "confidence", repel = TRUE, title = "") +
  scale_colour_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic()
print(p5a)
ggsave(file.path(out_dir, "Fig5a_PCA_biplot.png"), p5a, width = 8, height = 6, dpi = 300, bg = "white")
write.csv(pca_atlas$eig, file.path(out_dir, "Fig5a_PCA_eigenvalues.csv"))
write.csv(pca_atlas$var$coord, file.path(out_dir, "Fig5a_PCA_loadings.csv"))

## ---- Fig. 5b: within-treatment Pearson correlations (Hmisc::rcorr) ----
## For each elicitor, the 6 x 6 correlation matrix across its roots. p-values
## are BH-adjusted within each elicitor; the mosaic is descriptive.
cor_long <- data.frame()
for (tr in levels(d1$treatment)) {
  rc <- rcorr(as.matrix(d1[d1$treatment == tr, desc]), type = "pearson")
  ut <- which(upper.tri(rc$r), arr.ind = TRUE)
  tmp <- data.frame(treatment = tr, var1 = desc[ut[, 1]], var2 = desc[ut[, 2]],
                    r = rc$r[ut], p = rc$P[ut])
  tmp$p_adj_BH <- p.adjust(tmp$p, method = "BH")
  cor_long <- rbind(cor_long, tmp)
}
cor_long$treatment <- factor(cor_long$treatment, levels = levels(d1$treatment))

## Mosaic layout: each descriptor pair is one cell; inside it, one small tile per
## elicitor (5 columns x 4 rows).
cor_long$i  <- match(cor_long$var1, desc)
cor_long$j  <- match(cor_long$var2, desc)
cor_long$k  <- as.integer(cor_long$treatment) - 1
cor_long$x  <- cor_long$j - 0.5 + (cor_long$k %% 5 + 0.5) / 5
cor_long$y  <- cor_long$i - 0.5 + (3 - cor_long$k %/% 5 + 0.5) / 4

## Pairs whose coupling differs most across elicitors (largest SD of r)
pair_sd <- cor_long %>% group_by(var1, var2, i, j) %>%
  summarise(sd_r = sd(r, na.rm = TRUE), .groups = "drop") %>% arrange(desc(sd_r))
top_pairs <- head(pair_sd, 3)
pair_sd

p5b <- ggplot(cor_long, aes(x, y, fill = r)) +
  geom_tile(width = 1 / 5, height = 1 / 4, colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
  scale_x_continuous(breaks = 1:6, labels = desc, position = "top") +
  scale_y_reverse(breaks = 1:6, labels = desc) +
  coord_equal() + labs(x = NULL, y = NULL, fill = "Pearson r") + theme_minimal()
print(p5b)
ggsave(file.path(out_dir, "Fig5b_correlation_mosaic.png"), p5b, width = 7, height = 6, dpi = 300, bg = "white")
write.csv(cor_long[, c("treatment", "var1", "var2", "r", "p", "p_adj_BH")],
          file.path(out_dir, "Fig5b_within_treatment_correlations.csv"), row.names = FALSE)



################################################################################
### BLOCK 4. Fig. 6b: BC4 projected into the Ca2+tlas space; choice of K
################################################################################
## The 20 BC4 roots (Dataset 3) are supplementary individuals: they get
## coordinates in the atlas PCA but do not shape it. The coordinates of the
## 561 reference roots are identical to Block 3.

X_all   <- rbind(d1[, desc], d3[, desc])
ind_act <- seq_len(nrow(d1))                             # reference roots
ind_sup <- nrow(d1) + seq_len(nrow(d3))                  # BC4 roots

pca <- PCA(X_all, scale.unit = TRUE, ncp = 6, graph = FALSE, ind.sup = ind_sup)

## Number of PCs for inference: smallest K reaching >= 80 % variance
## (capped at min class size - 2 and at 10, as in the original analysis)
group  <- droplevels(d1$class_model)
classes <- levels(group)
cumvar <- cumsum(pca$eig[, "percentage of variance"]) / 100
K <- which(cumvar >= 0.80)[1]
K <- max(2, min(K, max(2, min(table(group)) - 2), 10))
cat(sprintf("\nK = %d PCs, %.1f %% of variance\n", K, 100 * cumvar[K]))

X_act <- as.data.frame(pca$ind$coord[, 1:K])             # reference roots in PC space
X_sup <- as.data.frame(pca$ind.sup$coord[, 1:K])         # BC4 roots in PC space
rownames(X_sup) <- paste0("BC4_", seq_len(nrow(X_sup)))

## Fig. 6b. Note: ellipse.type = "t" draws 95 % data ellipses (multivariate t),
## as in the original figure script. For centroid confidence ellipses use
## ellipse.type = "confidence".
p6b <- fviz_pca_ind(pca, habillage = group, invisible = "ind.sup", label = "none",
                    pointsize = 2.5, alpha.ind = 0, addEllipses = TRUE,
                    ellipse.type = "t", legend.title = "Elicitor class", title = "") +
  geom_point(data = as.data.frame(pca$ind.sup$coord), aes(Dim.1, Dim.2),
             inherit.aes = FALSE, colour = "magenta", shape = 8, size = 3) +
  scale_colour_viridis_d(option = "E") + scale_fill_viridis_d(option = "E") +
  theme_classic() +
  scale_x_reverse() + scale_y_reverse()                  # axis orientation as in the figure
print(p6b)
ggsave(file.path(out_dir, "Fig6b_PCA_BC4_supplementary.png"), p6b, width = 7, height = 6, dpi = 300, bg = "white")
write.csv(pca$ind.sup$coord, file.path(out_dir, "Fig6b_BC4_PC_coordinates.csv"))



################################################################################
### BLOCK 5. Class models and BC4 class assignment (Fig. 6c, S6b,c; Table S3)
### Needs Block 4
################################################################################
## One MVN model per class, in the first K PCs:
##   centroid   = mean of the class's reference roots
##   covariance = shrinkage estimate (cov.shrink), symmetrized
## Posterior for class c:  P(c | x) ~ prior(c) * MVN density(x | centroid_c, cov_c)

centroids <- lapply(classes, function(cl) colMeans(X_act[group == cl, ]))
names(centroids) <- classes
covs <- lapply(classes, function(cl) {
  S <- as.matrix(unclass(cov.shrink(as.matrix(X_act[group == cl, ]), verbose = FALSE)))
  (S + t(S)) / 2
})
names(covs) <- classes

priors_freq <- as.numeric(prop.table(table(group)));  names(priors_freq) <- classes
priors_unif <- rep(1 / length(classes), length(classes)); names(priors_unif) <- classes

## ---- 1. Euclidean distance to each class centroid ----
euclid_bc4 <- sapply(classes, function(cl)
  sqrt(rowSums(sweep(as.matrix(X_sup), 2, centroids[[cl]])^2)))

## ---- 2. Squared Mahalanobis distance (prior-free) ----
mahal_bc4 <- sapply(classes, function(cl) mahalanobis(as.matrix(X_sup), centroids[[cl]], covs[[cl]]))

## ---- 3. Posterior membership, frequency and uniform priors ----
loglik_bc4 <- sapply(classes, function(cl)
  dmvnorm(as.matrix(X_sup), mean = centroids[[cl]], sigma = covs[[cl]], log = TRUE))
post_freq <- exp(sweep(loglik_bc4, 2, log(priors_freq), "+"))
post_freq <- post_freq / rowSums(post_freq)
post_unif <- exp(sweep(loglik_bc4, 2, log(priors_unif), "+"))
post_unif <- post_unif / rowSums(post_unif)
rownames(euclid_bc4) <- rownames(mahal_bc4) <- rownames(post_freq) <- rownames(post_unif) <- rownames(X_sup)

## ---- 4. Typicality percentile: share of the class's own roots closer to its centroid ----
pctl_bc4 <- sapply(classes, function(cl) {
  d_ref <- mahalanobis(as.matrix(X_act[group == cl, ]), centroids[[cl]], covs[[cl]])
  sapply(mahal_bc4[, cl], function(d) mean(d_ref <= d))
})

## ---- 5. kNN: classes and elicitors among the k nearest reference roots ----
k_nn  <- min(max(5, floor(sqrt(nrow(X_act)))), nrow(X_act) - 1)      # 23
d_bc4 <- as.matrix(dist(rbind(X_act, X_sup)))[ind_sup, ind_act]
knn_class <- t(apply(d_bc4, 1, function(dd)
  prop.table(table(factor(group[order(dd)[1:k_nn]], levels = classes)))))
knn_elicitor <- t(apply(d_bc4, 1, function(dd)
  prop.table(table(factor(d1$treatment[order(dd)[1:k_nn]], levels = levels(d1$treatment))))))
rownames(knn_class) <- rownames(knn_elicitor) <- rownames(X_sup)

## ---- Per root, per class (long table) and per root (summary) ----
bc4_long <- data.frame(sup_id = rep(rownames(X_sup), times = length(classes)),
                       class = rep(classes, each = nrow(X_sup)),
                       euclid_dist = as.vector(euclid_bc4), mahal_dist2 = as.vector(mahal_bc4),
                       posterior_freq = as.vector(post_freq), posterior_unif = as.vector(post_unif),
                       typicality_pctl = as.vector(pctl_bc4), knn_share = as.vector(knn_class))

bc4_summary <- bc4_long %>% group_by(sup_id) %>%
  summarise(best_class_freq  = class[which.max(posterior_freq)],
            best_post_freq   = max(posterior_freq),
            margin_freq      = best_post_freq - sort(posterior_freq, decreasing = TRUE)[2],
            best_class_unif  = class[which.max(posterior_unif)],
            best_post_unif   = max(posterior_unif),
            nearest_mahal    = class[which.min(mahal_dist2)],
            nearest_euclid   = class[which.min(euclid_dist)],
            best_class_knn   = class[which.max(knn_share)],
            .groups = "drop")

## ---- Table S3 ----
tableS3 <- data.frame(
  class               = classes,
  prior_freq          = round(priors_freq, 3),
  mean_post_freq      = round(colMeans(post_freq), 3),
  mean_post_unif      = round(colMeans(post_unif), 3),
  assigned_freq       = as.integer(table(factor(classes[max.col(post_freq, "first")], levels = classes))),
  assigned_unif       = as.integer(table(factor(classes[max.col(post_unif, "first")], levels = classes))),
  median_percentile   = round(apply(pctl_bc4, 2, median), 2),
  within_central95    = colSums(pctl_bc4 <= 0.95),
  nearest_mahalanobis = as.integer(table(factor(classes[max.col(-mahal_bc4, "first")], levels = classes))),
  mean_knn_share      = round(colMeans(knn_class), 2),
  row.names = NULL)
tableS3 <- tableS3[order(-tableS3$mean_post_freq), ]
tableS3
round(sort(colMeans(knn_elicitor), decreasing = TRUE), 2)   # nearest reference elicitors of BC4

write.csv(tableS3, file.path(out_dir, "TableS3_BC4_class_assignment.csv"), row.names = FALSE)
write.csv(bc4_long, file.path(out_dir, "BC4_per_root_per_class.csv"), row.names = FALSE)
write.csv(bc4_summary, file.path(out_dir, "BC4_per_root_summary.csv"), row.names = FALSE)
write.csv(round(knn_elicitor, 2), file.path(out_dir, "BC4_knn_by_elicitor.csv"))

## ---- Plots: Fig. S6b (stacked bars), S6c (heatmap), Fig. 6c (cumulative
##      posterior mass), total mass and top-1 wins. Frequency priors as in the
##      figures; set `post_plot <- post_unif` to draw the uniform-prior versions.
post_plot <- post_freq
prior_tag <- "freq"

ord <- rownames(post_plot)[order(-apply(post_plot, 1, max),
                                 -(apply(post_plot, 1, max) -
                                   apply(post_plot, 1, function(v) sort(v, decreasing = TRUE)[2])))]
post_plot_long <- data.frame(sup_id = factor(rep(rownames(post_plot), times = length(classes)), levels = ord),
                             class = factor(rep(classes, each = nrow(post_plot)), levels = classes),
                             posterior = as.vector(post_plot))

pS6b <- ggplot(post_plot_long, aes(sup_id, posterior, fill = class)) +
  geom_col(width = 0.9) +
  scale_y_continuous(limits = c(0, 1.0001), expand = c(0, 0)) +
  scale_fill_viridis_d(option = "E") +
  labs(x = "BC4 root, ordered by best posterior", y = "Posterior membership probability", fill = "Class") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(pS6b)
ggsave(file.path(out_dir, paste0("FigS6b_stacked_posteriors_", prior_tag, ".png")), pS6b, width = 10, height = 4.5, dpi = 300, bg = "white")

pS6c <- ggplot(post_plot_long, aes(class, sup_id, fill = posterior)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "D", limits = c(0, 1)) +
  labs(x = "Class", y = "BC4 root", fill = "Posterior") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(pS6c)
ggsave(file.path(out_dir, paste0("FigS6c_posterior_heatmap_", prior_tag, ".png")), pS6c, width = 7, height = 5, dpi = 300, bg = "white")

cum_df <- post_plot_long %>% arrange(sup_id) %>% group_by(class) %>%
  mutate(idx = as.integer(sup_id), cum_posterior = cumsum(posterior)) %>% ungroup()
p6c <- ggplot(cum_df, aes(idx, cum_posterior, colour = class)) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(option = "E") +
  labs(x = "BC4 roots (ordered by best posterior)", y = "Cumulative posterior mass", colour = "Class") +
  theme_classic() + theme(aspect.ratio = 1)
print(p6c)
ggsave(file.path(out_dir, paste0("Fig6c_cumulative_posterior_", prior_tag, ".png")), p6c, width = 6, height = 5, dpi = 300, bg = "white")

mass_df <- post_plot_long %>% group_by(class) %>% summarise(total_posterior = sum(posterior), .groups = "drop")
p_mass <- ggplot(mass_df, aes(reorder(class, total_posterior), total_posterior, fill = class)) +
  geom_col() + coord_flip() + scale_fill_viridis_d(option = "E") +
  labs(x = NULL, y = "Total posterior mass across BC4 roots") +
  theme_classic() + theme(legend.position = "none")
print(p_mass)
ggsave(file.path(out_dir, paste0("BC4_total_posterior_mass_", prior_tag, ".png")), p_mass, width = 6, height = 4, dpi = 300, bg = "white")

top1_df <- data.frame(idx = seq_along(ord), class = classes[max.col(post_plot[ord, ], "first")])
top1_step <- expand.grid(idx = seq_along(ord), class = classes, stringsAsFactors = FALSE)
top1_step$win <- as.integer(paste(top1_step$idx, top1_step$class) %in% paste(top1_df$idx, top1_df$class))
top1_step <- top1_step %>% group_by(class) %>% arrange(idx) %>% mutate(cum_wins = cumsum(win)) %>% ungroup()
p_top1 <- ggplot(top1_step, aes(idx, cum_wins, colour = class)) +
  geom_step(linewidth = 1) + scale_colour_viridis_d(option = "E") +
  labs(x = "BC4 roots added (ordered by best posterior)", y = "Cumulative top-1 assignments", colour = "Class") +
  theme_classic()
print(p_top1)
ggsave(file.path(out_dir, paste0("BC4_top1_wins_", prior_tag, ".png")), p_top1, width = 7, height = 4, dpi = 300, bg = "white")



################################################################################
### BLOCK 6. Stratified 5-fold cross-validation (Fig. S6d, Table S2)
### Needs Block 4
################################################################################
## Every reference root is predicted once by a model fitted without it. Folds
## are drawn within each class so each fold keeps the class proportions.
## Class models and priors are refitted on the training folds only. The PCA
## itself is not refitted (see Block 7 for a test that refits it).

published$CV <- c(accuracy = 0.743, balanced_accuracy = 0.676)
published$recall <- c(MAMP = 0.903, Signaling = 0.797, "Nutrient/salt" = 0.780,
                      "Osmotic stress" = 0.667, "Temperature stress" = 0.545, DAMP = 0.365)

k_folds <- 5
set.seed(123)
fold_id <- integer(length(group))
for (cl in classes) {
  idx <- which(group == cl)
  fold_id[idx] <- sample(rep(1:k_folds, length.out = length(idx)))
}

pred_cv_freq <- character(length(group))
pred_cv_unif <- character(length(group))
for (f in 1:k_folds) {
  tr <- which(fold_id != f)
  te <- which(fold_id == f)
  g_tr <- droplevels(group[tr])
  cl_tr <- levels(g_tr)
  mu_f <- lapply(cl_tr, function(cl) colMeans(X_act[tr, ][g_tr == cl, ]))
  S_f  <- lapply(cl_tr, function(cl) {
    S <- as.matrix(unclass(cov.shrink(as.matrix(X_act[tr, ][g_tr == cl, ]), verbose = FALSE)))
    (S + t(S)) / 2
  })
  ll <- sapply(seq_along(cl_tr), function(j) dmvnorm(as.matrix(X_act[te, ]), mu_f[[j]], S_f[[j]], log = TRUE))
  ll <- matrix(ll, nrow = length(te))
  lp_freq <- sweep(ll, 2, log(as.numeric(prop.table(table(g_tr)))), "+")
  lp_unif <- sweep(ll, 2, log(rep(1 / length(cl_tr), length(cl_tr))), "+")
  pred_cv_freq[te] <- cl_tr[max.col(lp_freq, "first")]
  pred_cv_unif[te] <- cl_tr[max.col(lp_unif, "first")]
}

cm_freq <- table(True = group, Pred = factor(pred_cv_freq, levels = classes))
cm_unif <- table(True = group, Pred = factor(pred_cv_unif, levels = classes))
cm_freq

cv_summary <- data.frame(
  prior             = c("frequency", "uniform"),
  accuracy          = round(c(sum(diag(cm_freq)), sum(diag(cm_unif))) / length(group), 3),
  balanced_accuracy = round(c(mean(diag(cm_freq) / rowSums(cm_freq)), mean(diag(cm_unif) / rowSums(cm_unif))), 3),
  majority_baseline = round(max(prop.table(table(group))), 3),
  chance_balanced   = round(1 / length(classes), 3))
cv_summary

tableS2 <- data.frame(class = classes, n = as.integer(rowSums(cm_freq)),
                      recall_freq    = round(diag(cm_freq) / rowSums(cm_freq), 3),
                      precision_freq = round(diag(cm_freq) / colSums(cm_freq), 3),
                      recall_unif    = round(diag(cm_unif) / rowSums(cm_unif), 3),
                      row.names = NULL)
tableS2 <- tableS2[order(-tableS2$recall_freq), ]
tableS2

cm_df <- as.data.frame(cm_freq) %>% group_by(True) %>% mutate(row_fraction = Freq / sum(Freq)) %>% ungroup()
pS6d <- ggplot(cm_df, aes(Pred, True, fill = row_fraction)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_viridis_c(option = "D", limits = c(0, 1), name = "Row fraction") +
  labs(x = "Predicted class", y = "True class",
       title = sprintf("Stratified 5-fold CV: acc = %.2f, bal.acc = %.2f",
                       cv_summary$accuracy[1], cv_summary$balanced_accuracy[1])) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(pS6d)
ggsave(file.path(out_dir, "FigS6d_confusion_matrix.png"), pS6d, width = 7, height = 6, dpi = 300, bg = "white")
write.csv(tableS2, file.path(out_dir, "TableS2_per_class_recall.csv"), row.names = FALSE)
write.csv(cv_summary, file.path(out_dir, "TableS2_CV_summary.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(cm_freq), file.path(out_dir, "FigS6d_confusion_counts_freq.csv"))



################################################################################
### BLOCK 7. Leave-one-elicitor-out validation (Table S4)
### Needs Block 4
################################################################################
## Harder than Block 6: all roots of one elicitor are withheld, the PCA AND the
## class models are refitted on the other 18 elicitors, and the withheld
## elicitor is projected and classified. This asks whether an elicitor the atlas
## has never seen lands in its expected class. Classes with only one elicitor
## (osmotic, temperature) cannot be tested this way.

published$LOEO <- c(correct_freq = 12, correct_unif = 11, evaluable = 17,
                    root_acc_freq = 0.58, root_acc_unif = 0.56)

elic <- droplevels(d1$treatment)
loeo <- data.frame()

for (el in levels(elic)) {
  te <- which(elic == el)
  tr <- which(elic != el)
  true_cl <- as.character(group[te][1])
  g_tr <- droplevels(group[tr])

  if (!(true_cl %in% levels(g_tr))) {
    loeo <- rbind(loeo, data.frame(elicitor = el, true_class = true_cl, n = length(te), prior = c("frequency", "uniform"),
                                   predicted = "not evaluable", frac_correct = NA, mean_post_true = NA))
    next
  }

  pca_el <- PCA(d1[, desc], scale.unit = TRUE, ncp = K, ind.sup = te, graph = FALSE)
  Xtr <- pca_el$ind$coord[, 1:K]
  Xte <- pca_el$ind.sup$coord[, 1:K, drop = FALSE]
  cl_tr <- levels(g_tr)
  mu_e <- lapply(cl_tr, function(cl) colMeans(Xtr[g_tr == cl, ]))
  S_e  <- lapply(cl_tr, function(cl) {
    S <- as.matrix(unclass(cov.shrink(Xtr[g_tr == cl, ], verbose = FALSE)))
    (S + t(S)) / 2
  })
  ll <- matrix(sapply(seq_along(cl_tr), function(j) dmvnorm(Xte, mu_e[[j]], S_e[[j]], log = TRUE)),
               nrow = length(te))

  for (pr in c("frequency", "uniform")) {
    prior_vec <- if (pr == "uniform") rep(1 / length(cl_tr), length(cl_tr)) else as.numeric(prop.table(table(g_tr)))
    lp <- sweep(ll, 2, log(prior_vec), "+")
    P  <- exp(lp - apply(lp, 1, max)); P <- P / rowSums(P); colnames(P) <- cl_tr
    pred <- cl_tr[max.col(P, "first")]
    loeo <- rbind(loeo, data.frame(elicitor = el, true_class = true_cl, n = length(te), prior = pr,
                                   predicted = names(which.max(table(pred))),   # plurality vote
                                   frac_correct = round(mean(pred == true_cl), 3),
                                   mean_post_true = round(mean(P[, true_cl]), 3)))
  }
}

tableS4 <- pivot_wider(loeo, id_cols = c(elicitor, true_class, n), names_from = prior,
                       values_from = c(predicted, frac_correct, mean_post_true))
as.data.frame(tableS4)

loeo_summary <- loeo %>% filter(predicted != "not evaluable") %>% group_by(prior) %>%
  summarise(evaluable = n(), correct = sum(predicted == true_class),
            root_accuracy = round(weighted.mean(frac_correct, n), 2), .groups = "drop")
loeo_summary
write.csv(tableS4, file.path(out_dir, "TableS4_LOEO.csv"), row.names = FALSE)
write.csv(loeo_summary, file.path(out_dir, "TableS4_LOEO_summary.csv"), row.names = FALSE)



################################################################################
### BLOCK 8. Fig. S6e: BC4 descriptors with NaCl and eATP reference medians
################################################################################

ref_medians <- d1 %>% filter(treatment %in% c("NaCl", "eATP")) %>% group_by(treatment) %>%
  summarise(across(all_of(desc), median), .groups = "drop")
ref_medians

for (d in desc) {
  p <- ggplot(d3, aes(x = "BC4", y = .data[[d]])) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.12, seed = 1)) +
    geom_hline(data = ref_medians, aes(yintercept = .data[[d]], colour = treatment), linetype = "dashed") +
    scale_colour_manual(values = c(eATP = "#d7301f", NaCl = "#2c7fb8")) +
    expand_limits(y = range(c(d3[[d]], ref_medians[[d]]))) +
    labs(x = NULL, y = desc_label[d], colour = "Reference median") + theme_classic()
  if (d %in% c("vt", "vs")) p <- p + scale_y_log10()
  print(p)
  ggsave(file.path(out_dir, paste0("FigS6e_", d, ".png")), p, width = 3.2, height = 3.6, dpi = 300, bg = "white")
}



################################################################################
### BLOCK 9. Fig. S5c: GCaMP3 roots projected into the R-GECO1 eATP / L-Glu space
################################################################################
## PCA fitted on the 57 R-GECO1 roots (eATP 29, L-Glu 28); the GCaMP3 roots are
## supplementary. Use the % variance printed here on the axes.

published$S5c <- c(PC1 = 49.3, PC2 = 18.2)

s5_sup <- which(d2$biosensor == "GCaMP3")
pca_s5 <- PCA(d2[, desc], scale.unit = TRUE, ncp = 6, graph = FALSE, ind.sup = s5_sup)
pca_s5$eig[1:2, ]
table(d2$treatment, d2$biosensor)

s5_group <- factor(d2$treatment[-s5_sup])
s5_sup_xy <- data.frame(pca_s5$ind.sup$coord[, 1:2], treatment = d2$treatment[s5_sup])
pS5c <- fviz_pca_ind(pca_s5, habillage = s5_group, invisible = "ind.sup", label = "none",
                     addEllipses = TRUE, ellipse.type = "norm", ellipse.level = 0.95,
                     mean.point = TRUE, alpha.ind = 0.3, title = "") +
  geom_point(data = s5_sup_xy, aes(Dim.1, Dim.2, colour = treatment), inherit.aes = FALSE, shape = 4, size = 3) +
  theme_classic()
print(pS5c)
ggsave(file.path(out_dir, "FigS5c_GCaMP3_projection.png"), pS5c, width = 6, height = 5, dpi = 300, bg = "white")



################################################################################
### BLOCK 10. Fig. S5d,e: imaging platform (Nikon atlas vs Olympus mutant sessions)
################################################################################
## Nikon SMZ18: wild-type FLG22 and NaCl roots of the Ca2+tlas (Dataset 1).
## Olympus MVX10: wild-type FLG22 and NaCl roots recorded with the mutants
## (Col-0 rows of Dataset 4). Mann-Whitney per descriptor, BH-adjusted within
## each stimulus (6 tests). Then a PCA on the Nikon roots with the Olympus roots
## projected.

platform <- rbind(
  data.frame(platform = "Nikon",   treatment = as.character(d1$treatment), d1[, desc])[d1$treatment %in% c("FLG22", "NaCl"), ],
  data.frame(platform = "Olympus", treatment = d4$treatment, d4[, desc])[d4$group == "Col-0" & d4$treatment %in% c("FLG22", "NaCl") & d4$responder, ])
table(platform$platform, platform$treatment)

s5d_tests <- data.frame()
for (st in c("FLG22", "NaCl")) {
  for (d in desc) {
    sub <- platform[platform$treatment == st, ]
    wt <- wilcox.test(sub[[d]] ~ sub$platform, exact = FALSE)
    s5d_tests <- rbind(s5d_tests, data.frame(stimulus = st, descriptor = d,
                                             median_Nikon = median(sub[[d]][sub$platform == "Nikon"]),
                                             median_Olympus = median(sub[[d]][sub$platform == "Olympus"]),
                                             W = unname(wt$statistic), p = wt$p.value))
  }
}
s5d_tests$p_adj_BH <- ave(s5d_tests$p, s5d_tests$stimulus, FUN = function(p) p.adjust(p, "BH"))
s5d_tests$ratio_Olympus_Nikon <- round(s5d_tests$median_Olympus / s5d_tests$median_Nikon, 2)
s5d_tests
write.csv(s5d_tests, file.path(out_dir, "FigS5d_platform_MannWhitney.csv"), row.names = FALSE)

s5d_long <- pivot_longer(platform, all_of(desc), names_to = "descriptor", values_to = "value")
s5d_long$descriptor <- factor(s5d_long$descriptor, levels = desc)
pS5d <- ggplot(s5d_long, aes(treatment, value, colour = platform)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, seed = 1), alpha = 0.5) +
  facet_wrap(~ descriptor, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL) + theme_classic()
print(pS5d)
ggsave(file.path(out_dir, "FigS5d_platform_descriptors.png"), pS5d, width = 12, height = 3.5, dpi = 300, bg = "white")

s5e_sup <- which(platform$platform == "Olympus")
pca_s5e <- PCA(platform[, desc], scale.unit = TRUE, ncp = 6, graph = FALSE, ind.sup = s5e_sup)
s5e_df <- rbind(data.frame(pca_s5e$ind$coord[, 1:2],     platform[-s5e_sup, c("platform", "treatment")]),
                data.frame(pca_s5e$ind.sup$coord[, 1:2], platform[s5e_sup,  c("platform", "treatment")]))
pS5e <- ggplot(s5e_df, aes(Dim.1, Dim.2, colour = treatment, shape = platform, linetype = platform)) +
  geom_point() + stat_ellipse(type = "norm", level = 0.95) +
  scale_shape_manual(values = c(Nikon = 16, Olympus = 4)) +
  labs(x = sprintf("Dim1 (%.1f %%)", pca_s5e$eig[1, 2]), y = sprintf("Dim2 (%.1f %%)", pca_s5e$eig[2, 2])) +
  theme_classic()
print(pS5e)
ggsave(file.path(out_dir, "FigS5e_platform_PCA.png"), pS5e, width = 6, height = 5, dpi = 300, bg = "white")



################################################################################
### BLOCK 11. Fig. S8: FLG22 and NaCl responses of the reporter lines
################################################################################
## Response rates (responding/total) with Fisher's exact test across the five
## lines; descriptors of responding roots with Kruskal-Wallis and Dunn (BH).
## Lines with no responding root are left out of the descriptor tests.
## NOTE: response rates are only correct if Dataset 4 lists every imaged root
## (non-responders as NA).

s8_rates <- d4 %>% filter(treatment %in% c("FLG22", "NaCl")) %>% group_by(treatment, group) %>%
  summarise(responding = sum(responder), total = n(), .groups = "drop")
s8_rates

s8_fisher <- data.frame()
for (st in c("FLG22", "NaCl")) {
  sub <- d4[d4$treatment == st, ]
  ft <- fisher.test(table(droplevels(sub$group), factor(sub$responder, levels = c(TRUE, FALSE))))
  s8_fisher <- rbind(s8_fisher, data.frame(stimulus = st, p = signif(ft$p.value, 3)))
}
s8_fisher

s8_kw <- data.frame()
s8_letters <- list()
for (st in c("FLG22", "NaCl")) {
  sub <- d4[d4$treatment == st & d4$responder, ]
  sub$group <- droplevels(sub$group)
  for (d in desc) {
    kw <- kruskal.test(sub[[d]] ~ sub$group)
    s8_kw <- rbind(s8_kw, data.frame(stimulus = st, descriptor = d, H = round(unname(kw$statistic), 2),
                                     df = unname(kw$parameter), p = signif(kw$p.value, 3)))
    dunn_res <- dunnTest(as.formula(paste(d, "~ group")), data = sub, method = "bh")$res
    grp_order <- names(sort(tapply(sub[[d]], sub$group, median), decreasing = TRUE))
    p_mat <- matrix(1, length(grp_order), length(grp_order), dimnames = list(grp_order, grp_order))
    pairs <- strsplit(as.character(dunn_res$Comparison), " - ")
    for (i in seq_along(pairs)) {
      p_mat[pairs[[i]][1], pairs[[i]][2]] <- dunn_res$P.adj[i]
      p_mat[pairs[[i]][2], pairs[[i]][1]] <- dunn_res$P.adj[i]
    }
    s8_letters[[paste(st, d)]] <- multcompLetters(p_mat, compare = "<", threshold = 0.05)$Letters
  }

  s8_long <- pivot_longer(sub, all_of(desc), names_to = "descriptor", values_to = "value")
  s8_long$descriptor <- factor(s8_long$descriptor, levels = desc)
  p <- ggplot(s8_long, aes(group, value, colour = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.15, seed = 1), alpha = 0.6) +
    facet_wrap(~ descriptor, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = NULL, title = st) + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  print(p)
  ggsave(file.path(out_dir, paste0("FigS8_", st, "_descriptors.png")), p, width = 13, height = 4, dpi = 300, bg = "white")
}
s8_kw
s8_letters
write.csv(s8_rates, file.path(out_dir, "FigS8_response_rates.csv"), row.names = FALSE)
write.csv(s8_fisher, file.path(out_dir, "FigS8_Fisher.csv"), row.names = FALSE)
write.csv(s8_kw, file.path(out_dir, "FigS8_KruskalWallis.csv"), row.names = FALSE)



################################################################################
### BLOCK 12. Fig. 7b,c and S9d-g: BC4 descriptors by genotype
################################################################################
## BC4-treated roots: Col-0 (n = 20), moca1 #2-1 (13), cerk1-2 #1-3 (12).
## All roots responded. Kruskal-Wallis + Dunn (BH) per descriptor.

published$Fig7_H <- c(sss = 14.2, sd = 18.7, dfrt = 11.2)

bc4_geno <- d4[d4$treatment == "BC4", ]
bc4_geno$genotype <- factor(bc4_geno$genotype, levels = c("Col-0", "moca1", "cerk1-2"))
table(bc4_geno$genotype, bc4_geno$line)

fig7_kw <- data.frame()
fig7_letters <- list()
for (d in desc) {
  kw <- kruskal.test(bc4_geno[[d]] ~ bc4_geno$genotype)
  fig7_kw <- rbind(fig7_kw, data.frame(descriptor = d, H = round(unname(kw$statistic), 1),
                                       df = unname(kw$parameter), p = signif(kw$p.value, 3)))
  dunn_res <- dunnTest(as.formula(paste(d, "~ genotype")), data = bc4_geno, method = "bh")$res
  grp_order <- names(sort(tapply(bc4_geno[[d]], bc4_geno$genotype, median), decreasing = TRUE))
  p_mat <- matrix(1, 3, 3, dimnames = list(grp_order, grp_order))
  pairs <- strsplit(as.character(dunn_res$Comparison), " - ")
  for (i in seq_along(pairs)) {
    p_mat[pairs[[i]][1], pairs[[i]][2]] <- dunn_res$P.adj[i]
    p_mat[pairs[[i]][2], pairs[[i]][1]] <- dunn_res$P.adj[i]
  }
  fig7_letters[[d]] <- multcompLetters(p_mat, compare = "<", threshold = 0.05)$Letters

  lab <- data.frame(genotype = names(fig7_letters[[d]]), letter = fig7_letters[[d]],
                    y = max(bc4_geno[[d]]) * ifelse(d %in% c("vt", "vs"), 1.6, 1.08))
  p <- ggplot(bc4_geno, aes(genotype, .data[[d]], colour = genotype)) +
    geom_boxplot(outlier.alpha = 0, width = 0.6) +
    geom_point(position = position_jitter(width = 0.15, seed = 1), alpha = 0.6, size = 2) +
    geom_text(data = lab, aes(genotype, y, label = letter), colour = "black", size = 5) +
    scale_colour_manual(values = c("Col-0" = "magenta", moca1 = "#21908c", "cerk1-2" = "#e6b800")) +
    labs(x = NULL, y = desc_label[d]) + theme_classic() + theme(legend.position = "none")
  if (d %in% c("vt", "vs")) p <- p + scale_y_log10()
  print(p)
  fig_name <- c(sss = "Fig7b", sd = "Fig7c", sot = "FigS9d", dfrt = "FigS9e", vt = "FigS9f", vs = "FigS9g")[d]
  ggsave(file.path(out_dir, paste0(fig_name, "_BC4_", d, ".png")), p, width = 3.5, height = 4.5, dpi = 300, bg = "white")
}
fig7_kw
fig7_letters
write.csv(fig7_kw, file.path(out_dir, "Fig7_S9_BC4_genotype_KruskalWallis.csv"), row.names = FALSE)



################################################################################
### BLOCK 13. Fig. 7d: PCA of the BC4 roots of the three genotypes
################################################################################
## A separate PCA fitted on the 45 BC4 roots only (not the atlas space).

published$Fig7d <- c(PC1 = 31.9, PC2 = 24.7)

pca_7d <- PCA(bc4_geno[, desc], scale.unit = TRUE, ncp = 6, graph = FALSE)
pca_7d$eig[1:2, ]
p7d <- fviz_pca_ind(pca_7d, habillage = bc4_geno$genotype, label = "none", addEllipses = TRUE,
                    ellipse.type = "confidence", title = "") +
  scale_colour_manual(values = c("Col-0" = "magenta", moca1 = "#21908c", "cerk1-2" = "#e6b800")) +
  scale_fill_manual(values = c("Col-0" = "magenta", moca1 = "#21908c", "cerk1-2" = "#e6b800")) +
  theme_classic()
print(p7d)
ggsave(file.path(out_dir, "Fig7d_PCA_BC4_genotypes.png"), p7d, width = 6, height = 5, dpi = 300, bg = "white")



################################################################################
### BLOCK 14. Fig. S11 and Table S5: mutant BC4 roots in the Ca2+tlas space
### Needs Blocks 4 and 5
################################################################################
## The 45 BC4 roots are projected into the atlas PCA of Block 4 and evaluated
## with the class models of Block 5. For each root and class: squared
## Mahalanobis distance (prior-free) and posterior membership (both priors).
## The Col-0 roots are the same 20 roots as in Block 5 and must give the same
## numbers (Table S3). Question: do the mutants move away from nutrients/salts
## specifically, or from all classes?

published$TableS5_Col0_d2 <- c(DAMP = 0.96, "Nutrient/salt" = 3.68)

Z_mut <- predict(pca, newdata = bc4_geno[, desc])$coord[, 1:K]

d2_mut <- sapply(classes, function(cl) mahalanobis(Z_mut, centroids[[cl]], covs[[cl]]))
ll_mut <- sapply(classes, function(cl) dmvnorm(Z_mut, centroids[[cl]], covs[[cl]], log = TRUE))
pf_mut <- exp(sweep(ll_mut, 2, log(priors_freq), "+")); pf_mut <- pf_mut / rowSums(pf_mut)
pu_mut <- exp(sweep(ll_mut, 2, log(priors_unif), "+")); pu_mut <- pu_mut / rowSums(pu_mut)
knn_mut <- t(apply(Z_mut, 1, function(z) {
  dd <- colSums((t(as.matrix(X_act)) - z)^2)
  prop.table(table(factor(d1$treatment[order(dd)[1:k_nn]], levels = levels(d1$treatment))))
}))

mut_per_root <- data.frame(genotype = bc4_geno$genotype, line = bc4_geno$line,
                           setNames(as.data.frame(d2_mut), paste0("d2_", classes)),
                           setNames(as.data.frame(pf_mut), paste0("post_freq_", classes)),
                           setNames(as.data.frame(pu_mut), paste0("post_unif_", classes)),
                           check.names = FALSE)

tableS5_d2   <- aggregate(as.data.frame(d2_mut), list(genotype = bc4_geno$genotype), function(v) round(median(v), 2))
tableS5_pf   <- aggregate(as.data.frame(pf_mut), list(genotype = bc4_geno$genotype), function(v) round(mean(v), 2))
tableS5_pu   <- aggregate(as.data.frame(pu_mut), list(genotype = bc4_geno$genotype), function(v) round(mean(v), 2))
knn_mut_geno <- aggregate(as.data.frame(knn_mut), list(genotype = bc4_geno$genotype), function(v) round(mean(v), 2))
tableS5_d2
tableS5_pf
tableS5_pu
knn_mut_geno

## Kruskal-Wallis + Dunn on the nutrient/salt distance and posterior (Fig. S11b,c)
nut <- grep("^Nutrient", classes, value = TRUE)
s11_kw <- data.frame()
s11_letters <- list()
for (v in c(paste0("d2_", nut), paste0("post_freq_", nut))) {
  tmp <- data.frame(y = mut_per_root[[v]], genotype = mut_per_root$genotype)
  kw <- kruskal.test(y ~ genotype, data = tmp)
  s11_kw <- rbind(s11_kw, data.frame(variable = v, H = round(unname(kw$statistic), 2), p = signif(kw$p.value, 3)))
  dunn_res <- dunnTest(y ~ genotype, data = tmp, method = "bh")$res
  grp_order <- names(sort(tapply(tmp$y, tmp$genotype, median), decreasing = TRUE))
  p_mat <- matrix(1, 3, 3, dimnames = list(grp_order, grp_order))
  pairs <- strsplit(as.character(dunn_res$Comparison), " - ")
  for (i in seq_along(pairs)) {
    p_mat[pairs[[i]][1], pairs[[i]][2]] <- dunn_res$P.adj[i]
    p_mat[pairs[[i]][2], pairs[[i]][1]] <- dunn_res$P.adj[i]
  }
  s11_letters[[v]] <- multcompLetters(p_mat, compare = "<", threshold = 0.05)$Letters

  p <- ggplot(tmp, aes(genotype, y, colour = genotype)) +
    geom_boxplot(outlier.alpha = 0, width = 0.6) +
    geom_point(position = position_jitter(width = 0.15, seed = 1), alpha = 0.6, size = 2) +
    annotate("text", x = names(s11_letters[[v]]), y = max(tmp$y) * 1.08, label = s11_letters[[v]], size = 5) +
    scale_colour_manual(values = c("Col-0" = "magenta", moca1 = "#21908c", "cerk1-2" = "#e6b800")) +
    labs(x = NULL, y = v) + theme_classic() + theme(legend.position = "none")
  print(p)
  ggsave(file.path(out_dir, paste0(ifelse(grepl("^d2", v), "FigS11b_", "FigS11c_"), "BC4_genotypes.png")),
         p, width = 3.5, height = 4.5, dpi = 300, bg = "white")
}
s11_kw
s11_letters

## Fig. S11a: reference class ellipses with the three genotypes projected
mut_xy <- data.frame(Z_mut[, 1:2], genotype = bc4_geno$genotype)
geno_centroids <- aggregate(mut_xy[, 1:2], list(genotype = mut_xy$genotype), mean)
pS11a <- fviz_pca_ind(pca, habillage = group, invisible = "ind.sup", label = "none", alpha.ind = 0,
                      addEllipses = TRUE, ellipse.type = "t", title = "") +
  scale_colour_viridis_d(option = "E") + scale_fill_viridis_d(option = "E") +
  ## genotypes drawn with fixed colours (outside the class colour scale)
  geom_point(data = mut_xy[mut_xy$genotype == "Col-0", ],   aes(Dim.1, Dim.2), inherit.aes = FALSE, colour = "magenta", shape = 8,  size = 2.5) +
  geom_point(data = mut_xy[mut_xy$genotype == "moca1", ],   aes(Dim.1, Dim.2), inherit.aes = FALSE, colour = "#21908c", shape = 17, size = 2.5) +
  geom_point(data = mut_xy[mut_xy$genotype == "cerk1-2", ], aes(Dim.1, Dim.2), inherit.aes = FALSE, colour = "#e6b800", shape = 15, size = 2.5) +
  geom_point(data = geno_centroids, aes(Dim.1, Dim.2), inherit.aes = FALSE, shape = 23, size = 4,
             colour = "black", fill = c("magenta", "#21908c", "#e6b800")[match(geno_centroids$genotype, c("Col-0", "moca1", "cerk1-2"))]) +
  theme_classic() + scale_x_reverse() + scale_y_reverse()
print(pS11a)
ggsave(file.path(out_dir, "FigS11a_genotypes_in_atlas_space.png"), pS11a, width = 7, height = 6, dpi = 300, bg = "white")

write.csv(mut_per_root, file.path(out_dir, "TableS5_per_root.csv"), row.names = FALSE)
write.csv(tableS5_d2, file.path(out_dir, "TableS5b_median_d2.csv"), row.names = FALSE)
write.csv(tableS5_pf, file.path(out_dir, "TableS5a_mean_posterior_freq.csv"), row.names = FALSE)
write.csv(tableS5_pu, file.path(out_dir, "TableS5a_mean_posterior_unif.csv"), row.names = FALSE)
write.csv(knn_mut_geno, file.path(out_dir, "TableS5_knn_by_elicitor.csv"), row.names = FALSE)
write.csv(s11_kw, file.path(out_dir, "FigS11_KruskalWallis.csv"), row.names = FALSE)



################################################################################
### BLOCK 15. Fig. 7e and S10b: primary root growth
################################################################################
## Linear model: length ~ treatment x genotype + batch (batch = blocking factor).
## Estimated marginal means (emmeans): all pairwise comparisons among Control,
## Mock and BC4 within each genotype (3 tests per genotype), BH-adjusted. Type II ANOVA for the interaction. Repeated on
## log(length) (effects become ratios). Length is in cm; x 10 for mm.
## Fig. 7e: parental lines without the biosensor. Fig. S10b: reporter lines.

published$Fig7e <- c(BC4_vs_Control_mm = 2.8, t = 2.36, BC4_vs_Mock_mm = 3.6, t_Mock = 2.57, p_BH = 0.028,
                     F_interaction = 1.12, p_interaction = 0.35,
                     log_ratio = 1.14, log_ratio_Mock = 1.18, log_p_BH = 0.025,
                     log_F_interaction = 1.29, log_p_interaction = 0.27)

growth_parental <- droplevels(d10[d10$line == "parental", ])
table(growth_parental$genotype, growth_parental$treatment, growth_parental$batch)

fit_7e <- lm(length ~ treatment * genotype + batch, data = growth_parental)
anova_7e <- Anova(fit_7e, type = 2)
anova_7e
emm_7e <- emmeans(fit_7e, ~ treatment | genotype)
contr_7e <- pairs(emm_7e, reverse = TRUE, adjust = "BH")     # Mock-Control, BC4-Control, BC4-Mock
contr_7e
contr_7e_df <- as.data.frame(contr_7e)
contr_7e_df$estimate_mm <- round(contr_7e_df$estimate * 10, 2)

fit_7e_log <- lm(log(length) ~ treatment * genotype + batch, data = growth_parental)
anova_7e_log <- Anova(fit_7e_log, type = 2)
anova_7e_log
contr_7e_log <- pairs(emmeans(fit_7e_log, ~ treatment | genotype), reverse = TRUE,
                      adjust = "BH", type = "response")         # ratios
contr_7e_log
contr_7e_log_df <- as.data.frame(contr_7e_log)

p7e <- ggplot(growth_parental, aes(treatment, length * 10, colour = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, seed = 1), alpha = 0.5) +
  facet_wrap(~ genotype) + labs(x = NULL, y = "Primary root length (mm)") +
  scale_colour_manual(values = c(Control = "grey40", Mock = "black", BC4 = "magenta")) +
  theme_classic() + theme(legend.position = "none")
print(p7e)
ggsave(file.path(out_dir, "Fig7e_root_growth_parental.png"), p7e, width = 7, height = 4, dpi = 300, bg = "white")

## Fig. S10b: reporter lines
growth_reporter <- droplevels(d10[d10$line != "parental", ])
table(growth_reporter$genotype, growth_reporter$treatment, growth_reporter$batch)
fit_S10 <- lm(length ~ treatment * genotype + batch, data = growth_reporter)
Anova(fit_S10, type = 2)
contr_S10 <- pairs(emmeans(fit_S10, ~ treatment | genotype), reverse = TRUE, adjust = "BH")
contr_S10

pS10b <- ggplot(growth_reporter, aes(treatment, length * 10, colour = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, seed = 1), alpha = 0.5) +
  facet_wrap(~ paste(genotype, line)) + labs(x = NULL, y = "Primary root length (mm)") +
  scale_colour_manual(values = c(Control = "grey40", Mock = "black", BC4 = "magenta")) +
  theme_classic() + theme(legend.position = "none")
print(pS10b)
ggsave(file.path(out_dir, "FigS10b_root_growth_reporter.png"), pS10b, width = 6, height = 4, dpi = 300, bg = "white")

write.csv(as.data.frame(anova_7e), file.path(out_dir, "Fig7e_TypeII_ANOVA.csv"))
write.csv(contr_7e_df, file.path(out_dir, "Fig7e_pairwise_contrasts.csv"), row.names = FALSE)
write.csv(as.data.frame(contr_7e_log), file.path(out_dir, "Fig7e_pairwise_log_ratios.csv"), row.names = FALSE)
write.csv(as.data.frame(contr_S10), file.path(out_dir, "FigS10b_pairwise_contrasts.csv"), row.names = FALSE)



################################################################################
### BLOCK 16. Intensity traces: Fig. 1, S1 (Datasets 5, 6) and S7b,d, S9c (7-9)
################################################################################
## Datasets 5 and 6 are raw intensities. Here dF/F0 is computed with
## F0 = mean of all frames before stimulus (time < 0). CHECK that this is the
## F0 used for Fig. 1 and S1; change `f0_window` if not.
## Datasets 7-9 already contain mean and SD of dF/F0 across roots; they are
## plotted as deposited, with n_roots per phase shown in the facet title.

f0_window <- c(-Inf, 0)     # F0 = frames from the start of the recording up to t = 0

tr5 <- pivot_longer(d5, -time, names_to = "trace", values_to = "F")
tr5 <- tr5 %>% group_by(trace) %>%
  mutate(dFF0 = (F - mean(F[time >= f0_window[1] & time < f0_window[2]])) /
                 mean(F[time >= f0_window[1] & time < f0_window[2]])) %>% ungroup()
tr5$stimulus <- sub("_.*", "", tr5$trace)
tr5$roi      <- sub("^[^_]*_", "", tr5$trace)
p1 <- ggplot(tr5, aes(time, dFF0, colour = roi)) + geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~ stimulus) +
  labs(x = "Time (s)", y = "dF/F0", colour = "ROI") + theme_classic()
print(p1)
ggsave(file.path(out_dir, "Fig1_ROI_traces.png"), p1, width = 8, height = 3.5, dpi = 300, bg = "white")

tr6 <- pivot_longer(d6, -time, names_to = "trace", values_to = "F")
tr6 <- tr6 %>% group_by(trace) %>%
  mutate(dFF0 = (F - mean(F[time >= f0_window[1] & time < f0_window[2]])) /
                 mean(F[time >= f0_window[1] & time < f0_window[2]])) %>% ungroup()
tr6$stimulus <- tolower(sub("_.*", "", tr6$trace))
tr6$roi      <- sub("^[^_]*_", "", tr6$trace)
pS1 <- ggplot(tr6, aes(time, dFF0, colour = roi)) + geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~ stimulus) +
  labs(x = "Time (s)", y = "dF/F0", colour = "ROI") + theme_classic()
print(pS1)
ggsave(file.path(out_dir, "FigS1_centerline_vs_edge_traces.png"), pS1, width = 8, height = 3.5, dpi = 300, bg = "white")

## Datasets 7-9: phase = baseline / stimulus / eATP, taken from the `treatment` column
d7$phase <- ifelse(grepl("^baseline", d7$treatment), "baseline", ifelse(grepl("^atp", d7$treatment), "eATP", "stimulus"))
d8$phase <- ifelse(grepl("^baseline", d8$treatment), "baseline", ifelse(grepl("^atp", d8$treatment), "eATP", "stimulus"))
d9$phase <- ifelse(grepl("^baseline", d9$treatment), "baseline", ifelse(grepl("^atp", d9$treatment), "eATP", "stimulus"))

for (nm in c("d7", "d8", "d9")) {
  tr <- get(nm)
  p <- ggplot(tr, aes(time, mean)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
    geom_line(colour = "#b2182b") +
    geom_vline(xintercept = c(0, 900), linetype = "dashed") +
    facet_wrap(~ line_group) + labs(x = "Time (s)", y = "dF/F0 (mean +/- SD)") + theme_classic()
  print(p)
  fig_name <- c(d7 = "FigS7b_NaCl_traces", d8 = "FigS7d_FLG22_traces", d9 = "FigS9c_BC4_mock_traces")[nm]
  ggsave(file.path(out_dir, paste0(fig_name, ".png")), p, width = 11, height = 7, dpi = 300, bg = "white")
}

################################################################################
### BLOCK 17. EXPLORATORY, NOT REPORTED IN THE MANUSCRIPT
### Needs Blocks 4 and 5
################################################################################
## Two extra checks from the original analysis. They are kept for completeness
## but are not cited in the manuscript.
##  (a) Paired sign-flip permutation tests over the 20 BC4 roots: is the
##      nutrient/salt posterior higher, and the nutrient/salt Mahalanobis
##      distance lower, than for each other class? (B = 10,000, BH-adjusted)
##  (b) One-vs-rest logistic regression (nutrients/salts vs all other classes)
##      in the first K PCs, predicting P(nutrients/salts) for each BC4 root.

target_class <- nut
B <- 10000
set.seed(123)

perm_tests <- data.frame()
for (cl in setdiff(classes, target_class)) {
  for (what in c("posterior_freq", "mahalanobis")) {
    dvec <- if (what == "posterior_freq") post_freq[, target_class] - post_freq[, cl] else mahal_bc4[, cl] - mahal_bc4[, target_class]
    perm <- replicate(B, mean(sample(c(-1, 1), length(dvec), replace = TRUE) * dvec))
    perm_tests <- rbind(perm_tests, data.frame(metric = what, compare_to = cl, mean_diff = round(mean(dvec), 3),
                                               p_one_sided = mean(perm >= mean(dvec)),
                                               p_two_sided = mean(abs(perm) >= abs(mean(dvec)))))
  }
}
perm_tests$p_adj_BH_one_sided <- ave(perm_tests$p_one_sided, perm_tests$metric, FUN = function(p) p.adjust(p, "BH"))
perm_tests

y_bin <- factor(ifelse(group == target_class, "nutrients", "other"), levels = c("other", "nutrients"))
glm_fit <- glm(y ~ ., data = cbind(X_act, y = y_bin), family = binomial())
bc4_one_vs_rest <- data.frame(sup_id = rownames(X_sup),
                              P_nutrients = round(predict(glm_fit, newdata = X_sup, type = "response"), 3))
summary(bc4_one_vs_rest$P_nutrients)

write.csv(perm_tests, file.path(out_dir, "Exploratory_permutation_tests.csv"), row.names = FALSE)
write.csv(bc4_one_vs_rest, file.path(out_dir, "Exploratory_one_vs_rest_logistic.csv"), row.names = FALSE)



################################################################################
### BLOCK 18. Reproduction check and session information
### Needs Blocks 2-15
################################################################################
## Every number quoted in the manuscript, next to what this script produced.
## `match` is TRUE when they agree at the precision printed in the manuscript.

rep_contrast <- contr_7e_df[contr_7e_df$genotype == "Col-0" & contr_7e_df$contrast == "BC4 - Control", ]
rep_mock     <- contr_7e_df[contr_7e_df$genotype == "Col-0" & contr_7e_df$contrast == "BC4 - Mock", ]
rep_log      <- contr_7e_log_df[contr_7e_log_df$genotype == "Col-0" & contr_7e_log_df$contrast == "BC4 / Control", ]
rep_log_mock <- contr_7e_log_df[contr_7e_log_df$genotype == "Col-0" & contr_7e_log_df$contrast == "BC4 / Mock", ]
repro_check <- data.frame(
  item = c(paste("Fig. S4 H,", names(published$S4_H)),
           "Fig. 5a PC1 %", "Fig. 5a PC2 %", "Cumulative variance, K = 4 PCs",
           "CV accuracy", "CV balanced accuracy",
           paste("Table S2 recall,", names(published$recall)),
           "Table S3 BC4 assigned to nutrients/salts (freq)", "Table S3 BC4 assigned to DAMP (freq)",
           "Table S3 BC4 assigned to DAMP (uniform)",
           "LOEO elicitors correct (freq)", "LOEO elicitors correct (uniform)",
           "LOEO root accuracy (freq)", "LOEO root accuracy (uniform)",
           "Fig. S5c PC1 %", "Fig. S5c PC2 %",
           paste("Fig. 7 H,", names(published$Fig7_H)),
           "Fig. 7d PC1 %", "Fig. 7d PC2 %",
           "Table S5 Col-0 median d2, DAMP", "Table S5 Col-0 median d2, nutrients/salts",
           "Fig. 7e Col-0 BC4 vs Control (mm)", "Fig. 7e t, BC4 vs Control", "Fig. 7e Col-0 BC4 vs Mock (mm)",
           "Fig. 7e t, BC4 vs Mock", "Fig. 7e BH p (both)", "Fig. 7e interaction F", "Fig. 7e interaction p",
           "Fig. 7e log ratio BC4/Control", "Fig. 7e log ratio BC4/Mock", "Fig. 7e log BH p (both)",
           "Fig. 7e log interaction F", "Fig. 7e log interaction p"),
  manuscript = c(published$S4_H,
                 published$PC_var,
                 published$CV,
                 published$recall,
                 14, 4, 16,
                 published$LOEO[c("correct_freq", "correct_unif", "root_acc_freq", "root_acc_unif")],
                 published$S5c,
                 published$Fig7_H,
                 published$Fig7d,
                 published$TableS5_Col0_d2,
                 published$Fig7e),
  this_script = c(s4_kw$H[match(names(published$S4_H), s4_kw$descriptor)],
                  round(pca_atlas$eig[1, 2], 1), round(pca_atlas$eig[2, 2], 1), round(pca_atlas$eig[4, 3], 1),
                  cv_summary$accuracy[1], cv_summary$balanced_accuracy[1],
                  tableS2$recall_freq[match(names(published$recall), tableS2$class)],
                  tableS3$assigned_freq[tableS3$class == nut], tableS3$assigned_freq[tableS3$class == "DAMP"],
                  tableS3$assigned_unif[tableS3$class == "DAMP"],
                  loeo_summary$correct[loeo_summary$prior == "frequency"], loeo_summary$correct[loeo_summary$prior == "uniform"],
                  loeo_summary$root_accuracy[loeo_summary$prior == "frequency"], loeo_summary$root_accuracy[loeo_summary$prior == "uniform"],
                  round(pca_s5$eig[1, 2], 1), round(pca_s5$eig[2, 2], 1),
                  fig7_kw$H[match(names(published$Fig7_H), fig7_kw$descriptor)],
                  round(pca_7d$eig[1, 2], 1), round(pca_7d$eig[2, 2], 1),
                  tableS5_d2$DAMP[tableS5_d2$genotype == "Col-0"], tableS5_d2[[nut]][tableS5_d2$genotype == "Col-0"],
                  round(rep_contrast$estimate * 10, 1), round(rep_contrast$t.ratio, 2), round(rep_mock$estimate * 10, 1),
                  round(rep_mock$t.ratio, 2), round(rep_contrast$p.value, 3),
                  round(anova_7e["treatment:genotype", "F value"], 2), round(anova_7e["treatment:genotype", "Pr(>F)"], 2),
                  round(rep_log$ratio, 2), round(rep_log_mock$ratio, 2), round(rep_log$p.value, 3),
                  round(anova_7e_log["treatment:genotype", "F value"], 2), round(anova_7e_log["treatment:genotype", "Pr(>F)"], 2)))
## tolerance = half of the last digit printed in the manuscript (e.g. 2.37 -> 0.005)
decimals <- sapply(repro_check$manuscript, function(v) { f <- format(v, trim = TRUE); if (grepl("\\.", f)) nchar(sub(".*\\.", "", f)) else 0 })
decimals <- pmax(decimals, ifelse(repro_check$manuscript %% 1 == 0, 0, 2))   # "0.30" is stored as 0.3; use at least 2 decimals
repro_check$match <- abs(repro_check$manuscript - repro_check$this_script) <= 0.5 * 10^(-decimals) + 1e-9
repro_check
write.csv(repro_check, file.path(out_dir, "Reproduction_check.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
sessionInfo()



################################################################################
################################################################################
##
##  OUTPUT MAP (all files in `out_dir`)
##  ---------------------------------------------------------------------------
##  Figure / table        File(s)                                    Block
##  Fig. S4               FigS4_<descriptor>.png, FigS4_*.csv          2
##  Fig. 5a               Fig5a_PCA_biplot.png, eigenvalues, loadings  3
##  Fig. 5b               Fig5b_correlation_mosaic.png, *.csv          3
##  Fig. 6b               Fig6b_PCA_BC4_supplementary.png              4
##  Fig. 6c               Fig6c_cumulative_posterior_freq.png          5
##  Fig. S6b, S6c         FigS6b_*.png, FigS6c_*.png                   5
##  Table S3              TableS3_BC4_class_assignment.csv             5
##  Fig. S6d, Table S2    FigS6d_confusion_matrix.png, TableS2_*.csv   6
##  Table S4              TableS4_LOEO*.csv                            7
##  Fig. S6e              FigS6e_<descriptor>.png                      8
##  Fig. S5c              FigS5c_GCaMP3_projection.png                 9
##  Fig. S5d,e            FigS5d_*.png/.csv, FigS5e_platform_PCA.png   10
##  Fig. S8               FigS8_*.png/.csv                             11
##  Fig. 7b,c, S9d-g      Fig7b_*, Fig7c_*, FigS9d-g_*.png, *.csv      12
##  Fig. 7d               Fig7d_PCA_BC4_genotypes.png                  13
##  Fig. S11, Table S5    FigS11a-c_*.png, TableS5*.csv                14
##  Fig. 7e, S10b         Fig7e_*.png/.csv, FigS10b_*.png/.csv         15
##  Fig. 1, S1, S7b,d, S9c  Fig1_*, FigS1_*, FigS7b_*, FigS7d_*, FigS9c_*  16
##  (not in manuscript)   Exploratory_*.csv                            17
##  Check                 Reproduction_check.csv, sessionInfo.txt      18
##
##  METHOD CHOICES TO KNOW WHEN READING THE RESULTS
##  ---------------------------------------------------------------------------
##  - The PCA uses untransformed, standardized descriptors. Velocities are
##    right-skewed; log10 is used for display only.
##  - Class priors: frequency priors reflect how many roots were recorded per
##    class, not how likely an unknown stimulus is to belong to it. BC4 is
##    therefore reported under both priors (Table S3). Under frequency priors
##    BC4 leans to nutrients/salts; under uniform priors to DAMP. The
##    prior-free Mahalanobis distance puts almost all BC4 roots nearest DAMP.
##  - Cross-validation (Block 6) keeps the PCA fixed and resamples roots of
##    known elicitors; leave-one-elicitor-out (Block 7) is the stricter test of
##    how an unseen stimulus is placed and gives lower accuracy.
##  - With class_scheme <- "relabel" (Block 1) the numbers of Blocks 4-8 and 14
##    change; the reproduction check then no longer applies to them.
##  - Fig. 6b and S11a draw 95 % data ellipses (ellipse.type = "t"), as in the
##    original figure script; Fig. 5a draws confidence ellipses of the mean.
##
##  RESULT OF THE REPRODUCTION CHECK (run on the Sep 29 dataset file)
##  ---------------------------------------------------------------------------
##  Reproduced exactly: Fig. S4 H values; PCA variance (46.0 / 17.2 / 87.4 %);
##  CV accuracy 0.743 and balanced accuracy 0.676; all six recalls (Table S2);
##  all of Table S3; LOEO 12/17 (0.58) and 11/17 (0.56) (Table S4); Fig. S5c
##  axes 49.3 / 18.2 %; Fig. 7 H values (14.2, 18.7, 11.2); Fig. 7d 31.9 /
##  24.7 %; all of Table S5; all Fig. 7e growth statistics (all pairwise
##  comparisons within genotype, BH-adjusted).
##
##  NOTE ON COMPACT LETTERS
##  ---------------------------------------------------------------------------
##  Groupings are what matter; the letter names produced here can differ from
##  those in the figures.
##
##  NOT INCLUDED
##  ---------------------------------------------------------------------------
##  - The raw per-root trace processing that produced Datasets 7-9 (the
##    dF/F0 computation and the exclusion of distorted traces).
##
################################################################################
################################################################################
