## *****************************************************************************
##
## Script name: balloonplot_split_example.R
##
## Purpose of script: Base script used to generate split p-value balloon plots
##   for Felgines et al. 2024 in Fig. 2B.
##
## Author: Calvin Matteoli
##
## Last Modified: 29-Aug-2024
##
## *****************************************************************************
##
## Notes: This is the script used to make advanced balloon plots where the WT
##   have split balloons representing two different p-values for two
##   comparisons.
##   Using different data there are a few things that need to be changed, they
##   are noted as "Edit according to dataset". This script is more advanced than
##   the base example, please consider using the other one first. For further
##   explanation please contact the corresponding authors.
##
## *****************************************************************************

library(rstudioapi)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(scatterpie)
library(viridis)
library(ggnewscale)
library(cowplot)
library(scales)

#### SETUP ####

# Changing wd to current script's path
setwd(dirname(getActiveDocumentContext()$path))

# Data loading
# Edit according to dataset
balloon.data <- fread("balloonplot_example_data_WT_vs_two_mutants_Fig2B.csv")

#### DATA PREPARATION ####
# Setting the sample order for the plot
# Edit according to dataset
sample.order <-
  c(
    "NRPD1-WT-1A",
    "NRPD1-WT-1B",
    "NRPD1-WT-2A",
    "NRPD1-WT-2B",
    "NRPD1-WT-3A",
    "NRPD1-WT-3B",
    "NRPD1-MUT1-1A",
    "NRPD1-MUT1-1B",
    "NRPD1-MUT1-2A",
    "NRPD1-MUT1-2B",
    "NRPD1-MUT1-3A",
    "NRPD1-MUT1-3B",
    "NRPD1-MUT2-1A",
    "NRPD1-MUT2-1B",
    "NRPD1-MUT2-2A",
    "NRPD1-MUT2-2B",
    "NRPD1-MUT2-3A",
    "NRPD1-MUT2-3B"
  )
# Setting the gene order for the plot
# Edit according to dataset
protein.order <-
  c(
    "NRPD1",
    "NRPD2",
    "NRPD3a",
    "NRPD3b",
    "NRPD4",
    "NRPD5a",
    "NRPD5b",
    "NRPD6b",
    "NRPD7a",
    "NRPD7b",
    "NRPD8b",
    "NRPD9a",
    "NRPD9b",
    "NRPD10",
    "NRPD11",
    "NRPD12",
    "RDR2",
    "CLSY1",
    "CLSY2",
    "CLSY3",
    "CLSY4",
    "SHH1"
  )

# Preparing for Fold Change calculation
# Edit according to dataset

WT.names <- c("NRPD1-WT-1A",
              "NRPD1-WT-1B",
              "NRPD1-WT-2A",
              "NRPD1-WT-2B",
              "NRPD1-WT-3A",
              "NRPD1-WT-3B")

MUT1.names <- c("NRPD1-MUT1-1A",
              "NRPD1-MUT1-1B",
              "NRPD1-MUT1-2A",
              "NRPD1-MUT1-2B",
              "NRPD1-MUT1-3A",
              "NRPD1-MUT1-3B")

MUT2.names <- c("NRPD1-MUT2-1A",
              "NRPD1-MUT2-1B",
              "NRPD1-MUT2-2A",
              "NRPD1-MUT2-2B",
              "NRPD1-MUT2-3A",
              "NRPD1-MUT2-3B")

# Calculating sums for FC
balloon.data$sumWT <- balloon.data[, rowSums(.SD), .SDcols = WT.names]
balloon.data$sumMUT1 <- balloon.data[, rowSums(.SD), .SDcols = MUT1.names]
balloon.data$sumMUT2 <- balloon.data[, rowSums(.SD), .SDcols = MUT2.names]

# FC calculations
# alpha is used to avoid log(0)
alpha = 1E-10
balloon.data[, log2FC_M9 := log2(sumWT+alpha) - log2(sumMUT1+alpha)]
balloon.data[, log2FC_M10 := log2(sumWT+alpha) - log2(sumMUT2+alpha)]
balloon.data[, sumWT := NULL]
balloon.data[, sumMUT1 := NULL]
balloon.data[, sumMUT2 := NULL]

# Converting to long format for ggplot
balloon.data.long <-
  as.data.table(
    gather(
      balloon.data,
      genotype,
      peptide,
      "NRPD1-WT-1A":"NRPD1-MUT2-3B",
      factor_key = TRUE
    )
  )

# Setting the order for the samples and genes using factors
balloon.data.long[, genotype := factor(as.character(genotype),
                                       levels = sample.order)]
balloon.data.long[, `Gene name` := factor(`Gene name`,
                                          levels = protein.order)]

# Log calculation for padj and number of peptides
balloon.data.long[, neg_log10_adj_pvalue_M9 := -log10(adj.pvalue_M9)]
balloon.data.long[, neg_log10_adj_pvalue_M10 := -log10(adj.pvalue_M10)]
balloon.data.long[, log_peptide := log10(peptide)]

##### Labelling #####
# Setting Main and Sub categories for the facet_grid function later
# This depends on the way the genotypes are labeled and should be changed
# accordingly
# Edit according to dataset
balloon.data.long[, Cross := substring(genotype, 1,
                                       nchar(as.character(`genotype`)) -
                                         3)]
balloon.data.long[, c("Main", "Sub") := tstrsplit(Cross, "-")]

# Setting the Main and Sub values as factors to order them
# Edit according to the dataset
balloon.data.long[, Main := factor(Main,
                                   levels = c("NRPD1"))]
balloon.data.long[, Sub := factor(Sub,
                                  levels = c("WT", "MUT1", "MUT2"))]

##### Setting the p-value thresholds #####
# Defining padj thresholds (>0.05, <0.05, <0.01, <0.001)
thresholds <- c(0.05, 0.01, 0.001)

# Binning the padj for plot coloring and ordering
# Edit according to dataset
balloon.data.long[, bin_neg_log10_adj_pvalue_M9 :=
                    cut(neg_log10_adj_pvalue_M9, breaks =
                          c(-Inf,-log10(thresholds), Inf))]
balloon.data.long[, bin_neg_log10_adj_pvalue_M10 :=
                    cut(neg_log10_adj_pvalue_M10, breaks =
                          c(-Inf,-log10(thresholds), Inf))]
balloon.data.long[, bin_neg_log10_adj_pvalue_M9 := paste0("M09_", bin_neg_log10_adj_pvalue_M9)]
balloon.data.long[, bin_neg_log10_adj_pvalue_M10 := paste0("M10_", bin_neg_log10_adj_pvalue_M10)]
balloon.data.long[, IDs := .I]

# Setting |FC|<2 rows to the same as adj.p-value>0.05
balloon.data.long[abs(log2FC_M9) < 1, bin_neg_log10_adj_pvalue_M9 := "M09_(-Inf,1.3]"]
balloon.data.long[abs(log2FC_M10) < 1, bin_neg_log10_adj_pvalue_M10 := "M10_(-Inf,1.3]"]

# Adding new columns corresponding to the possible padj categories for later
value <- 0.5
balloon.data.long <- balloon.data.long %>%
  count(bin_neg_log10_adj_pvalue_M9, IDs) %>%
  spread(bin_neg_log10_adj_pvalue_M9,n,fill = 0) %>%
  mutate_at(vars(-IDs),funs(.*value)) %>%
  left_join(balloon.data.long)
balloon.data.long <- balloon.data.long %>%
  count(bin_neg_log10_adj_pvalue_M10, IDs) %>%
  spread(bin_neg_log10_adj_pvalue_M10,n,fill = 0) %>%
  mutate_at(vars(-IDs),funs(.*value)) %>%
  left_join(balloon.data.long)
balloon.data.long <- as.data.table(balloon.data.long)

# Setting character representations of the log10 thresholds
# They should correspond to the thresholds set earlier and be ordered
log.thresholds.cut <- c("(-Inf,1.3]",
                       "(1.3,2]",
                       "(2,3]",
                       "(3, Inf]")
# Making an ordered factor making sure the values for the MUT10 are first
# This makes them appear on the right of the pie charts as they fill from right
# to left.
# Edit according to dataset
ltc_M9 <- paste0("M09_", log.thresholds.cut)
ltc_M10 <- paste0("M10_", log.thresholds.cut)
ltc_order <- factor(c(ltc_M10, ltc_M9), levels = c(ltc_M10, ltc_M9))
ltc_order <- ltc_order[ltc_order %in% colnames(balloon.data.long)]

# Setting variables for the ggforce functions
balloon.data.long[, xpos := rowid(`Gene name`)]
balloon.data.long[, ypos := as.numeric(`Gene name`)]

# Melting and cleaning the table
# This is only for the WT to have two rows with the two different padj
# categories
balloon.data.hacked <- melt(balloon.data.long,
                            id.vars = c("IDs"),
                            measure.vars = as.character(ltc_order),
                            variable.name = "type", value.name = "value"
                            )[value == 0.5] %>%
  left_join(balloon.data.long)
cols <- c("type", "Gene name", "genotype", "log_peptide", "peptide", "Cross",
          "Main", "Sub", "xpos", "ypos", "bin_neg_log10_adj_pvalue_M9",
          "bin_neg_log10_adj_pvalue_M10")
balloon.data.hacked <- balloon.data.hacked[, ..cols]
balloon.data.hacked[, type := factor(type, levels = ltc_order)]
balloon.data.hacked[, xpos := as.numeric(xpos)]
# Adding some margins between the genotypes
# Edit accotding to dataset
balloon.data.hacked[Sub == "MUT1", xpos := xpos + .25]
balloon.data.hacked[Sub == "MUT2", xpos := xpos + .5]


# Plot color and label setup
genPal <- colorRampPalette(colors = c("grey90", "darkorchid"))
# As we have two set of factors for MUT9 and 10, we repeat the colors to make
# the illusion of of having only 5 colors where in reality we have 10
palette <- rep(genPal(4), 2)
names(palette) <- ltc_order


# The biggest log_peptide value is 2.371068 and I want it to be 12pt on the
# graph
# This gives a size of 10.12202pt to a value of 2, giving the size for 100
# peptides
# This number works for the resolution I am aiming for
# Edit according to dataset
magic.number <- 12/2.371068

# Setting stuff up for the fill palette
# Edit according to dataset
pvalue.labels <- c(">0.05", "<0.05", "<0.01", "<0.001")
legend.palette <- genPal(4)
names(legend.palette) <- pvalue.labels

# This is used to set the linewidth of the pies/circles and to scale the legend
# font size
# Edit according to dataset
thickness.size.factor <- .5
balloon.data.hacked$rescaled <- rescale(sqrt(balloon.data.hacked$log_peptide),
                                        to = c(0, 2.371068))

#PLOT####

p <- ggplot() +
  # Plotting only the WT values as scatterpies as they need to be pies to be
  # split in two
  geom_scatterpie(
    data = as.data.frame(balloon.data.hacked[Sub == "WT"]),
    mapping = aes(
      y = ypos,
      x = xpos,
      r = rescaled/magic.number,
      linewidth = I(0.5*thickness.size.factor)
    ), long_format = TRUE, cols = "variable"
  ) +
  scale_fill_manual(values = palette,
                    labels = rep(pvalue.labels, 2),
                    guide = "none") +
  # Plotting the MUT data, each using its own padj category
  # We use geom circle else the scatterpie leaves a bar in the shape even when
  # filled at 100%
  geom_circle(data = as.data.frame(balloon.data.hacked[Sub == "MUT1"]),
              aes(x0 = xpos, y0 = ypos, r = rescaled/magic.number, fill = bin_neg_log10_adj_pvalue_M9, linewidth = I(0.5*thickness.size.factor))) +
  geom_circle(data = as.data.frame(balloon.data.hacked[Sub == "MUT2"]),
              aes(x0 = xpos, y0 = ypos, r = rescaled/magic.number, fill = bin_neg_log10_adj_pvalue_M10, linewidth = I(0.5*thickness.size.factor))) +
  theme_light() +
  # Setting the breaks only where we have data
  scale_x_continuous(breaks=unique(balloon.data.hacked$xpos), labels=rep("", 18), minor_breaks = NULL) +
  scale_y_reverse(breaks=c(1:22), labels=levels(balloon.data.long$`Gene name`[1]), minor_breaks = NULL) +
  # Zooming on both axes
  # We need coord_fixed or else we don't have circles but ellipses
  coord_fixed(ylim = c(21.75, 1.25), xlim = c(1.25, 18.25)) +
  guides(
    size = guide_legend(label.position = "bottom", title.position =
                          "top")
  ) +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
  )  ;p

# Creating an independent plot to use as legend as the original legend has a
# double scale
# Edit according to dataset
l.size <- 2/magic.number
custom.legend.data <- data.table(
  "text" = c("100", "10", "1", "<0.001", "<0.01", "<0.05", ">0.05 or\n|FC|<2"),
  "size" = c(l.size, 1/magic.number, 0.1/magic.number, rep(l.size, 4)),
  "xpos" = c(1, 1.85, 2.40, 4.5, 5.5, 6.5, 7.5),
  "ypos" = rep(1, 7))
custom.legend.data$text <- factor(custom.legend.data$text,
                                  levels = custom.legend.data$text)
custom.legend.palette <- c(rep("white", 3), rev(legend.palette))
names(custom.legend.palette) <- custom.legend.data$text
# Creating the legend plot
legend.plot <- ggplot(custom.legend.data) +
  geom_circle(mapping = aes(x0 = xpos, y0 = ypos, r = size, fill = text),
              linewidth = I(0.5*thickness.size.factor)) +
  geom_text(mapping = aes(x = xpos, y = ypos - .75, label = text, fontface = 2),
            size = 4*thickness.size.factor) +
  geom_text(mapping = aes(x = 1.7, y = 1.75, label = "spectral count",
                          fontface = 2), size = 5.5*thickness.size.factor) +
  geom_text(mapping = aes(x = 6.5, y = 1.75, label = "adj. p-value",
                          fontface = 2), size = 5.5*thickness.size.factor) +
  theme_light() +
  scale_fill_manual(values = custom.legend.palette, guide = "none") +
  coord_fixed(ylim = c(0, 2), xlim = c(-3.75, 13.25)) +
  scale_x_continuous(name = NULL, minor_breaks = NULL, labels = NULL, breaks = c(-5:13)) +
  scale_y_continuous("", minor_breaks = NULL, labels = rep("NRPD/E7a", 5), breaks = c(-2, -1, 0, 1, 2)) +
  geom_rect(mapping=aes(xmin=0, xmax=9.5, ymin=0, ymax=2), fill = NA, color="black", size = .5*thickness.size.factor) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "white"))

# Setting size values
# Edit according to dataset
wdt <- 600
hgt <- 700
rel.factor <- 1 #To change if changing width or height

# Making a plot with background grids to double check that the legend plot has
# the same scale as the original plot (check only after save)
plot.final <- plot_grid(
  NULL,
  p + theme(plot.background = element_rect(color = "black")),
  NULL,
  legend.plot + theme(plot.background = element_rect(color = "black")),
  ncol = 1, rel_heights = c(-1.5*rel.factor, 10, -2*rel.factor, 1)
); plot.final


#SAVE####
# Save the background grid check plot
ggsave(
  file = "example_WT_vs_MUTs_grid.pdf",
  plot = plot.final,
  units = "px",
  width = wdt,
  height = hgt,
  dpi = 105,
  bg = "white"
)

# Create the final plot where the legend doesn't have a background grid
plot.final.gridless <- plot_grid(
  NULL,
  p + theme(panel.grid.major.x = element_blank()) ,
  NULL,
  legend.plot + theme(panel.grid.major = element_blank(),
                      panel.border = element_blank()),
  ncol = 1, rel_heights = c(-1.5*rel.factor, 10, -2*rel.factor, 1)
); plot.final.gridless

# Saving the plot as PDF, SVG and PNG
# Edit according to dataset
ggsave(
  file = "example_WT_vs_MUTs_gridless.pdf",
  plot = plot.final.gridless,
  units = "px",
  width = wdt,
  height = hgt,
  dpi = 105,
  bg = "white"
)
ggsave(
  file = "example_WT_vs_MUTs_gridless.svg",
  plot = plot.final.gridless,
  units = "px",
  width = wdt,
  height = hgt,
  dpi = 105,
  bg = "white"
)
ggsave(
  file = "example_WT_vs_MUTs_gridless.png",
  plot = plot.final.gridless,
  units = "px",
  width = wdt,
  height = hgt,
  dpi = 105,
  bg = "white"
)

