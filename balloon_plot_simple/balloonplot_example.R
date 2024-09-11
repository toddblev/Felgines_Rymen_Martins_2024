## *****************************************************************************
##
## Script name: balloonplot_example.R
##
## Purpose of script: Base script used to generate balloon plots for 
##   Felgines et al. 2024 in Fig. 3B and 3D.
##
## Author: Calvin Matteoli
##
## Last Modified: 29-Aug-2024
##
## *****************************************************************************
##
## Notes: This is the base script used for most of the balloon plots in the 
##   paper. It uses the data from figure 3D as an example. Note that the 
##   "WT Col-0" is a negative control not shown in the figure of the paper, the 
##   statistics refer to the other two samples. Using different data there are a
##   few things that need to be changed, they are noted as 
##   "Edit according to dataset". For further explanation please contact the 
##   corresponding authors.
##
## *****************************************************************************

library(rstudioapi)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)

#### SETUP ####

# Changing wd to current script's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data loading
# Edit according to dataset
balloon.data <- fread("balloon_plot_example_data_3xHA-CLSY1_Fig3D.txt")

#### DATA PREPARATION ####
# Setting the sample order for the plot
# Edit according to dataset
sample.order <-
  c(
    "Col_0_A",
    "Col_0_B",
    "Col_0_C",
    "NRPD1_WT_A",
    "NRPD1_WT_B",
    "NRPD1_WT_C",
    "NRPD1_MUT_A",
    "NRPD1_MUT_B",
    "NRPD1_MUT_C"
  )
# Setting the gene order for the plot
# Edit according to dataset
gene.order <-
  c(
    "CLSY1",
    "CLSY2",
    "CLSY3",
    "CLSY4",
    "SHH1",
    "RDR2",
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
    "NRPD12"
  )

# Preparing for Fold Change calculation
# Edit according to dataset

WT.names <- c("NRPD1_WT_A",
              "NRPD1_WT_B",
              "NRPD1_WT_C")

MUT.names <- c("NRPD1_MUT_A",
               "NRPD1_MUT_B",
               "NRPD1_MUT_C")

# Calculating sums for FC
balloon.data$sumWT <- balloon.data[, rowSums(.SD), .SDcols = WT.names]
balloon.data$sumMUT <- balloon.data[, rowSums(.SD), .SDcols = MUT.names]

# FC calculations
# alpha is used to avoid log(0)
alpha = 1E-10
balloon.data[, log2FC := log2(sumWT+alpha) - log2(sumMUT+alpha)]
balloon.data[, sumWT := NULL]
balloon.data[, sumMUT := NULL]

# Converting to long format for ggplot
balloon.data.long <-
  as.data.table(gather(
    balloon.data,
    genotype,
    peptide,
    all_of(sample.order),
    factor_key = TRUE
  ))

# Setting the order for the samples and genes using factors
balloon.data.long[, genotype := factor(as.character(genotype),
                                       levels = sample.order)]
balloon.data.long[, `Gene name` := factor(`Gene name`,
                                          levels = gene.order)]

# Log calculation for padj and number of peptides
balloon.data.long[, neg_log10_adj_pvalue := -log10(adj_pvalue)]
balloon.data.long[, log_peptide := log10(peptide)]

##### Labelling #####
# Setting Main and Sub categories for the facet_grid function later
# This depends on the way the genotypes are labeled and should be changed 
# accordingly
# Edit according to dataset
balloon.data.long[, Cross :=
                    substring(genotype, 1, nchar(as.character(`genotype`)) - 2)]
balloon.data.long[Cross == "Col_0", Main := "none\nWT"]
balloon.data.long[Cross == "Col_0", Sub := "Col-0"]
balloon.data.long[Cross == "NRPD1_WT", Main := "3xHA-CLSY1\nNRPD1-3xF"]
balloon.data.long[Cross == "NRPD1_WT", Sub := "(CYC-YPMF)"]
balloon.data.long[Cross == "NRPD1_MUT", Main := "3xHA-CLSY1\nNRPD1-3xF"]
balloon.data.long[Cross == "NRPD1_MUT", Sub := "(AAA-YPMF)"]


# Setting the Main and Sub values as factors to order them
# We can use unique directly as all variables are ordered in the table like we
# want them for the plot
# If this is not your case please edit according to the dataset
balloon.data.long[, Main := factor(Main,
                                   levels = unique(balloon.data.long$Main))]
balloon.data.long[, Sub := factor(Sub,
                                  levels = unique(balloon.data.long$Sub))]

##### Setting the p-value thresholds #####
# Defining padj thresholds (>0.05, <0.05, <0.01, <0.001)
thresholds <- c(0.05, 0.01, 0.001)
# Set breaks for the plot
t.breaks <- c(-Inf, -log10(thresholds), Inf)

# Binning the padj for plot coloring and ordering
balloon.data.long[, bin_neg_log10_adj_pvalue :=
                    cut(neg_log10_adj_pvalue, breaks =
                          t.breaks)]
# Setting character representations of the log10 thresholds
# They should correspond to the thresholds set earlier and be ordered
log.thresholds.cut <- c("(-Inf,1.3]",
                        "(1.3,2]",
                        "(2,3]",
                        "(3, Inf]")
# Setting thresholds order according to the previous variable
balloon.data.long[, bin_neg_log10_adj_pvalue := factor(bin_neg_log10_adj_pvalue,
                                                       levels = log.thresholds.cut)]

# Setting |FC|<2 rows to the same as adj.p-value>0.05
balloon.data.long[abs(log2FC) < 1, bin_neg_log10_adj_pvalue := "(-Inf,1.3]"]

# Plot color and label setup
pvalue.labels <- c("<0.001", "<0.01", "<0.05", ">0.05\n|log2FC| < 1")
genPal <- colorRampPalette(colors = c("grey90", "darkorchid"))
palette <- genPal(4)
names(palette) <- log.thresholds.cut

#### PLOT ####
# Processing plot
p <-
  ggplot(data = balloon.data.long) +
           # size.range = c(0, 12)) +
  geom_point(mapping = aes(x = genotype,
                           y = `Gene name`,
                           size = log_peptide,
                           fill = bin_neg_log10_adj_pvalue),
             shape = 21,
             show.legend=TRUE) +
  scale_size_continuous(
    breaks = c("100" = 2, "10" = 1, "1" = 0),
    range = c(0, 10),
    limits = c(0, 2.8)
  ) +
  scale_fill_manual(values = palette,
                    labels = pvalue.labels,
                    drop = FALSE,
                    na.translate = F,
                    guide = guide_legend(
                      override.aes = list(size = 10, fill = rev(palette)),
                      label.position = "bottom",
                      title.position = "top",
                      reverse = FALSE,
                      order = 2
                    )
                    ) +
  scale_y_discrete(limits = rev) +
  labs(size = "spectral count", fill = "adj. p-value") +
  facet_wrap(Main ~ Sub, nrow = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.title.x = ggplot2::element_blank(), 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    legend.position = "bottom",
    legend.justification = c("left", "bottom"),
    legend.box.margin = margin(l = -30),
    legend.spacing.x = unit(0, 'cm'),
    legend.margin = margin(t = 0, l = 10),
    legend.text = element_text(face = "bold"),
    legend.title.align = 0.5,
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold")
  ) +
  guides(
    size = guide_legend(
      label.position = "bottom",
      title.position =
        "top",
      order = 1
    )
  ); p

# Saving the plot as SVG, PNG and PDF
# Edit according to dataset
ggsave(
  file = "example_3xHA-CLSY1_Fig3D.svg",
  plot = p,
  units = "px",
  width = 390,
  height = 764.5,
  dpi = 105,
  bg = "white"
)
ggsave(
  file = "example_3xHA-CLSY1_Fig3D.png",
  plot = p,
  units = "px",
  width = 390,
  height = 764.5,
  dpi = 105,
  bg = "white"
)
ggsave(
  file = "example_3xHA-CLSY1_Fig3D.pdf",
  plot = p,
  units = "px",
  width = 390,
  height = 764.5,
  dpi = 105,
  bg = "white"
)
