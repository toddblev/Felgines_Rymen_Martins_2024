## *****************************************************************************
##
## Script name: volcano_example.R
##
## Purpose of script: Base script used to generate volcano plots for 
##   Felgines et al. 2024 in Fig. 2C, 2D, 3A and 3C.
##
## Author: Calvin Matteoli
##
## Last Modified: 29-Aug-2024
##
## *****************************************************************************
##
## Notes: This is the base script used for all the volcano plots in the paper. 
##   It uses the data from figure 2C as an example. Using different data there 
##   are a few things that need to be changed, they are noted as 
##   "Edit according to dataset". For further explanation please contact the 
##   corresponding authors.
##
## *****************************************************************************

library(ggplot2)
library(data.table)
library(ggrepel)

#### SETUP ####

# Changing wd to current script's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data loading
# Edit according to dataset
volcano.data <- fread("volcano_example_data_Fig2C.csv")
# Annotations loading
# They are used to set the which points will have names and what are the
# different categories. This information is used in the plots later for labels,
# coloring and size of the points
annotations <- fread("volcano_annotations.csv")
# Setting the legend order in the order of the annotations for the categories
# All non-annotated proteins will be set as "Other protein", we need to draw 
# them first to make the other categories pop but we want to show them last in 
# the legend
draw.order <- c("Other protein", unique(annotations$Category))
legend.order <- c(unique(annotations$Category), "Other protein")
annotations$Category <- factor(annotations$Category, levels = draw.order)

# Adding the annotations to the data
volcano.data.annot <- merge.data.table(volcano.data, annotations, 
                                       by = "accession", all.x = TRUE)
# Adding the missing Category data
volcano.data.annot[is.na(Category), Category := "Other protein"]
# Ordering to draw the categories in order
setorder(volcano.data.annot, "Category")

# Setting colors for the categories
# Edit according to dataset
palette <- c("#b99aff", "#649468", "#f2f32f", "#bebebe")
names(palette) <- legend.order

# Adding outline to specific points
# Edit according to dataset
volcano.data.annot[, outline := "No"]
volcano.data.annot[Category == "RDR2", outline := "Yes"]
# Setting outline palette
outline.palette <- c("black", alpha("red", 0))
names(outline.palette) <- c("Yes", "No")

# Setting size values
volcano.data.annot[, size := "Small"]
volcano.data.annot[Category != "Other protein", size := "Big"]
size.values <- c(3, 1.5)
names(size.values) <- c("Big", "Small")

# Volcano volcano.data
# Edit according to dataset
p <- ggplot(data = volcano.data.annot, 
            aes(
              x = Log2FC, 
              y = -log10(adj_pvalue), 
              fill = Category,
              label = `Plot Name`,
              color = outline,
              size = size)
            ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black") +
  geom_point(shape = 21) +
  scale_fill_manual(values = palette, breaks = legend.order) +
  scale_color_manual(values = outline.palette, guide="none") +
  scale_size_manual(values = size.values, guide="none") +
  geom_text_repel(box.padding = 4, max.overlaps = Inf, color = "black", 
                  nudge_y = 0.05, nudge_x = 0.05, force = 4,
                  angle        = 0,
                  hjust        = 0,
                  segment.size = 0.2,
                  min.segment.length = 0) +
  xlim(c(-8, 8)) +
  ylim(c(0, 18)) +
  theme_classic(); p

# Saving the plot as SVG, PNG and PDF
# Edit according to dataset
ggsave(
  file = "example_Fig2_AAA-YPMF_vs_WT.svg", plot = p, 
  width = 600, height = 450, units = "px", dpi = 100,
  bg = "white"
)
ggsave(
  file = "example_Fig2_AAA-YPMF_vs_WT.png", plot = p, 
  width = 600, height = 450, units = "px", dpi = 100,
  bg = "white"
)
ggsave(
  file = "example_Fig2_AAA-YPMF_vs_WT.pdf", plot = p, 
  width = 600, height = 450, units = "px", dpi = 100,
  bg = "white"
)

