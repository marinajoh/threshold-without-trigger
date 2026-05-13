# =============================================================================
# dag.R
#
# Builds the directed acyclic graph (DAG) figure used in the Methods chapter
# to visualise the causal structure of the analysis. The dagitty package is
# used to specify the graph formally and to run causal-graph analyses
# (paths, d-separation, adjustment sets); the figure itself is constructed
# manually with ggplot and ggforce to give precise control over node
# layout, colour coding, and edge styling.
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggforce)
library(dagitty)

# ---- Palette and base theme -------------------------------------------------

# Palette and base theme are duplicated across descriptive.R, modeltest.R,
# and dag.R so that each script can be run independently of the others.
# Any change here should be mirrored in those files.

palette_full <- c(
  "#6C8FD4",
  "#A78FD6",
  "#7FC8B2",
  "#B8D986",
  "#E6A86C",
  "#D89595",
  "#C7B486")

base_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(color = "black", face = "bold", size = 12),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"))

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

# ---- Formal DAG specification -----------------------------------------------

# Specifies the causal graph using dagitty syntax. Nodes are abbreviated
# names that map to the substantive concepts as follows:
#   cyber = Article 5 cyber recognition (the exposure of interest)
#   beh   = observed cyber behaviour (the outcome)
#   cred  = perceived NATO credibility (latent mediator)
#   crim  = Crimea annexation (co-occurring with the cutoff)
#   wales = other Wales Summit commitments (co-occurring with the cutoff)
#   trend = long-term trends (absorbed by the time-trend spline)
#
# The hypothesised pathway is cyber -> cred -> beh: the recognition
# affects perceived credibility, which affects cyber behaviour. The other
# arrows capture co-occurring developments that the post-2014 exposure
# indicator absorbs jointly with the cyber recognition.

nato_dag <- dagitty('dag {
  cyber [exposure]
  beh   [outcome]
  cred  [latent]

  cyber -> cred -> beh
  crim  -> beh
  wales -> beh
  trend -> beh
}')

# ---- Causal-graph diagnostics -----------------------------------------------

# Three checks validate the causal structure encoded in the DAG.
#
# paths() lists all directed and undirected paths between cyber and beh,
# confirming that the hypothesised mediator pathway is the only directed
# route from exposure to outcome through cred.
#
# dseparated() tests whether cyber and beh are d-separated (statistically
# independent) given that cred is held fixed. A TRUE result confirms that
# the model treats cred as a complete mediator.
#
# adjustmentSets() identifies which variables must be controlled for to
# estimate the total causal effect of cyber on beh. Used to verify that
# no observable confounders are omitted from the analytical model.

print(paths(nato_dag, from = "cyber", to = "beh"))
print(dseparated(nato_dag, "cyber", "beh", "cred"))
print(adjustmentSets(nato_dag, "cyber", "beh", effect = "total"))

# ---- Visual styling for nodes and edges -------------------------------------

# Node colours by category. Different categories indicate the role of each
# node in the analysis: the focal exposure, the latent mediator, the
# outcome, the co-occurring developments bundled with the cutoff, and
# the long-term trend absorbed by the spline.

node_colors <- setNames(
  c(palette_full[2], "#C8B7E2", palette_full[6], palette_full[7], palette_full[3]),
  c("Focal", "Latent", "Outcome", "Co-occurring", "Absorbed by spline"))

# Box dimensions for node rectangles in plot coordinates.
box_w <- 2.0
box_h <- 1.15

# ---- Node positions and labels ----------------------------------------------

# Each node has hardcoded x/y coordinates that place it visually in the
# figure. Top row (y = 4) contains the cutoff-period developments bundled
# in the post-2014 indicator. Middle (y = 2.5) is the latent mediator.
# Bottom row (y = 1) is the outcome and the long-term trend.

nodes <- tibble(
  name = c("crim", "cyber", "wales", "cred", "beh", "trend"),
  title = c(
    "Crimea annexation",
    "Article 5 cyber recognition",
    "Other Wales commitments",
    "Perceived credibility",
    "Cyber behaviour",
    "Long-term trends"),
  subtitle = c(
    "Feb to Mar 2014",
    "Wales para 72, 5 Sep 2014",
    "RAP and NRF, 4 to 5 Sep 2014",
    "Latent (M)",
    "Count, composition (Y)",
    "Absorbed by spline"),
  x = c(1.2, 3.5, 5.8, 3.5, 3.5, 6.3),
  y = c(4, 4, 4, 2.5, 1, 1),
  category = c(
    "Co-occurring", "Focal", "Co-occurring",
    "Latent", "Outcome", "Absorbed by spline")) %>%
  mutate(
    category = factor(
      category,
      levels = c("Focal", "Latent", "Outcome", "Co-occurring", "Absorbed by spline")),
    xmin = x - box_w / 2,
    xmax = x + box_w / 2,
    ymin = y - box_h / 2,
    ymax = y + box_h / 2,
    title_y = y + 0.16,
    subtitle_y = y - 0.20)

# ---- Helper: build rounded box polygons -------------------------------------

# Converts each node's bounding-box coordinates into a four-corner polygon
# tibble. ggforce's geom_shape later renders these with rounded corners.

rounded_box <- function(xmin, xmax, ymin, ymax, id) {
  tibble(
    x = c(xmin, xmax, xmax, xmin),
    y = c(ymin, ymin, ymax, ymax),
    id = id)
}

box_polygons <- nodes %>%
  rowwise() %>%
  do(rounded_box(.$xmin, .$xmax, .$ymin, .$ymax, .$name)) %>%
  ungroup() %>%
  left_join(nodes %>% dplyr::select(name, category), by = c("id" = "name"))

# ---- Edges between nodes ----------------------------------------------------

# Each row is one arrow in the DAG. Starting coordinates (x, y) and ending
# coordinates (xend, yend) are tuned manually so arrows terminate at box
# boundaries rather than node centres. The kind column distinguishes
# hypothesised mechanisms (cyber to cred to beh) from co-occurring
# developments and from the absorbed long-term trend.

edges <- tibble(
  x =    c(3.5,   3.5,   1.2,   5.8,   5.3),
  y =    c(3.425, 1.925, 3.425, 3.425, 1),
  xend = c(3.5,   3.5,   2.6,   4.4,   4.5),
  yend = c(3.075, 1.575, 1.575, 1.575, 1),
  kind = c(
    "Hypothesised", "Hypothesised", "Co-occurring",
    "Co-occurring", "Absorbed")) %>%
  mutate(
    kind = factor(
      kind,
      levels = c("Hypothesised", "Co-occurring", "Absorbed")))

# ---- DAG figure -------------------------------------------------------------

# Builds the figure layer by layer: the dashed bundle box at the top, the
# arrows between nodes, the node rectangles with rounded corners, and the
# title and subtitle text on each node. The latent mediator is drawn with
# a dashed border to mark it as unobserved.

p_dag <- ggplot() +
  annotate(
    "rect",
    xmin = 0.1, xmax = 6.9, ymin = 3.35, ymax = 4.65,
    fill = NA, color = "grey55",
    linetype = "dashed", linewidth = 0.4) +
  annotate(
    "text",
    x = 3.5, y = 4.83,
    label = "Bundled in post-2014 indicator (X)",
    color = "grey40", size = 3.7) +
  geom_segment(
    data = edges,
    aes(
      x = x, y = y, xend = xend, yend = yend,
      linetype = kind, linewidth = kind),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "grey25",
    show.legend = FALSE) +
  scale_linetype_manual(
    values = c(
      "Hypothesised" = "solid",
      "Co-occurring" = "solid",
      "Absorbed" = "dashed")) +
  scale_linewidth_manual(
    values = c(
      "Hypothesised" = 0.8,
      "Co-occurring" = 0.4,
      "Absorbed" = 0.4)) +
  geom_shape(
    data = box_polygons %>% filter(category != "Latent"),
    aes(x = x, y = y, group = id, fill = category),
    color = "grey20", linewidth = 0.4, radius = unit(2.5, "mm")) +
  geom_shape(
    data = box_polygons %>% filter(category == "Latent"),
    aes(x = x, y = y, group = id, fill = category),
    color = "grey20", linewidth = 0.5, linetype = "dashed", radius = unit(2.5, "mm")) +
  geom_text(
    data = nodes,
    aes(x = x, y = title_y, label = title),
    size = 3.7, color = "grey10", fontface = "bold",
    vjust = 0.5) +
  geom_text(
    data = nodes,
    aes(x = x, y = subtitle_y, label = subtitle),
    size = 3.2, color = "grey25",
    vjust = 0.5) +
  scale_fill_manual(values = node_colors, name = "Node") +
  coord_cartesian(xlim = c(0, 7.5), ylim = c(0.2, 5.05)) +
  base_theme +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("Figures/dag.pdf", p_dag, width = 10, height = 5.8)