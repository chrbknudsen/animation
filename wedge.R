library(tidyverse)

width  <- 10
height <- 6

left_x <- 1.2
top_y  <- 5.0
bot_y  <- 1.0
mid_y  <- (top_y + bot_y) / 2

n_bands <- 6

band_cols <- c(
  "#F94144", # 1 rød
  "#F3722C", # 2 orange
  "#F9C74F", # 3 gul
  "#90BE6D", # 4 grøn
  "#43AAE8", # 5 blå
  "#8E6FD1"  # 6 lilla
)

point_len  <- 3.2
band_width <- 0.75

orange_top_fraction <- 1 / 3

boundary <- tibble(
  id = seq_len(n_bands),
  base_x = left_x + (id - 1) * band_width,
  tip_x  = left_x + (id - 1) * band_width + point_len
)

# Orange ligger mellem base_x[2] og base_x[3].
# Cuttet lægges inde i orange, så orange får 1/3 af topbredden tilbage.
cut_x <- boundary$base_x[2] + (1 - orange_top_fraction) * band_width

make_band <- function(i) {
  outer <- boundary |> filter(id == i)
  
  if (i == 1) {
    return(tibble(
      band = i,
      x = c(outer$base_x, outer$tip_x, outer$base_x),
      y = c(top_y, mid_y, bot_y)
    ))
  }
  
  inner <- boundary |> filter(id == i - 1)
  
  tibble(
    band = i,
    x = c(
      inner$base_x,
      outer$base_x,
      outer$tip_x,
      outer$base_x,
      inner$base_x,
      inner$tip_x
    ),
    y = c(
      top_y,
      top_y,
      mid_y,
      bot_y,
      bot_y,
      mid_y
    )
  )
}

clip_right <- function(poly, cut_x) {
  pts <- poly |> select(x, y)
  n <- nrow(pts)
  out <- list()
  
  add_point <- function(x, y) {
    out[[length(out) + 1]] <<- tibble(x = x, y = y)
  }
  
  for (i in seq_len(n)) {
    p1 <- pts[i, ]
    p2 <- pts[ifelse(i == n, 1, i + 1), ]
    
    in1 <- p1$x >= cut_x
    in2 <- p2$x >= cut_x
    
    if (in1 && in2) {
      add_point(p2$x, p2$y)
    }
    
    if (in1 && !in2) {
      t <- (cut_x - p1$x) / (p2$x - p1$x)
      add_point(cut_x, p1$y + t * (p2$y - p1$y))
    }
    
    if (!in1 && in2) {
      t <- (cut_x - p1$x) / (p2$x - p1$x)
      add_point(cut_x, p1$y + t * (p2$y - p1$y))
      add_point(p2$x, p2$y)
    }
  }
  
  if (length(out) == 0) return(NULL)
  
  bind_rows(out) |>
    mutate(band = poly$band[1])
}

bands <- map_dfr(seq_len(n_bands), make_band)

clipped_bands <- bands |>
  group_by(band) |>
  group_split() |>
  map_dfr(clip_right, cut_x = cut_x)

p <- ggplot() +
  geom_polygon(
    data = clipped_bands,
    aes(x, y, group = band, fill = factor(band)),
    colour = "black",
    linewidth = 1.1
  ) +
  scale_fill_manual(values = band_cols, guide = "none") +
  coord_fixed(
    xlim = c(cut_x - 0.25, max(boundary$tip_x) + 0.5),
    ylim = c(0, height),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  filename = "decorated_wedge_cut.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

p