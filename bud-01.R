library(tidyverse)
library(gganimate)
library(gifski)

# -----------------------------
# Parametre
# -----------------------------

n_stripes <- 6 #opr 16
n_points  <- 120
n_frames  <- 180

width  <- 10
height <- 6

dt <- 0.035

k_spring <- 65       # fjederstyrke langs striberne
k_anchor <- 18       # træk tilbage mod oprindelig y-position
damping  <- 0.88     # hastighedsdæmpning

triangle_height <- 1.8
triangle_width  <- 1.5

# -----------------------------
# Startpositioner
# -----------------------------

state <- expand_grid(
  stripe = seq_len(n_stripes),
  point  = seq_len(n_points)
) |>
  mutate(
    x  = seq(0, width, length.out = n_points)[point],
    y0 = stripe / (n_stripes + 1) * height,
    y  = y0,
    vx = 0,
    vy = 0
  )

rest_dx <- width / (n_points - 1)

# -----------------------------
# Hjælpefunktion: trekant
# -----------------------------

triangle_at <- function(frame) {
  tip_x <- -1.5 + frame / n_frames * 8
  tip_y <- height / 2
  
  tibble(
    x = c(tip_x, tip_x - triangle_width, tip_x - triangle_width),
    y = c(tip_y, tip_y - triangle_height / 2, tip_y + triangle_height / 2)
  )
}

# Punkt inde i trekant med spids mod højre
inside_triangle <- function(x, y, tri) {
  tip_x <- tri$x[1]
  tip_y <- tri$y[1]
  base_x <- tri$x[2]
  half_h <- triangle_height / 2
  
  inside_x <- x >= base_x & x <= tip_x
  
  # Trekanten bliver smallere mod spidsen
  relative_width <- (tip_x - x) / (tip_x - base_x)
  allowed_dy <- half_h * relative_width
  
  inside_x & abs(y - tip_y) <= allowed_dy
}

# Skub punkt ud af trekanten
push_out_triangle <- function(x, y, tri) {
  tip_x <- tri$x[1]
  tip_y <- tri$y[1]
  base_x <- tri$x[2]
  half_h <- triangle_height / 2
  
  relative_width <- (tip_x - x) / (tip_x - base_x)
  boundary_y <- tip_y + sign(y - tip_y) * half_h * relative_width
  
  # Hvis punktet ligger præcis på midten, vælg retning
  boundary_y[y == tip_y] <- tip_y + half_h * relative_width[y == tip_y]
  
  boundary_y
}

# -----------------------------
# Simulering
# -----------------------------

frames <- vector("list", n_frames)
triangles <- vector("list", n_frames)

for (f in seq_len(n_frames)) {
  
  tri <- triangle_at(f)
  
  state <- state |>
    group_by(stripe) |>
    arrange(point, .by_group = TRUE) |>
    mutate(
      x_left  = lag(x),
      y_left  = lag(y),
      x_right = lead(x),
      y_right = lead(y),
      
      # Fjederkraft fra venstre nabo
      fx_left = if_else(
        is.na(x_left),
        0,
        k_spring * ((x_left + rest_dx) - x)
      ),
      fy_left = if_else(
        is.na(y_left),
        0,
        k_spring * (y_left - y)
      ),
      
      # Fjederkraft fra højre nabo
      fx_right = if_else(
        is.na(x_right),
        0,
        k_spring * ((x_right - rest_dx) - x)
      ),
      fy_right = if_else(
        is.na(y_right),
        0,
        k_spring * (y_right - y)
      ),
      
      # Ankerkraft tilbage mod oprindelig højde
      fy_anchor = k_anchor * (y0 - y),
      
      fx = fx_left + fx_right,
      fy = fy_left + fy_right + fy_anchor,
      
      vx = damping * (vx + fx * dt),
      vy = damping * (vy + fy * dt),
      
      x = x + vx * dt,
      y = y + vy * dt
    ) |>
    ungroup()
  
  hit <- inside_triangle(state$x, state$y, tri)
  
  if (any(hit)) {
    state$y[hit] <- push_out_triangle(state$x[hit], state$y[hit], tri)
    state$vy[hit] <- 0
  }
  
  frames[[f]] <- state |>
    mutate(frame = f)
  
  triangles[[f]] <- tri |>
    mutate(frame = f)
}

sim_data <- bind_rows(frames)
tri_data <- bind_rows(triangles)

# -----------------------------
# Animation
# -----------------------------

p <- ggplot() +
  geom_path(
    data = sim_data,
    aes(x, y, group = interaction(frame, stripe)),
    linewidth = 0.55
  ) +
  geom_polygon(
    data = tri_data,
    aes(x, y, group = frame),
    fill = "grey75",
    colour = "black",
    linewidth = 0.8
  ) +
  coord_fixed(
    xlim = c(0, width),
    ylim = c(0, height),
    expand = FALSE
  ) +
  theme_void() +
  transition_manual(frame)

animate(
  p,
  nframes = n_frames,
  fps = 30,
  width = 900,
  height = 550,
  renderer = gifski_renderer("triangle_stripes_sim.gif")
)