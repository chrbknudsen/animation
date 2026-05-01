library(tidyverse)
library(gganimate)
library(gifski)

nx <- 140
ny <- 45
n_frames <- 150

# Antal bånd mellem linjerne
n_bands <- 6

# Antal linjer = antal bånd + 1
n_stripes <- n_bands + 1

width  <- 10
height <- 6

dt <- 0.025

k_spring <- 55
k_anchor <- 0.6
damping  <- 0.94

triangle_height <- 2.4
triangle_width  <- 2.2

grid <- expand_grid(
  row = seq_len(ny),
  col = seq_len(nx)
) |>
  mutate(
    x0 = (col - 1) / (nx - 1) * width,
    y0 = (row - 1) / (ny - 1) * height,
    x = x0,
    y = y0,
    vx = 0,
    vy = 0
  )

dx0 <- width / (nx - 1)
dy0 <- height / (ny - 1)

stripe_rows <- round(seq(4, ny - 3, length.out = n_stripes))

triangle_at <- function(frame) {
  tip_x <- -1.0 + (frame - 1) / (n_frames - 1) * (triangle_width + 1.0)
  tip_y <- height / 2
  
  tibble(
    x = c(tip_x, tip_x - triangle_width, tip_x - triangle_width),
    y = c(tip_y, tip_y - triangle_height / 2, tip_y + triangle_height / 2)
  )
}

inside_triangle <- function(x, y, tri) {
  tip_x  <- tri$x[1]
  tip_y  <- tri$y[1]
  base_x <- tri$x[2]
  half_h <- triangle_height / 2
  
  inside_x <- x >= base_x & x <= tip_x
  
  rel <- (tip_x - x) / (tip_x - base_x)
  allowed_dy <- half_h * rel
  
  inside_x & abs(y - tip_y) <= allowed_dy
}

frames <- vector("list", n_frames)
triangles <- vector("list", n_frames)

state <- grid

for (f in seq_len(n_frames)) {
  
  tri <- triangle_at(f)
  
  state <- state |>
    arrange(row, col) |>
    group_by(row) |>
    mutate(
      x_l = lag(x),
      y_l = lag(y),
      x_r = lead(x),
      y_r = lead(y)
    ) |>
    ungroup() |>
    arrange(col, row) |>
    group_by(col) |>
    mutate(
      x_d = lag(x),
      y_d = lag(y),
      x_u = lead(x),
      y_u = lead(y)
    ) |>
    ungroup() |>
    mutate(
      fx_l = if_else(is.na(x_l), 0, k_spring * ((x_l + dx0) - x)),
      fy_l = if_else(is.na(y_l), 0, k_spring * (y_l - y)),
      
      fx_r = if_else(is.na(x_r), 0, k_spring * ((x_r - dx0) - x)),
      fy_r = if_else(is.na(y_r), 0, k_spring * (y_r - y)),
      
      fx_d = if_else(is.na(x_d), 0, k_spring * (x_d - x)),
      fy_d = if_else(is.na(y_d), 0, k_spring * ((y_d + dy0) - y)),
      
      fx_u = if_else(is.na(x_u), 0, k_spring * (x_u - x)),
      fy_u = if_else(is.na(y_u), 0, k_spring * ((y_u - dy0) - y)),
      
      fx_anchor = k_anchor * (x0 - x),
      fy_anchor = k_anchor * (y0 - y),
      
      fx = fx_l + fx_r + fx_d + fx_u + fx_anchor,
      fy = fy_l + fy_r + fy_d + fy_u + fy_anchor,
      
      vx = damping * (vx + fx * dt),
      vy = damping * (vy + fy * dt),
      
      x = x + vx * dt,
      y = y + vy * dt
    )
  
  hit <- inside_triangle(state$x, state$y, tri)
  
  if (any(hit)) {
    tip_x  <- tri$x[1]
    tip_y  <- tri$y[1]
    base_x <- tri$x[2]
    half_h <- triangle_height / 2
    
    xh <- state$x[hit]
    yh <- state$y[hit]
    
    rel <- (tip_x - xh) / (tip_x - base_x)
    rel <- pmax(0, pmin(1, rel))
    
    side <- sign(yh - tip_y)
    side[side == 0] <- 1
    
    boundary_y <- tip_y + side * half_h * rel
    penetration <- boundary_y - yh
    
    state$y[hit] <- state$y[hit] + 0.35 * penetration
    state$vy[hit] <- state$vy[hit] + 7.5 * penetration
    state$vx[hit] <- state$vx[hit] + 1.8 * (1 - rel)
    
    state$vx[hit] <- pmin(state$vx[hit], 3.5)
    state$vy[hit] <- pmax(pmin(state$vy[hit], 4), -4)
  }
  
  # Kileformet massefortrængning
  tip_x  <- tri$x[1]
  tip_y  <- tri$y[1]
  base_x <- tri$x[2]
  half_h <- triangle_height / 2
  
  wedge <- state$x > base_x & state$x < tip_x + 2.4
  
  if (any(wedge)) {
    xw <- state$x[wedge]
    yw <- state$y[wedge]
    
    rel <- (tip_x - xw) / (tip_x - base_x)
    
    rel_extended <- pmax(-0.9, pmin(1, rel))
    wedge_half_h <- half_h * pmax(rel_extended, 0.05)
    
    dy <- yw - tip_y
    side <- sign(dy)
    side[side == 0] <- 1
    
    inside_wedge_band <- abs(dy) < wedge_half_h + 0.7
    
    influence <- exp(-((abs(dy) - wedge_half_h)^2 / 0.35)) *
      exp(-pmax(xw - tip_x, 0)^2 / 1.2)
    
    idx <- which(wedge)[inside_wedge_band]
    
    state$vx[idx] <- state$vx[idx] + 0.85 * influence[inside_wedge_band]
    state$vy[idx] <- state$vy[idx] + side[inside_wedge_band] * 1.15 * influence[inside_wedge_band]
  }
  
  # Rendering:
  # fjern punkter inde i trekanten
  visible <- state |>
    filter(row %in% stripe_rows) |>
    filter(!inside_triangle(x, y, tri))
  
  # Tilføj kontaktpunkter på trekantens side
  contact_points <- map_dfr(stripe_rows, function(r) {
    y_base <- grid |>
      filter(row == r) |>
      slice(1) |>
      pull(y0)
    
    if (abs(y_base - tip_y) > half_h) {
      return(NULL)
    }
    
    side <- sign(y_base - tip_y)
    if (side == 0) side <- 1
    
    slope <- half_h / triangle_width
    
    x_contact <- tip_x - abs(y_base - tip_y) / slope
    y_contact <- y_base
    
    if (x_contact < base_x || x_contact > tip_x) {
      return(NULL)
    }
    
    tibble(
      row = r,
      col = 0,
      x0 = x_contact,
      y0 = y_contact,
      x = x_contact,
      y = y_contact,
      vx = 0,
      vy = 0
    )
  })
  
  visible <- bind_rows(visible, contact_points) |>
    arrange(row, x)
  
  frames[[f]] <- visible |>
    mutate(frame = f)
  
  triangles[[f]] <- tri |>
    mutate(frame = f)
}

sim_data <- bind_rows(frames)
tri_data <- bind_rows(triangles)

p <- ggplot() +
  geom_path(
    data = sim_data,
    aes(x, y, group = interaction(frame, row)),
    linewidth = 0.7
  ) +
  geom_polygon(
    data = tri_data,
    aes(x, y, group = frame),
    fill = "grey75",
    colour = "black",
    linewidth = 0.9
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
  renderer = gifski_renderer("triangle_material_sim_6_bands.gif")
)