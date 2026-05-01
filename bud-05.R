library(tidyverse)
library(gganimate)
library(gifski)

# ------------------------------------------------------------
# Parametre
# ------------------------------------------------------------

nx <- 220
n_frames <- 150

width  <- 10
height <- 6

n_bands <- 6
n_lines <- n_bands + 1

band_cols <- c(
  "#F94144",
  "#F3722C",
  "#F9C74F",
  "#90BE6D",
  "#43AA8B",
  "#577590"
)

triangle_width  <- 2.2
triangle_height <- 2.4

dt <- 0.08
damping <- 0.86

# Elastisk opførsel
k_smooth  <- 0.85   # udglatning langs x
k_restore <- 0.05   # svag tilbagevenden til oprindelig position
k_contact <- 2.6    # hvor hårdt trekanten skubber direkte
k_pressure <- 1.7   # trykfelt foran trekanten

# Hindrer at lagene krydser hinanden
min_gap <- height / (n_bands + 2) * 0.35
constraint_iterations <- 18

x_grid <- seq(0, width, length.out = nx)

# 7 linjer = grænser mellem 6 bånd
y0_lines <- seq(0.55, height - 0.55, length.out = n_lines)

state <- expand_grid(
  line = seq_len(n_lines),
  ix = seq_len(nx)
) |>
  mutate(
    x = x_grid[ix],
    y0 = y0_lines[line],
    y = y0,
    vy = 0
  )

triangle_at <- function(frame) {
  tip_x <- -1.0 + (frame - 1) / (n_frames - 1) * (triangle_width + 1.0)
  tip_y <- height / 2
  
  tibble(
    x = c(tip_x, tip_x - triangle_width, tip_x - triangle_width),
    y = c(tip_y, tip_y - triangle_height / 2, tip_y + triangle_height / 2)
  )
}

triangle_half_height_at_x <- function(x, tri) {
  tip_x  <- tri$x[1]
  base_x <- tri$x[2]
  half_h <- triangle_height / 2
  
  rel <- (tip_x - x) / (tip_x - base_x)
  
  ifelse(
    x >= base_x & x <= tip_x,
    half_h * pmax(0, pmin(1, rel)),
    0
  )
}

enforce_spacing <- function(state) {
  state <- state |> arrange(ix, line)
  
  y_mat  <- matrix(state$y,  nrow = n_lines, ncol = nx)
  vy_mat <- matrix(state$vy, nrow = n_lines, ncol = nx)
  
  for (iter in seq_len(constraint_iterations)) {
    for (j in 2:n_lines) {
      gap <- y_mat[j, ] - y_mat[j - 1, ]
      bad <- gap < min_gap
      
      if (any(bad)) {
        correction <- (min_gap - gap[bad]) / 2
        
        y_mat[j, bad]     <- y_mat[j, bad] + correction
        y_mat[j - 1, bad] <- y_mat[j - 1, bad] - correction
        
        vy_mat[j, bad]     <- vy_mat[j, bad] + 0.25 * correction / dt
        vy_mat[j - 1, bad] <- vy_mat[j - 1, bad] - 0.25 * correction / dt
      }
    }
  }
  
  state$y  <- as.vector(y_mat)
  state$vy <- as.vector(vy_mat)
  
  state
}

frames <- vector("list", n_frames)
triangles <- vector("list", n_frames)

for (f in seq_len(n_frames)) {
  
  tri <- triangle_at(f)
  
  tip_x  <- tri$x[1]
  tip_y  <- tri$y[1]
  base_x <- tri$x[2]
  
  material_contact <- tip_x >= 0
  
  state <- state |>
    arrange(line, ix) |>
    group_by(line) |>
    mutate(
      y_l = lag(y),
      y_r = lead(y),
      curvature = coalesce(y_l, y) + coalesce(y_r, y) - 2 * y
    ) |>
    ungroup()
  
  half_h_here <- triangle_half_height_at_x(state$x, tri)
  
  inside_x <- state$x >= base_x & state$x <= tip_x
  dy <- state$y - tip_y
  
  inside_tri <- inside_x & abs(dy) < half_h_here
  
  force <- k_smooth * state$curvature +
    k_restore * (state$y0 - state$y)
  
  if (any(inside_tri)) {
    side <- sign(dy[inside_tri])
    side[side == 0] <- ifelse(state$y0[inside_tri] >= tip_y, 1, -1)
    
    target_y <- tip_y + side * (half_h_here[inside_tri] + min_gap * 0.45)
    penetration <- target_y - state$y[inside_tri]
    
    force[inside_tri] <- force[inside_tri] + k_contact * penetration
  }
  
  if (material_contact) {
    ahead <- state$x >= tip_x & state$x <= tip_x + 5.0
    
    if (any(ahead)) {
      xdist <- state$x[ahead] - tip_x
      ydist <- state$y0[ahead] - tip_y
      
      side <- sign(ydist)
      side[side == 0] <- 1
      
      progress <- pmin(1, pmax(0, tip_x / width * 3.2))
      
      pressure <- progress *
        exp(-(xdist^2) / 8.5) *
        exp(-(ydist^2) / 1.6)
      
      force[ahead] <- force[ahead] + side * k_pressure * pressure
    }
  }
  
  state <- state |>
    mutate(
      vy = damping * (vy + force * dt),
      y = y + vy * dt,
      y = pmin(pmax(y, 0.05), height - 0.05)
    ) |>
    select(line, ix, x, y0, y, vy)
  
  state <- enforce_spacing(state)
  
  frames[[f]] <- state |>
    mutate(frame = f)
  
  triangles[[f]] <- tri |>
    mutate(frame = f)
}

line_data <- bind_rows(frames)
tri_data <- bind_rows(triangles)

band_data <- map_dfr(seq_len(n_frames), function(f) {
  dat <- line_data |>
    filter(frame == f)
  
  map_dfr(seq_len(n_bands), function(b) {
    lower <- dat |>
      filter(line == b) |>
      arrange(x)
    
    upper <- dat |>
      filter(line == b + 1) |>
      arrange(x)
    
    bind_rows(
      lower |> mutate(band = b),
      upper |> arrange(desc(x)) |> mutate(band = b)
    )
  })
})

p <- ggplot() +
  geom_polygon(
    data = band_data,
    aes(
      x = x,
      y = y,
      group = interaction(frame, band),
      fill = factor(band)
    ),
    colour = NA,
    alpha = 0.78
  ) +
  geom_path(
    data = line_data,
    aes(
      x = x,
      y = y,
      group = interaction(frame, line)
    ),
    linewidth = 0.75,
    colour = "black"
  ) +
  geom_polygon(
    data = tri_data,
    aes(x = x, y = y, group = frame),
    fill = "grey75",
    colour = "black",
    linewidth = 0.9
  ) +
  scale_fill_manual(values = band_cols, guide = "none") +
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
  renderer = gifski_renderer("triangle_layered_material_sim.gif")
)