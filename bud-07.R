library(tidyverse)
library(gganimate)
library(gifski)

# ============================================================
# Domæne
# ============================================================

nx <- 240
n_frames <- 170

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

# ============================================================
# Materialets vertikale udstrækning
# ============================================================

material_top    <- height - 0.55
material_bottom <- 0.55

material_height <- material_top - material_bottom

# ============================================================
# Kileparametre
# ============================================================

# Kilehøjde som andel af materialets højde
# 1.0 = hele højden
# 0.5 = halvdelen
wedge_height_fraction <- 1.0

# Parallel del som andel af materialets bredde
# 0.4 = 40% af bredden
parallel_fraction <- 0.25

# 1 = 90 graders spids
# >1 = spidsere
# <1 = stumpere
wedge_tip_factor <- 1.0

# ============================================================
# Udledte kiledimensioner
# ============================================================

wedge_height <- wedge_height_fraction * material_height

half_h <- wedge_height / 2

# Længden af selve spidsdelen
point_len <- wedge_tip_factor * half_h

# Parallel del
parallel_length <- parallel_fraction * width

# Samlet kilebredde
wedge_width <- parallel_length + point_len

# ============================================================
# Fysik
# ============================================================

dt <- 0.08
damping <- 0.86

k_smooth   <- 0.95
k_restore  <- 0.035
k_contact  <- 3.4
k_pressure <- 1.9

min_gap <- height / (n_bands + 2) * 0.32
constraint_iterations <- 22

# ============================================================
# Initialisering
# ============================================================

x_grid <- seq(0, width, length.out = nx)

y0_lines <- seq(
  material_bottom,
  material_top,
  length.out = n_lines
)

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

# ============================================================
# Kilegeometri
# ============================================================

wedge_at <- function(frame) {
  
  tip_x <- -0.8 +
    (frame - 1) / (n_frames - 1) *
    (wedge_width + 3.4)
  
  tip_y <- height / 2
  
  half_h <- wedge_height / 2
  
  point_len <- wedge_tip_factor * half_h
  
  shoulder_x <- tip_x - point_len
  left_x     <- tip_x - wedge_width
  
  tibble(
    x = c(
      left_x,
      shoulder_x,
      tip_x,
      shoulder_x,
      left_x
    ),
    y = c(
      tip_y + half_h,
      tip_y + half_h,
      tip_y,
      tip_y - half_h,
      tip_y - half_h
    )
  )
}

wedge_half_height_at_x <- function(x, wedge) {
  
  tip_x      <- wedge$x[3]
  shoulder_x <- wedge$x[2]
  left_x     <- wedge$x[1]
  
  half_h <- abs(wedge$y[1] - wedge$y[3])
  
  case_when(
    
    x >= left_x &
      x <= shoulder_x ~ half_h,
    
    x > shoulder_x &
      x <= tip_x ~
      half_h *
      (tip_x - x) /
      (tip_x - shoulder_x),
    
    TRUE ~ 0
  )
}

# ============================================================
# Hindrer lagene i at krydse hinanden
# ============================================================

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
        
        y_mat[j, bad] <-
          y_mat[j, bad] + correction
        
        y_mat[j - 1, bad] <-
          y_mat[j - 1, bad] - correction
        
        vy_mat[j, bad] <-
          vy_mat[j, bad] +
          0.25 * correction / dt
        
        vy_mat[j - 1, bad] <-
          vy_mat[j - 1, bad] -
          0.25 * correction / dt
      }
    }
  }
  
  state$y  <- as.vector(y_mat)
  state$vy <- as.vector(vy_mat)
  
  state
}

# ============================================================
# Simulation
# ============================================================

frames <- vector("list", n_frames)
wedges <- vector("list", n_frames)

for (f in seq_len(n_frames)) {
  
  wedge <- wedge_at(f)
  
  tip_x      <- wedge$x[3]
  tip_y      <- wedge$y[3]
  shoulder_x <- wedge$x[2]
  left_x     <- wedge$x[1]
  
  material_contact <- tip_x >= 0
  
  # ----------------------------------------------------------
  # Lokal krumning
  # ----------------------------------------------------------
  
  state <- state |>
    arrange(line, ix) |>
    group_by(line) |>
    mutate(
      y_l = lag(y),
      y_r = lead(y),
      curvature =
        coalesce(y_l, y) +
        coalesce(y_r, y) -
        2 * y
    ) |>
    ungroup()
  
  half_h_here <- wedge_half_height_at_x(state$x, wedge)
  
  inside_x <-
    state$x >= left_x &
    state$x <= tip_x
  
  dy <- state$y - tip_y
  
  inside_wedge <-
    inside_x &
    abs(dy) < half_h_here
  
  # ----------------------------------------------------------
  # Basisfysik
  # ----------------------------------------------------------
  
  force <-
    k_smooth * state$curvature +
    k_restore * (state$y0 - state$y)
  
  # ----------------------------------------------------------
  # Kontakt med kile
  # ----------------------------------------------------------
  
  if (any(inside_wedge)) {
    
    side <- sign(dy[inside_wedge])
    
    side[side == 0] <-
      ifelse(
        state$y0[inside_wedge] >= tip_y,
        1,
        -1
      )
    
    target_y <-
      tip_y +
      side * (
        half_h_here[inside_wedge] +
          min_gap * 0.5
      )
    
    penetration <-
      target_y -
      state$y[inside_wedge]
    
    force[inside_wedge] <-
      force[inside_wedge] +
      k_contact * penetration
  }
  
  # ----------------------------------------------------------
  # Trykbølge foran kilen
  # ----------------------------------------------------------
  
  if (material_contact) {
    
    ahead <-
      state$x >= tip_x &
      state$x <= tip_x + 5.2
    
    if (any(ahead)) {
      
      xdist <- state$x[ahead] - tip_x
      ydist <- state$y0[ahead] - tip_y
      
      side <- sign(ydist)
      side[side == 0] <- 1
      
      progress <- pmin(
        1,
        pmax(0, tip_x / 2.8)
      )
      
      pressure <-
        progress *
        exp(-(xdist^2) / 9.5) *
        exp(-(ydist^2) / 1.8)
      
      force[ahead] <-
        force[ahead] +
        side * k_pressure * pressure
    }
    
    # --------------------------------------------------------
    # Parallel del af kilen
    # --------------------------------------------------------
    
    around_parallel_part <-
      state$x >= pmax(left_x, 0) &
      state$x <= shoulder_x + 0.35 &
      abs(state$y0 - tip_y) <
      wedge_height * 1.25
    
    if (any(around_parallel_part)) {
      
      ydist <-
        state$y0[around_parallel_part] - tip_y
      
      side <- sign(ydist)
      side[side == 0] <- 1
      
      xdist_from_left <-
        state$x[around_parallel_part] -
        pmax(left_x, 0)
      
      body_pressure <-
        exp(-(xdist_from_left^2) / 6.0) *
        exp(-(ydist^2) / 2.2)
      
      force[around_parallel_part] <-
        force[around_parallel_part] +
        side *
        0.85 *
        k_pressure *
        body_pressure
    }
  }
  
  # ----------------------------------------------------------
  # Integrer bevægelse
  # ----------------------------------------------------------
  
  state <- state |>
    mutate(
      vy = damping * (vy + force * dt),
      y  = y + vy * dt,
      y  = pmin(
        pmax(y, 0.05),
        height - 0.05
      )
    ) |>
    select(line, ix, x, y0, y, vy)
  
  state <- enforce_spacing(state)
  
  frames[[f]] <- state |>
    mutate(frame = f)
  
  wedges[[f]] <- wedge |>
    mutate(frame = f)
}

# ============================================================
# Saml data
# ============================================================

line_data  <- bind_rows(frames)
wedge_data <- bind_rows(wedges)

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
      lower |>
        mutate(band = b),
      
      upper |>
        arrange(desc(x)) |>
        mutate(band = b)
    )
  })
})

# ============================================================
# Plot
# ============================================================

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
    data = wedge_data,
    aes(
      x = x,
      y = y,
      group = frame
    ),
    fill = "white",
    colour = "black",
    linewidth = 1.4
  ) +
  
  scale_fill_manual(
    values = band_cols,
    guide = "none"
  ) +
  
  coord_fixed(
    xlim = c(0, width),
    ylim = c(0, height),
    expand = FALSE
  ) +
  
  theme_void() +
  transition_manual(frame)

# ============================================================
# Render GIF
# ============================================================

animate(
  p,
  nframes = n_frames,
  fps = 30,
  width = 900,
  height = 550,
  renderer = gifski_renderer(
    "wedge_layered_material_sim.gif"
  )
)