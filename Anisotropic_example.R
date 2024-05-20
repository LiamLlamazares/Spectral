source("S2C.R")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
logkappa <- 0

# Defining vectors for anisotropy
v1 <- c(0.001, 0) # Isotropic Matern
v2 <- c(1, 0) # Anisotropic
v3 <- c(-1, 0)

# Defining spectrum for different values of v

S_fun1 <- function(omega) {
  S_aniso_matrix(
    logkappa = logkappa,
    v = v1,
    omega = omega,
    sigma = 1
  )
}
S_fun2 <- function(omega) {
  S_aniso_matrix(
    logkappa = logkappa,
    v = v2,
    omega = omega,
    sigma = 1
  )
}
S_fun3 <- function(omega) {
  S_aniso_matrix(
    logkappa = logkappa,
    v = v3,
    omega = omega,
    sigma = 1
  )
}

# Defining discretisation points
n <- c(1024, 1024)
L <- c(10, 10)
h <- L / n
x_ <- make_x(n, L)
omega_ <- make_omega(n, L)
omega <- as.matrix(expand.grid(omega_))

# Defining spectrum at each frequency
S_1 <- S_fun1(omega = omega)
S_2 <- S_fun2(omega = omega)
S_3 <- S_fun3(omega = omega)

# Defining data frame for spectrum and covariance for each field
S_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(omega_))),
    data.frame(s_1 = S_1, s_2 = S_2, s_3 = S_3)
  ) %>%
  pivot_longer(
    cols = c(s_1, s_2, s_3),
    names_to = "Type",
    values_to = "Value"
  )

C_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      s_1 = as.vector(S2C(S_1, n, h)),
      s_2 = as.vector(S2C(S_2, n, h)),
      s_3 = as.vector(S2C(S_3, n, h))
    )
  ) %>%
  pivot_longer(
    cols = c("s_1", "s_2", "s_3"),
    names_to = "Type",
    values_to = "Value"
  )


# Plotting spectrum and covariance
ggplot(S_df) +
  geom_raster(aes(w1, w2, fill = Value)) +
  scale_fill_gradientn(colours = brewer.pal(7, "Reds")) +
  # scale_fill_viridis(option="inferno")+
  coord_equal() +
  facet_wrap(vars(Type), nrow = 1) +
  xlab("omega1") +
  ylab("omega2") +
  labs(fill = "Spectrum")
txt_size <- 40
p <- ggplot(C_df) +
  geom_raster(aes(x1, x2, fill = Value)) +
  geom_contour(aes(x1, x2, z = Value)) +
  scale_fill_gradientn(name = "K",
    colours = c("#FFFFFFFF", "#FF0000FF"),
    limits = c(min(C_df$Value), max(abs(C_df$Value)))
  ) +
  facet_wrap(vars(Type), ncol = 3, as.table = FALSE) +
  labs(x = expression(x[1]), y = expression(x[2]), fill = "Covariance") + # create facets
  coord_equal() +
  theme(
    text = element_text(size = txt_size), plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
    strip.text = element_blank(), legend.key.height = unit(1.5, "cm")
  )+
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_y_continuous(breaks = c(-4, 0, 4))
p
ggsave(
  "Images/Spectral_covariance.pdf",
  p,
  height = 5,
  dpi = 300,
  width = 15.25
)



# Sampling
# Calculate samples
sample_1 <- S2sample(
  S = S_1,
  dim = n,
  h = h,
  seed = NULL,
  conjugate = FALSE
)
sample_2 <- S2sample(
  S = S_2,
  dim = n,
  h = h,
  seed = NULL,
  conjugate = FALSE
)
sample_3 <- S2sample(
  S = S_3,
  dim = n,
  h = h,
  seed = NULL,
  conjugate = FALSE
)

# Calculate spatial grid and assign samples to their location
x_ <- lapply(make_x_sampling(n, h), function(x) x-5)
samples_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      u1 = as.vector(sample_1),
      u2 = as.vector(sample_2),
      u3 = as.vector(sample_3)
    )
  ) %>%
  pivot_longer(
    cols = c("u1", "u2", "u3"),
    names_to = "Type",
    values_to = "Value"
  )

# Plot samples over grid
# Plot samples over grid
p_sample <- ggplot(samples_df) +
  geom_raster(aes(x1, x2, fill = Re(Value))) +
  scale_fill_gradientn(
    name = expression(u),
    colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"),
    limits = c(-max(abs(Re(samples_df$Value))), max(abs(Re(samples_df$Value))))
  ) +
  scale_x_continuous(breaks = c(-4, 0,4)) +
  scale_y_continuous(breaks = c( -4, 0,4)) +
  facet_wrap(vars(Type), ncol = 3, as.table = FALSE) +
  labs(x = expression(x[1]), y = expression(x[2])) +
  coord_equal() +
  theme(
    text = element_text(size = txt_size), plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
    strip.text = element_blank(), legend.key.height = unit(1.5, "cm")
  )
p_sample
ggsave(
  "Images/Spectral_sample_u.pdf",
  p_sample,
  height = 5,
  width = 15.25,
  dpi = 300
)
