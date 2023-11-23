source("S2C.R")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
# Example usages
v <- c(1,2)
logkappa <- 0

new <- function(omega) S_aniso(logkappa=logkappa, v=v, omega=omega)

S_aniso_matrix <- function(logkappa, v, omega_matrix) {
  # Apply S_aniso to each row of omega_matrix
  new <-function(omega) S_aniso(logkappa=logkappa, v=v, omega=omega)
  apply(omega_matrix, 1, new)
}

v1 <- c(0.001,0)
v2 <- c(1,0)
v3 <- c(-1,0)

S_fun1 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v1, omega=omega)
S_fun2 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v2, omega=omega)
S_fun3 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v3, omega=omega)

#Number of frequencies taken in each dimension
n <- c(256, 256)
L <- c(1, 1)
h <- L / n

#Subsivisions of x
x_ <- make_x(n, L)

#Build omega, a matrix of size n^2x2
omega_ <- make_omega(n, L)
omega <- as.matrix(expand.grid(omega_))

#Applies S to each omega to get vector of size n^2
S_1 <- S_fun1(omega = omega)
S_2 <- S_fun2(omega = omega)
S_3 <- S_fun3(omega = omega)

S_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(omega_))),
    data.frame(
      s1 = S_1,
      s2 = S_2,
      s3 = S_3
    )
  ) %>%
  pivot_longer(
    cols = c(s1, s2, s3),
    names_to = "Type", values_to = "Value"
  )

C_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      s1 = as.vector(S2C(S_truncated, n, h)),
      s2 = as.vector(S2C(S_folded, n, h)),
      s3 = as.vector(S2C(S_folded_cell, n, h))
    )
  ) %>%
  pivot_longer(
    cols = c("s1", "s2", "s3"),
    names_to = "Type", values_to = "Value"
  )
max(S_df["Value"])
ggplot(S_df) +
  geom_raster(aes(w1, w2, fill = Value)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"))+
  # scale_fill_viridis(option="inferno")+
  coord_equal()+
  facet_wrap(vars(Type),nrow=1) +
  xlab("omega1") + ylab("omega2") +
  labs(fill = "Spectrum")
ggplot(C_df) +
  geom_raster(aes(x1, x2, fill = Value)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"))+
  coord_equal()+
  facet_wrap(vars(Type),nrow =1) +
  xlab("u1") + ylab("u2") +
  labs(fill = "Covariance")




```{r}
S_df %>% group_by(Type) %>% summarise(N = sum(Value > 0))
```
```{r}
n <- c(256, 256)
L <- c(1, 1)
h <- L / n
x_ <- make_x_sampling(n, L)
omega_ <- make_omega_sampling(n, L)
omega <- as.matrix(expand.grid(omega_))
S_truncated <- S_fun1(omega = omega)
S_folded <- S_fun2(omega = omega)
S_folded_cell <- S_fun3(omega = omega)
```

```{r}
C_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      s1 = as.vector(S2sample(S_truncated, n, h)),
      s2 = as.vector(S2sample(S_folded, n, h)),
      s3 = as.vector(S2sample(S_folded_cell, n, h))
    )
  ) %>%
  pivot_longer(
    cols = c("s1", "s2", "s3"),
    names_to = "Type", values_to = "Value"
  )
```

```{r}
ggplot(C_df) +
  coord_equal()+
  geom_raster(aes(x1, x2, fill = Value)) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),limits=c(-20, 20))+
  facet_wrap(vars(Type),nrow =1) +
  xlab("u1") + ylab("u2") +
  labs(fill = "Sample")



























S_fun <- function(omega) {
  rho <- 10
  exp(-rowSums(omega^2) / rho^2 / 2) / rho^2 / (2 * pi)
}
n <- c(256, 256)
L <- c(1, 1)
h <- L / n
x_ <- make_x(n, L)
omega_ <- make_omega(n, L)
omega <- as.matrix(expand.grid(omega_))
S_truncated <- S_fun(omega = omega)
S_folded <- fold_spectrum(omega = omega, S_fun = S_fun, h)
S_folded_cell <- fold_spectrum(omega = omega, S_fun = S_fun, h, pointwise = FALSE)



S_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(omega_))),
    data.frame(
      Truncated = S_truncated,
      Folded = S_folded,
      "Folded_cell" = S_folded_cell
    )
  ) %>%
  pivot_longer(
    cols = c("Truncated", "Folded", "Folded_cell"),
    names_to = "Type", values_to = "Value"
  )

C_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      Truncated = as.vector(S2C(S_truncated, n, h)),
      Folded = as.vector(S2C(S_folded, n, h)),
      "Folded_cell" = as.vector(S2C(S_folded_cell, n, h))
    )
  ) %>%
  pivot_longer(
    cols = c("Truncated", "Folded", "Folded_cell"),
    names_to = "Type", values_to = "Value"
  )



ggplot(S_df) +
  geom_raster(aes(w1, w2, fill = Value)) +
  facet_grid(vars(Type)) +
  xlab("omega1") + ylab("omega2") +
  labs(fill = "Spectral density")



ggplot(C_df) +
  geom_raster(aes(x1, x2, fill = Value)) +
  facet_grid(vars(Type)) +
  xlab("x1") + ylab("x2") +
  labs(fill = "Covariance")



S_df %>% group_by(Type) %>% summarise(N = sum(Value > 0))


## 1D example


#S_fun <- function(omega) {
#  rho <- 2
#  exp(-rowSums(omega^2) * rho^2 / 2) / sqrt(2 * pi) * rho
#}
S_fun <- function(omega) {
  nu <- 1.5
  rho <- 3 / sqrt(8 * nu)
  v <- rho^(2 * nu) * gamma(nu) / gamma(nu + 0.5) / sqrt(4 * pi)
  rho^(2 * nu + 1) / (2 * pi) /
    (1 + rowSums(omega^2) * rho^2)^(nu + 0.5)
}
#S_fun <- function(omega) {
#  rho <- 0.1
#  1/(1-0.9*2*(rowSums(omega^2) * rho^2)^2+(rowSums(omega^2) * rho^2)^4)^2
#}
n <- c(512)
L <- c(20)
h <- L / n
x_ <- make_x(n, L)
omega_ <- make_omega(n, L)
omega <- as.matrix(expand.grid(omega_))
S_truncated <- fold_spectrum(omega = omega, S_fun = S_fun, h, fold = FALSE)
S_folded <- fold_spectrum(omega = omega, S_fun = S_fun, h)
S_folded_cell <- fold_spectrum(omega = omega, S_fun = S_fun, h, pointwise = FALSE)


S_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(omega_))),
    data.frame(
      Truncated = S_truncated,
      Folded = S_folded,
      "Folded_cell" = S_folded_cell
    )
  ) %>%
  pivot_longer(
    cols = c("Truncated", "Folded", "Folded_cell"),
    names_to = "Type", values_to = "Value"
  )

C_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      Truncated = as.vector(S2C(S_truncated, n, h)),
      Folded = as.vector(S2C(S_folded, n, h)),
      "Folded_cell" = as.vector(S2C(S_folded_cell, n, h))
    )
  ) %>%
  pivot_longer(
    cols = c("Truncated", "Folded", "Folded_cell"),
    names_to = "Type", values_to = "Value"
  )



ggplot(S_df, aes(w1, Value)) +
  geom_line(aes(col = Type)) +
  xlab("omega") + ylab("Spectral density")



ggplot(C_df, aes(x1, Value)) +
  geom_line(aes(col = Type)) +
  coord_cartesian(xlim = c(0, L / 2)) +
  xlab("x") + ylab("Covariance")



S_df %>% group_by(Type) %>% summarise(N = sum(Value > 0))


## Sampling example


n_sampling <- n * 2
x_ <- make_x_sampling(n_sampling, h)
omega_ <- make_omega_sampling(n_sampling, h)
omega <- as.matrix(expand.grid(omega_))
S_truncated <- S_fun(omega = omega)
S_folded <- fold_spectrum(omega = omega, S_fun = S_fun, h)
S_folded_cell <- fold_spectrum(omega = omega, S_fun = S_fun, h, pointwise = FALSE)
samples_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      Truncated = as.vector(S2sample(S_truncated, n_sampling, h, seed = 123L)),
      Folded = as.vector(S2sample(S_folded, n_sampling, h, seed = 123L)),
      "Folded_cell" = as.vector(S2sample(S_folded_cell, n_sampling, h, seed = 123L))
    )
  ) %>%
  pivot_longer(
    cols = c("Truncated", "Folded", "Folded_cell"),
    names_to = "Type", values_to = "Value"
  )

dim(samples_df)

ggplot(samples_df %>% filter(x1 < L), aes(x1)) +
  geom_line(aes(y = Re(Value), col = Type, lty="Real part")) +
  geom_line(aes(y = Im(Value), col = Type, lty="Imag part")) +
  coord_cartesian(xlim = c(0, L)) +
  xlab("x") + ylab("Sample")






if (is.matrix(samples_df)) {
  samples_df <- as.data.frame(samples_df)
  # Rename columns appropriately
  # names(samples_df) <- c("x1", "Type", "ReValue", "ImValue") # Adjust as per your matrix structure
}

# Extract real and imaginary parts if 'Value' is a complex number
# If 'Value' is not a complex number, replace 'Re(Value)' and 'Im(Value)' with the appropriate column names
samples_df$ReValue <- Re(samples_df$Value)
samples_df$ImValue <- Im(samples_df$Value)

# Filter and plot
L <- max(samples_df$x1) # Replace with the actual value of L if it's different
ggplot(samples_df %>% filter(x1 < L), aes(x = x1)) +
  geom_line(aes(y = ReValue, color = Type, linetype = "Real part")) +
  geom_line(aes(y = ImValue, color = Type, linetype = "Imag part")) +
  coord_cartesian(xlim = c(0, L)) +
  xlab("x") + ylab("Sample") +
  scale_linetype_manual(values = c("Real part" = "solid", "Imag part" = "dashed"))


