source("S2C.R")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
logkappa <- 0

#Defining vectors for anisotopy
v1 <- c(0.001,0) #Isotropic matern
v2 <- c(1,0)     #Anisotropic
v3 <- c(-1,0)

#Defining spectrum for different values of v

S_fun1 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v1, omega=omega)
S_fun2 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v2, omega=omega)
S_fun3 <- function(omega) S_aniso_matrix(logkappa=logkappa, v=v3, omega=omega)

#Defining discretisation points
n <- c(512, 512)
L <- c(2, 2)
h <- L / n
x_ <- make_x(n, L)
omega_ <- make_omega(n, L)
omega <- as.matrix(expand.grid(omega_))

#Defining spectrum at each frequency
S_1 <- S_fun1(omega = omega)
S_2 <- S_fun2(omega = omega)
S_3 <- S_fun3(omega = omega)

#Defining data frame for spectrum and covariance for each field
S_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(omega_))),
    data.frame(
      s_1 = S_1,
      s_2 = S_2,
      s_3 = S_3
    )
  ) %>%
  pivot_longer(
    cols = c(s_1, s_2, s_3),
    names_to = "Type", values_to = "Value"
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
    names_to = "Type", values_to = "Value"
  )


#Plotting spectrum and covariance
ggplot(S_df) +
  geom_raster(aes(w1, w2, fill = Value)) +
  scale_fill_gradientn(colours=brewer.pal(7,"Reds"))+
  # scale_fill_viridis(option="inferno")+
  coord_equal()+
  facet_wrap(vars(Type),nrow=1) +
  xlab("omega1") + ylab("omega2") +
  labs(fill = "Spectrum")
ggplot(C_df) +
  geom_raster(aes(x1, x2, fill = Value)) +
  scale_fill_gradientn(colours=brewer.pal(7,"Reds"))+
  coord_equal()+
  facet_wrap(vars(Type),nrow =1) +
  xlab("u1") + ylab("u2") +
  labs(fill = "Covariance")



#Sampling
#Calculate samples
sample_1 <- S2sample(S=S_1, dim=n, h=h, seed = NULL, conjugate = FALSE)
sample_2 <- S2sample(S=S_2, dim=n, h=h, seed = NULL, conjugate = FALSE)
sample_3 <- S2sample(S=S_3, dim=n, h=h, seed = NULL, conjugate = FALSE)

#Calculate spatial grid and assign samples to their location
x_ <- make_x_sampling(n, h)
samples_df <-
  cbind(
    as.data.frame(as.matrix(expand.grid(x_))),
    data.frame(
      u1 = as.vector(sample_1),u2 = as.vector(sample_2),u3 = as.vector(sample_3))
    )%>%
  pivot_longer(
    cols = c("u1", "u2", "u3"),
    names_to = "Type", values_to = "Value"
  )

#Plot samples over grid
ggplot(samples_df) +
  geom_raster(aes(x1, x2, fill = Re(Value))) +
  #scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"))+
  scale_fill_gradientn(colours=brewer.pal(7,"Reds"))+
  coord_equal()+
  facet_wrap(vars(Type),nrow=1) +
  xlab("x1") + ylab("x2") +
  labs(fill = "Field")
