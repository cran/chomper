## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, message=FALSE-------------------------------------------------
library(chomper)

data(italy)

# Select variables for the estimation
#   1. Multinomial
#     SEX: sex
#     ANASC: year of birth
#     CIT: whether Italian or not
#     NASCREG: place of birth
#     STUDIO: education level
#   2. Gaussian (standardized)
#     NETINCOME: net annual income
#     VALABIT: estimated price of residence
discrete_colnames <- c("SEX", "ANASC", "CIT", "NASCREG", "STUDIO")
continuous_colnames <- c("NETINCOME", "VALABIT")
common_fields <- c(discrete_colnames, continuous_colnames)

# Set the types of variables.
# If one type is missing, please set it to 0
discrete_fields <- 1:5
n_discrete_fields <- 5
continuous_fields <- 6:7 # If no continuous fields are used, replace them with 0,
n_continuous_fields <- 2 # and replace this argument with 0

k <- length(italy) # number of files, 2 in this case
n <- numeric(k) # number of records in each file

M <- numeric(n_discrete_fields) # number of levels in each categorical variable
for (l in seq_along(discrete_colnames)) {
  for (i in 1:k) {
    if (i == 1) {
      categories <- italy[[i]][, discrete_colnames[l]]
    } else {
      categories <- c(categories, italy[[i]][, discrete_colnames[l]])
    }
  }
  M[l] <- length(unique(categories))
}

dat <- list()
for (i in 1:k) {
  # Select records from the prespecified region;
  #   IREG: current living region of a respondent
  dat[[i]] <-
    italy[[i]][italy[[i]][, "IREG"] == 7, ]
  n[i] <- nrow(dat[[i]])
}
N <- sum(n) # total number of records
N_max <- N # total number of records, which is used for setting a prior

p <- length(common_fields) # number of variables

# Define the data to link and confirm that each data has a form of matrix
x <- list()
for (i in 1:k) {
  x[[i]] <- as.matrix(dat[[i]][, common_fields])
}

## ----echo=TRUE, message=FALSE-------------------------------------------------
# hyperparameter for categorical fields (Multinomial distribution)
hyper_phi <- rep(2.0, n_discrete_fields)
hyper_tau <- rep(0.01, n_discrete_fields)

hyper_epsilon_discrete <-
  c(rep(0, n_discrete_fields - 1), 1) # hitting range for categorical fields.
# As education level is assumed to have multiple truths,
# we set epsilon of that field as 1.
hyper_epsilon_continuous <-
  rep(0.1, n_continuous_fields) # hitting range for continuous fields.
# Continuous fields are assumed to have multiple truths, we set epsilon as 0.1.

## ----echo=TRUE, message=FALSE-------------------------------------------------
hyper_sigma <- matrix(
  rep(c(0.01, 0.01), n_continuous_fields),
  ncol = 2, byrow = TRUE
) # prior distribution for the variance of Gaussian, setting non-informative

hyper_beta <- matrix(
  rep(c(N_max * 0.1 * 0.01, N_max * 0.1), p),
  ncol = 2, byrow = TRUE
) # prior distribution for distortion rate

## ----echo=TRUE, message=FALSE, eval=FALSE-------------------------------------
#  res <- chomperMCMC(
#    x = x, k = k, n = n, N = N, p = p, M = M,
#    discrete_fields = discrete_fields,
#    continuous_fields = continuous_fields,
#    hyper_beta = hyper_beta,
#    hyper_phi = hyper_phi,
#    hyper_tau = hyper_tau,
#    hyper_epsilon_discrete = hyper_epsilon_discrete,
#    hyper_epsilon_continuous = hyper_epsilon_continuous,
#    hyper_sigma = hyper_sigma,
#    n_burnin = 4000, # number of burn-in samples
#    n_gibbs = 1000, # number of MCMC samples to store
#    n_split_merge = 1000, # number of split and merge procedures in each Gibbs iteration
#    max_time = 86400 # maximum time in second for MCMC
#  )
#  
#  # Example code for using `chomperEVIL`
#  # res <- chomperEVIL(
#  #     x = x, k = k, n = n, N = N, p = p, M = M,
#  #     discrete_fields = discrete_fields,
#  #     continuous_fields = continuous_fields,
#  #     hyper_beta = hyper_beta,
#  #     hyper_phi = hyper_phi,
#  #     hyper_tau = hyper_tau,
#  #     hyper_delta = hyper_delta,
#  #     hyper_sigma = hyper_sigma,
#  #     overlap_prob = 0.5, # probability for the initialization.
#  #     # This probability of records will be initially considered as linked.
#  #     n_parents = 50, # number of parents, N_P
#  #     n_children = 100, # number of offspring, N_O
#  #     tol_cavi = 1e-5, # stopping criterion for a single CAVI
#  #     max_iter_cavi = 10, # maximum number of CAVI update for a single member
#  #     tol_evi = 1e-5, # stopping criterion for an evolution
#  #     max_iter_evi = 50, # maximum number of evolution, N_E
#  #     verbose = 1,
#  #     n_threads = 7 # number of cores to be parallelized
#  # )

## ----echo=TRUE, message=FALSE-------------------------------------------------
library(blink)
library(salso)

for (i in 1:k) { # `k` represents the number of files
  if (i == 1) {
    key_true <- dat[[i]][, "ID"]
  } else {
    key_true <- c(key_true, dat[[i]][, "ID"])
  }
}

# Find the true linkage structure.
# An external R package, `blink`, is used, especially, `blink::links`
truth_binded <- matrix(key_true, nrow = 1)
linkage_structure_true <- links(truth_binded, TRUE, TRUE)
linkage_truth <- matrix(linkage_structure_true)

head(linkage_truth)

## ----echo=FALSE, message=FALSE------------------------------------------------
res <- readRDS("vignette_lambda.rds")
res$lambda <- res$lambda[4001:5000]

## ----echo=TRUE, message=FALSE-------------------------------------------------
# In order to use MCMC samples:
#   1. slice the actual MCMC samples
#   2. reshape a list into a matrix form
lambda_binded <-
  flatten_posterior_samples(
    res$lambda, k, N
  )

psm_lambda <- psm_mcmc(lambda_binded)

# # Approximate posterior distribution is used directly
# # to calculate a posterior similarity matrix
# psm_lambda <- psm_vi(result$nu)

# Calculate a Bayes estimate that minimizes Binder's loss
salso_estimate <- salso(round(psm_lambda, nchar(N) + 1), # rounding is implemented to avoid floating-point error
  loss = binder(),
  maxZealousAttempts = 0, probSequentialAllocation = 1
)

# Convert a Bayes estimate to have an appropriate structure to evaluate the performance
linkage_structure <- list()
for (ll in seq_along(salso_estimate)) {
  linkage_structure[[ll]] <- which(salso_estimate == salso_estimate[ll])
}
linkage_estimation <- matrix(linkage_structure)

## ----echo=TRUE, message=FALSE-------------------------------------------------
# Calculate the performance of the estimation
perf <- performance(linkage_estimation, linkage_truth, N, FALSE)

# Print the performance
performance_matrix <-
  c(perf$fnr, perf$fpr, perf$tp, perf$fp, perf$tn, perf$fn)
performance_matrix <-
  matrix(
    c(2 * performance_matrix[3] /
      (2 * performance_matrix[3] +
        performance_matrix[4] +
        performance_matrix[6]), performance_matrix),
    nrow = 1
  )
colnames(performance_matrix) <- c("F1", "FNR", "FPR", "TP", "FP", "TN", "FN")
rownames(performance_matrix) <- c("Region 7")
print(performance_matrix)

## ----echo=FALSE, message=FALSE------------------------------------------------
res <- readRDS("vignette_lambda.rds")

## ----echo=TRUE, message=FALSE, eval=FALSE-------------------------------------
#  res <- chomperMCMC(
#    x = x, k = k, n = n, N = N, p = p, M = M,
#    discrete_fields = discrete_fields,
#    continuous_fields = continuous_fields,
#    hyper_beta = hyper_beta,
#    hyper_phi = hyper_phi,
#    hyper_tau = hyper_tau,
#    hyper_epsilon_discrete = hyper_epsilon_discrete,
#    hyper_epsilon_continuous = hyper_epsilon_continuous,
#    hyper_sigma = hyper_sigma,
#    n_burnin = 0, # number of burn-in samples
#    n_gibbs = 5000, # number of MCMC samples to store
#    n_split_merge = 1000, # number of split and merge procedures in each Gibbs iteration
#    max_time = 86400 # maximum time in second for MCMC
#  )

## ----echo=TRUE, message=FALSE, warning=FALSE, fig.width=8, fig.height=3, out.width="95%"----
library(ggplot2)
library(patchwork)

unique_entities <- numeric(5000)
for (i in 1:5000) {
  unique_entities[i] <-
    length(unique(c(res$lambda[[i]][[1]], res$lambda[[i]][[2]])))
}

df_trace <- data.frame(idx = 1:5000, unique_entities = unique_entities)

plot_trace <- ggplot(df_trace, aes(x = idx)) +
  geom_line(aes(y = unique_entities)) +
  labs(title = "Traceplot of MCMC Samples", x = NULL, y = NULL) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = .2),
    axis.line.y = element_line(linewidth = .2)
  )

true_unique_entities <- length(unique(key_true))

df_density <- data.frame(unique_entities = unique_entities[4001:5000])

plot_density <- ggplot(df_density) +
  geom_density(aes(x = unique_entities, colour = "Estimate"),
    adjust = 2, key_glyph = "path"
  ) +
  geom_vline(aes(xintercept = length(unique(key_true)), colour = "Truth")) +
  scale_colour_manual(values = c("Estimate" = "black", "Truth" = "red")) +
  xlim(500, 1000) +
  labs(title = "Estimated Number of Unique Entities", x = NULL, y = NULL) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.9),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = .2),
    axis.line.y = element_line(linewidth = .2)
  )

plot_trace + plot_density + plot_layout(ncol = 2)

