library(devtools)
library(MGDrivE)
library(MGDrivE2)

rm(list = ls())
gc()

source("cube.R")
options(scipen = 999)

# external data
b_matrix <- read.csv("b-matrix-coluzii-low.csv")
b_genos <- b_matrix[, 1]


# simulation parameters
tmax <- 365 * 15
dt <- 1
dt_stoch <- 0.01

mating.comp <- c(
  BB = 1.0,
  BH = 0.9,
  BR = 1.0,
  BW = 1.0,
  HH = 0.9,
  HR = 0.9,
  HW = 0.9,
  RR = 1.0,
  RW = 1.0,
  WW = 1.0
)
eta <-
  split(cbind(names(mating.comp), mating.comp), seq(length(mating.comp)))
cube <- cubeTP13(
  gtype = c("BB", "BH", "BR", "BW", "HH", "HR", "HW", "RR", "RW", "WW"),
  p1_germline_M = 0,
  p2_germline_M = 0.9787776157,
  p3_germline_M = 0,
  p1_germline_F = 0,
  p2_germline_F = 0.9850854772,
  p3_germline_F = 0,
  p1_MDHH = 0.0437169674,
  p2_MDHH = 0,
  p3_MDHH = 0.0127526596,
  p1_MDH = 0.0017235774,
  p2_MDH = 0,
  p3_MDH = 0.0006749124,
  fc = c(
    BB = 1.0,
    BH = 0.9,
    BR = 1.0,
    BW = 1.0,
    HH = 0.9,
    HR = 0.9,
    HW = 0.9,
    RR = 1.0,
    RW = 1.0,
    WW = 1.0
  ),
  # genotype-specific modifiers
  eta = eta,
  phi = NULL,
  omega = NULL,
  xiF = NULL,
  xiM = NULL,
  s  = NULL
)

STP_Pfpr <- 0.02
NH <- 1000
theta <- imperial_model_param_list_create(NH = NH)
age_vector <-
  c(0,
    5,
    17,
    40,
    60,
    99)

ft <- 0.02
IRS_cov <- 0.665
LLIN_cov <- 0.62
theta <- add_interventions(theta, IRS_cov, LLIN_cov)

SPN_P <- spn_P_epi_decoupled_node(params = theta, cube = cube)
SPN_T <-
  spn_T_epi_decoupled_node(spn_P = SPN_P,
                           params = theta,
                           cube = cube)

# Stoichiometry matrix
S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)
prev <- STP_Pfpr
eir <- convert_prevalence_to_eir(prev, age_vector, ft, theta)

# Hazards

eqm <-
  equilibrium_Imperial_decoupled(age_vector, ft, eir, theta, cube, SPN_P)

# extract updated theta and full set of initial conditions
theta <- eqm$theta
cube <- eqm$cube
initialCons <- eqm$initialCons

# set up time varying stuff
# set up time varying carrying capacity
carry <-
  read.csv2(file = "./data/carrying_capacity_stp.csv",
            sep = ",",
            stringsAsFactors = FALSE)
K_mean <- mean(as.numeric(carry$K))
K_eq <- theta$K
adjustment <- K_eq / K_mean

K_ts <- as.numeric(carry$K)
K_ts <- K_ts * adjustment

step_K <-
  stats::stepfun(
    x = 1:nrow(carry),
    y = c(K_ts[1], K_ts),
    f = 0,
    right = FALSE
  )
theta$K <- c(step_K)
b_names <- colnames(b_matrix)[-1]
thresh <- 5
b_fitted <- b_matrix[, b_names[thresh]]
names(b_fitted) <- b_genos
theta$b0 <- b_fitted

make_larvae_mort_haz_log_inhom <-
  function(trans,
           u,
           l_ix,
           node,
           cube,
           params,
           exact = TRUE,
           tol = 1e-8) {
    # rate constants
    muL <- params$muL
    K <- params$K[[node]]
    if (typeof(K) != "closure") {
      stop(
        "Inhomogeneous hazard 'make_larvae_mort_haz_log', ",
        "value 'K' in 'params' list needs to be a function"
      )
    }
    
    
    # which places have input arcs to this transition
    s <- trans$s
    
    # weights of those arcs
    w <- trans$s_w
    
    # assign here so that each newly generated closure has the right indices
    l_ix <- l_ix
    
    # return the hazard function
    if (exact) {
      # EXACT hazards (check enabling degree: for discrete simulation only)
      return(function(t, M) {
        if (w <= M[s]) {
          L <- sum(M[l_ix])
          return(muL * (1 + (L / K(t))) * M[s])
          # return(muL*(1 + (L/K_mean))*M[s])
        } else {
          return(0)
        }
      })
      
    } else {
      # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
      return(function(t, M) {
        # get total males
        L <- sum(M[l_ix])
        haz <- muL * (1 + (L / K(t))) * M[s]
        # haz <- muL*(1 + (L/K_mean))*M[s]
        # check and return
        if (haz < tol) {
          return(0)
        } else {
          return(haz)
        }
      })
      
    }
    # end of function
  }

# make hazards by hand
make_hazards <-
  function(spn_P,
           spn_T,
           cube,
           par,
           log_dd = TRUE,
           exact = TRUE,
           tol = 1e-12,
           verbose = TRUE) {
    if (tol > 1e-6 & !exact) {
      cat(
        "warning: hazard function tolerance ",
        tol,
        " is large; consider tolerance < 1e-6 for sufficient accuracy\n"
      )
    }
    
    if (log_dd) {
      if (!("K" %in% names(par))) {
        stop(
          "if using logistic (carrying capacity) based density-dependent larval mortality, please specify parameter 'K' in par"
        )
      }
    } else {
      if (!("gamma" %in% names(par))) {
        stop(
          "if using Lotka-Volterra based density-dependent larval mortality, please specify parameter 'gamma' in par"
        )
      }
    }
    
    # transitions and places
    v <- spn_T$v
    u <- spn_P$u
    
    n <- length(v)
    if (verbose) {
      pb <- txtProgressBar(min = 1,
                           max = n,
                           style = 3)
      pp <- 1
    }
    
    # the hazard functions
    h <- vector("list", n)
    h <- setNames(h, v)
    
    # get male and larvae indices
    l_ix <- as.vector(spn_P$ix[[1]]$larvae)
    m_ix <- spn_P$ix[[1]]$males
    
    # human indices
    h_ix <- spn_P$ix[[1]]$humans
    
    cat(" --- generating hazard functions for SPN --- \n")
    
    # make the hazards
    for (t in 1:n) {
      type <- spn_T$T[[t]]$class
      
      # make the correct type of hazard
      
      # MOSQUITO HAZARDS
      if (type == "oviposit") {
        h[[t]] <-
          MGDrivE2:::make_oviposit_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "egg_adv") {
        h[[t]] <-
          MGDrivE2:::make_egg_adv_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "egg_mort") {
        h[[t]] <-
          MGDrivE2:::make_egg_mort_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "larvae_adv") {
        h[[t]] <-
          MGDrivE2:::make_larvae_adv_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
        # INHOMOGENEOUS
      } else if (type == "larvae_mort") {
        h[[t]] <-
          make_larvae_mort_haz_log_inhom(
            t = spn_T$T[[t]],
            u = u,
            l_ix = l_ix,
            node = 1,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "pupae_adv") {
        h[[t]] <-
          MGDrivE2:::make_pupae_adv_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "pupae_mort") {
        h[[t]] <-
          MGDrivE2:::make_pupae_mort_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "pupae_2m") {
        h[[t]] <-
          MGDrivE2:::make_pupae_2male_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "pupae_2f") {
        h[[t]] <-
          MGDrivE2:::make_pupae_2female_haz(
            t = spn_T$T[[t]],
            u = u,
            m_ix = m_ix,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "pupae_2unmated") {
        h[[t]] <-
          MGDrivE2:::make_pupae_2unmated_haz(
            t = spn_T$T[[t]],
            u = u,
            m_ix = m_ix,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "female_unmated_mate") {
        h[[t]] <-
          MGDrivE2:::make_unmated_2female_haz(
            t = spn_T$T[[t]],
            u = u,
            m_ix = m_ix,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
        # INHOMOGENEOUS
      } else if (type == "male_mort") {
        h[[t]] <-
          MGDrivE2:::make_male_mort_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
        # INHOMOGENEOUS
      } else if (type %in% c("female_mort", "female_unmated_mort")) {
        h[[t]] <-
          MGDrivE2:::make_female_mort_haz(
            t = spn_T$T[[t]],
            u = u,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "female_inf") {
        h[[t]] <-
          MGDrivE2:::make_female_inf_epi_haz_decoupled_Imperial(
            t = spn_T$T[[t]],
            u = u,
            h_ix = h_ix,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "female_eip") {
        h[[t]] <-
          MGDrivE2:::make_female_eip_epi_haz(
            t = spn_T$T[[t]],
            u = u,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "female_inc") {
        # can reuse above hazard because transition hazard is the same
        h[[t]] <-
          MGDrivE2:::make_female_eip_epi_haz(
            t = spn_T$T[[t]],
            u = u,
            par = par,
            exact = exact,
            tol = tol
          )
        # HUMAN HAZARDS
      } else if (type == "H_birth") {
        h[[t]] <-
          MGDrivE2:::make_human_birth_sis_haz(
            t = spn_T$T[[t]],
            u = u,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "H_mort") {
        h[[t]] <-
          MGDrivE2:::make_human_death_sis_haz(
            t = spn_T$T[[t]],
            u = u,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "H_infection") {
        h[[t]] <-
          MGDrivE2:::make_human_inf_sis_haz(
            t = spn_T$T[[t]],
            u = u,
            h_ix = h_ix,
            cube = cube,
            par = par,
            exact = exact,
            tol = tol
          )
      } else if (type == "H_recovery") {
        h[[t]] <-
          MGDrivE2:::make_human_rec_sis_haz(
            t = spn_T$T[[t]],
            u = u,
            par = par,
            exact = exact,
            tol = tol
          )
      } else {
        stop(paste0(
          "error in making hazard function for unknown class type: ",
          type
        ))
      }
      
      if (verbose) {
        setTxtProgressBar(pb, t)
      }
    }
    
    if (verbose) {
      close(pb)
    }
    
    cat(" --- done generating hazard functions for SPN --- \n")
    
    return(list("hazards" = h, "flag" = exact))
  }

# hazard vector
hazards <- make_hazards(
  spn_P = SPN_P,
  spn_T = SPN_T,
  cube = cube,
  par = theta,
  log_dd = TRUE,
  exact = TRUE,
  verbose = TRUE
)

# release strategy
r_times <- seq(from = 365 * 10,
               length.out = 8,
               by = 7)
r_size <- 1e4 * 2

events <- data.frame(
  "var" = paste0("M_", cube$releaseType),
  "time" = r_times,
  "value" = r_size,
  "method" = "add",
  stringsAsFactors = FALSE
)


tau_out <- sim_trajectory_R_decoupled(
  x0 = initialCons$M0,
  h0 = initialCons$H,
  SPN_P = SPN_P,
  theta = theta,
  tmax = tmax,
  inf_labels = SPN_T$inf_labels,
  dt = dt,
  dt_stoch = dt_stoch,
  S = S,
  hazards = hazards,
  sampler = "tau-decoupled",
  events = events,
  verbose = FALSE,
  human_ode = "Imperial",
  cube = cube,
  maxhaz = 1e12
)
