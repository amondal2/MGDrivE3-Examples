library(MGDrivE)
library(devtools)

setwd("~/MGDrivE/MGDrivE-2/MGDrivE2")
load_all(".")
setwd("~/MGDrivE/Main/SoftwarePaper3")

rm(list=ls()); gc()

cube <- MGDrivE::cubeHoming1RA()
cube <- MGDrivE::cubeRIDL()
# entomological parameters
theta <- list(
  qE = 1 / 4,
  nE = 2,
  qL = 1 / 3,
  nL = 3,
  qP = 1 / 6,
  nP = 2,
  muE = 0.05,
  muL = 0.05,
  muP = 0.05,
  muF = 0.09,
  muM = 0.09,
  beta = 16,
  nu = 1 / (4 / 24)
)

# simulation parameters
tmax <- 365 * 5
dt <- 1

# get number of nodes from PTS file
pts <- read.csv('./data/STP_grid_Sites.csv')
move_probs_raw <-
  read.csv('./data/STP_grid_Migration.csv', header = F)

# create binary adjacency matrix from rates
adj <- move_probs_raw
adj[adj > 0] <- 1
diag(adj) <- 0
adj <- as.matrix(adj, sparse=T)
adj <- as(adj, 'ngCMatrix')
n <- nrow(pts)
adj <- unname(adj)

# Places and transitions
SPN_P <-
  spn_P_lifecycle_network(num_nodes = n,
                          params = theta,
                          cube = cube)
SPN_T <-
  spn_T_lifecycle_network(
    spn_P = SPN_P,
    params = theta,
    cube = cube,
    m_move = adj
  )


# Stoichiometry matrix
S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)

# now that we have a network size, set adult females in each node
NF <- rep(x = 1000, times = n)

# calculate equilibrium and setup initial conditions
#  outputs required parameters in the named list "params"
#  outputs intial equilibrium for adv users, "init
#  outputs properly filled initial markings, "M0"
initialCons <-
  equilibrium_lifeycle(
    params = theta,
    NF = NF,
    phi = 0.5,
    log_dd = TRUE,
    spn_P = SPN_P,
    cube = cube
  )

# extract indices of trap nodes
types <- pts$TrapsType
idx <- which(types == 0)

move_probs_raw <- Matrix::as.matrix(move_probs_raw, sparse = T)
move_probs <- move_probs_raw
# remove probabilities associated with traps
move_probs <-
  move_probs[1:(nrow(move_probs) - length(idx)), 1:(ncol(move_probs) - length(idx))]

# solve Prob = 1-exp(-Rate)
P <- diag(move_probs)
move_rates <- log(1 / (1 - P))

# set move_rates for traps to 0
move_rates <- c(move_rates, rep(0, length(idx)))

# renormalize matrix
move_probs_renorm <- move_probs_raw
diag(move_probs_renorm) <- 0
move_probs_renorm <-
  t(scale(
    t(move_probs_renorm),
    center = F,
    scale = rowSums(move_probs_renorm)
  ))

# replace NA with 0 for trap rows - these rows don't matter since the movement rate is 0 anyway
move_probs_renorm[is.na(move_probs_renorm)] <- 0

# put rates and probs into the parameter list
initialCons$params$mosquito_move_rates <- move_rates
initialCons$params$mosquito_move_probs <- move_probs_renorm

# approximate hazards for continous approximation
approx_hazards <-
  spn_hazards(
    spn_P = SPN_P,
    spn_T = SPN_T,
    cube = cube,
    params = initialCons$params,
    type = "life",
    log_dd = TRUE,
    exact = FALSE,
    tol = 1e-6,
    verbose = FALSE
  )

# release strategy
r_times <- seq(from = 365 * 2,
               length.out = 8,
               by = 7)
r_size <- 1e4 * 2
events <- data.frame(
  "var" = paste0("M_", cube$releaseType, "_1"),
  "time" = r_times,
  "value" = r_size,
  "method" = "add",
  stringsAsFactors = FALSE
)

dt_stoch <- 0.05

# tau leaping simulation
PTS_out <-
  sim_trajectory_R(
    x0 = initialCons$M0,
    tmax = tmax,
    dt = dt,
    dt_stoch = dt_stoch,
    S = S,
    hazards = approx_hazards,
    sampler = "tau",
    events = events,
    verbose = TRUE,
    maxhaz=1e14
  )

# summarize females/males by genotype
PTS_female <- summarize_females(out = PTS_out$state, spn_P = SPN_P)
PTS_male <- summarize_males(out = PTS_out$state)

# add sex for plotting
PTS_female$sex <- "Female"
PTS_male$sex <- "Male"

# plot
ggplot(data = rbind(PTS_female, PTS_male)) +
  geom_line(aes(x = time, y = value, color = genotype)) +
  facet_grid(node ~ sex, scales = "fixed") +
  theme_bw() +
  ggtitle("SPN: Tau-leaping Approximation")
