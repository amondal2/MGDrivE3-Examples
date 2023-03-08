#' Generate and Modify Default Genotype-specific Parameters
#'
#' This is an internal function for cubes.
#'
#' @param gtype character vector of genotypes
#' @param eta genotype-specific mating fitness, handles assortative mating as well
#' @param phi genotype-specific sex ratio at emergence
#' @param omega genotype-specific multiplicative modifier of adult mortality
#' @param xiF genotype-specific female pupatory success
#' @param xiM genotype-specific male pupatory success
#' @param s genotype-specific fractional reduction(increase) in fertility
#'
cubeModifiers <- function(gtype, eta = NULL, phi = NULL, omega = NULL,
                          xiF = NULL, xiM = NULL, s = NULL){

  # check all numeric arguments to have proper bounds
#  if(any(eta < 0)) stop("eta values must be positive [X>0]")
  if(any(phi > 1) || any(phi < 0)) stop("phi values must be between [0,1]")
  if(any(omega < 0)) stop("omega values must be positive [X>0]")
  if(any(xiF > 1) || any(xiF < 0)) stop("xiF values must be between [0,1]")
  if(any(xiM > 1) || any(xiM < 0)) stop("xiM values must be between [0,1]")
  if(any(s < 0)) stop("s values must be positive [X>0]")

  # add parameters in the right place
  set <- function(gtype,vector,vectorNew){
    if(is.null(vectorNew)){
      return(vector)
    } else {
      if(any(!names(vectorNew) %in% gtype)){
        stop("genotype(s) do not match genotypes in cube; please check names of input genotype-specific parameters")
      }
      if(length(vectorNew)==1L){
        if(is.null(names(vectorNew))){
          vector[1:length(vector)] = unname(vectorNew)
        } else {
          vector[names(vectorNew)] = vectorNew
        }
      } else {
        vector[names(vectorNew)] = vectorNew
      }
      return(vector)
    }
  }

  ## genotype-specific modifiers
  size = length(gtype)
#  etaN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific mating fitness
  phiN = setNames(object = rep.int(x = 0.5, times = size), nm = gtype)     # genotype-specific sex ratio at emergence
  omegaN = setNames(object = rep.int(x = 1, times = size), nm = gtype)    # genotype-specific multiplicative modifier of adult mortality
  xiFN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific female pupatory success
  xiMN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific male pupatory success
  sN = setNames(object = rep.int(x = 1, times = size), nm = gtype)        # genotype-specific fractional reduction(increase) in fertility

  # add the user parameters
#  etaO = set(gtype,etaN,eta)
  phiO = set(gtype,phiN,phi)
  omegaO = set(gtype,omegaN,omega)
  xiFO = set(gtype,xiFN,xiF)
  xiMO = set(gtype,xiMN,xiM)
  sO = set(gtype,sN,s)

  ##########
  # generate eta matrix
  ##########
  if(is.null(eta)){
    # default behaviour
    # This generates a mating matrix with all mating equal
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))

  } else if(all(lengths(eta) == 2) ) {
    # set mating weights over specific male genotypes and all female genotypes
    # takes input list of mate genotype, and weight

    # pull out values
    mGeno <- unlist(lapply(X = eta, FUN = "[[", 1))
    weight <- as.numeric(lapply(X = eta, FUN = "[[", 2))

    # check for correctly specified genotypes
    if(!all(mGeno %in% gtype)){
      stop(paste0("Genotypes specified for male mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(any(weight < 0)){
      stop(paste0("Mating ability weights must be greater than or equal to zero."))
    }

    # setup/fill mating matrix
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))
    etaO[ ,mGeno] <- rep(x = weight, each = size)

  } else if(all(lengths(eta) == 3) ) {
    # mating ability, depends on genotype of male and female
    # takes input list of female genotype, mate genotype, and weight

    # pull out values
    fGeno <- unlist(lapply(X = eta, FUN = "[[", 1))
    mGeno <- unlist(lapply(X = eta, FUN = "[[", 2))
    weight <- as.numeric(lapply(X = eta, FUN = "[[", 3))

    # check for correctly specified genotypes
    if(!all(fGeno %in% gtype)){
      stop(paste0("Genotypes specified for female mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(!all(mGeno %in% gtype)){
      stop(paste0("Genotypes specified for male mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(any(weight < 0)){
      stop(paste0("Mating ability weights must be greater than or equal to zero."))
    }

    # setup/fill mating matrix
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))
    etaO[cbind(fGeno,mGeno)] <- weight

  } else {
    stop("eta is incorrectly formatted.\n
         Check the examples or set as NULL for default behaviour.")
  }

  # return named list
  return(
    list(
      eta=etaO,
      phi=phiO,
      omega=omegaO,
      xiF=xiFO,
      xiM=xiMO,
      s=sO
    )
  )
}

###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   TP13
#
###############################################################################

#' Inheritance Cube: Dual effector population modification gene-drive strains of
#' the African malaria mosquitoes, Anopheles gambiae and Anopheles coluzzii
#'
#' Two strains have a coupled Cas9/gRNA-based, autonomous gene-drive system with
#' anti-Plasmodium effector genes comprising dual single-chain variable fragment
#' monoclonal antibodies (scFvs) targeting the parasite ookinetes and sporozoites
#' in the African malaria mosquitoes Anopheles gambiae (AgTP13) and An. coluzzii
#' (AcTP13).
#'
#' @param gtype         Vector of genotypes, with the wild-type in the *LAST* position
#' @param p1_germline_M Probability of homing 1 - male
#' @param p2_germline_M Probability of homing 2 - male
#' @param p3_germline_M Probability of homing 3 - male
#' @param p1_germline_F Probability of homing 1 - female
#' @param p2_germline_F Probability of homing 2 - female
#' @param p3_germline_F Probability of homing 3 - female
#' @param p1_MDHH       Probability of maternal deposition 1 - homozygous mother HH
#' @param p2_MDHH       Probability of maternal deposition 2 - homozygous mother HH
#' @param p3_MDHH       Probability of maternal deposition 3 - homozygous mother HH
#' @param p1_MDH        Probability of maternal deposition 1 - hemizygous mother H_
#' @param p2_MDH        Probability of maternal deposition 2 - hemizygous mother H_
#' @param p3_MDH        Probability of maternal deposition 3 - hemizygous mother H_
#' @param fc            Genotype-specific fitness cost
#' @param eta           Genotype-specific mating fitness
#'
#' @return            Named list containing the inheritance cube, transition
#' matrix, genotypes, wild-type allele and all genotype-specific parameters.

cubeTP13 <-
  function(gtype        = c("BB", "BH", "BR", "BW", "HH", "HR", "HW", "RR", "RW", "WW"),
           p1_germline_M = 0,
           p2_germline_M = 0.997,
           p3_germline_M = 0,
           p1_germline_F = 0.004,
           p2_germline_F = 0.924 / (1 - 0.004),
           p3_germline_F = 0,
           p1_MDHH       = 0.0248,
           p2_MDHH       = 0,
           p3_MDHH       = 0.0220 / (1 - 0.0248),
           p1_MDH        = 0.0012,
           p2_MDH        = 0,
           p3_MDH        = 0.0026 / (1 - 0.0012),
           fc            = c(
             BB = 1.0,
             BH = 1.0,
             BR = 1.0,
             BW = 1.0,
             HH = 1.0,
             HR = 1.0,
             HW = 1.0,
             RR = 1.0,
             RW = 1.0,
             WW = 1.0
           ),
           eta           = NULL,
           phi           = NULL,
           omega         = NULL,
           xiF           = NULL,
           xiM           = NULL,
           s             = NULL) {
    # gtype: sorted by alphabetical order
    gtype = sort(vapply(strsplit(gtype, split = ''),
                        function(x) {
                          paste0(sort(x, method = 'auto'), collapse = '')
                        },
                        FUN.VALUE = character(1)))
    
    # size: number of genotypes
    size = length(gtype)
    
    # N: number of different possible alleles
    N = length(unique(unlist(strsplit(gtype, split = ''))))
    
    # Matrix Dimensions Key: [femaleGenotype, maleGenotype, offspringGenotype]
    tMatrix = array(
      data = 0,
      dim = c(size, size, size),
      dimnames = list(gtype, gtype, gtype)
    )
    
    #############################################################################
    # safety checks
    #############################################################################
    
    if (any(c(p1_germline_M, p2_germline_M, p3_germline_M, p1_germline_F, p2_germline_F, p3_germline_F, p1_MDHH, p2_MDHH, p3_MDHH, p1_MDH, p2_MDH, p3_MDH) >
            1) ||
        any(c(p1_germline_M, p2_germline_M, p3_germline_M, p1_germline_F, p2_germline_F, p3_germline_F, p1_MDHH, p2_MDHH, p3_MDHH, p1_MDH, p2_MDH, p3_MDH) <
            0) ||
        is.null(c(p1_germline_M)) ||
        is.null(c(p2_germline_M)) ||
        is.null(c(p3_germline_M)) ||
        is.null(c(p1_germline_F)) ||
        is.null(c(p2_germline_F)) ||
        is.null(c(p3_germline_F)) ||
        is.null(c(p1_MDHH)) ||
        is.null(c(p2_MDHH)) ||
        is.null(c(p3_MDHH)) ||
        is.null(c(p1_MDH)) ||
        is.null(c(p2_MDH)) ||
        is.null(c(p3_MDH))) {
      stop("
         Probabilities must be defined between 0 and 1.
         ")
    }
    
    #############################################################################
    # crosses
    #############################################################################
    
    # mating competition added via cubeModifiers
    # eta: genotype-specific mating fitness
    
    # homing probabilities in germline of converting W into B, H, R or remaining W
    ## male
    pB_M = p1_germline_M
    pH_M = (1 - p1_germline_M) * p2_germline_M
    pR_M = (1 - p1_germline_M) * (1 - p2_germline_M) * p3_germline_M
    pW_M = (1 - p1_germline_M) * (1 - p2_germline_M) * (1 - p3_germline_M)
    ## female
    pB_F = p1_germline_F
    pH_F = (1 - p1_germline_F) * p2_germline_F
    pR_F = (1 - p1_germline_F) * (1 - p2_germline_F) * p3_germline_F
    pW_F = (1 - p1_germline_F) * (1 - p2_germline_F) * (1 - p3_germline_F)
    
    # aux vector to count genotypes of offspring
    testVec = setNames(object = numeric(size), nm = gtype)
    
    # fill tMatrix with probabilities following Mendelian inheritance after homing
    for (i in gtype) {
      for (j in gtype) {
        # male and females gametes
        male.gamete   = strsplit(j, split = '')[[1]]
        female.gamete = strsplit(i, split = '')[[1]]
        
        # homing
        if ((female.gamete[1] == "H" &
             female.gamete[2] == "W") &
            (male.gamete[1] == "H" & male.gamete[2] == "W")) {
          female.gamete            = c("H", "B", "H", "R", "W")
          female.gamete.proportion = c(0.5, pB_F * 0.5, pH_F * 0.5, pR_F *
                                         0.5, pW_F * 0.5)
          male.gamete              = c("H", "B", "H", "R", "W")
          male.gamete.proportion   = c(0.5, pB_M * 0.5, pH_M * 0.5, pR_M *
                                         0.5, pW_M * 0.5)
          offspring                = as.vector(outer(male.gamete, female.gamete, paste0, sep =
                                                       ''))
          offspring.proportion     = as.vector(male.gamete.proportion %o% female.gamete.proportion)
          
        } else if (female.gamete[1] == "H" &
                   female.gamete[2] == "W") {
          female.gamete            = c("H", "B", "H", "R", "W")
          female.gamete.proportion = c(0.5, pB_F * 0.5, pH_F * 0.5, pR_F *
                                         0.5, pW_F * 0.5)
          male.gamete.proportion   = c(0.5, 0.5)
          offspring                = as.vector(outer(male.gamete, female.gamete, paste0, sep =
                                                       ''))
          offspring.proportion     = as.vector(male.gamete.proportion %o% female.gamete.proportion)
          
        } else if (male.gamete[1] == "H" & male.gamete[2] == "W") {
          male.gamete              = c("H", "B", "H", "R", "W")
          male.gamete.proportion   = c(0.5, pB_M * 0.5, pH_M * 0.5, pR_M *
                                         0.5, pW_M * 0.5)
          female.gamete.proportion = c(0.5, 0.5)
          offspring                = as.vector(outer(male.gamete, female.gamete, paste0, sep =
                                                       ''))
          offspring.proportion     = as.vector(male.gamete.proportion %o% female.gamete.proportion)
          
        } else {
          male.gamete.proportion   = c(0.5, 0.5)
          female.gamete.proportion = c(0.5, 0.5)
          offspring                = as.vector(outer(male.gamete, female.gamete, paste0, sep =
                                                       ''))
          offspring.proportion     = as.vector(male.gamete.proportion %o% female.gamete.proportion)
          
        }
        
        # reorder all offspring alleles to match allele order in gtype
        offspring = vapply(strsplit(offspring, split = ''),
                           function(x) {
                             paste0(sort(x, method = 'auto'), collapse = '')
                           },
                           FUN.VALUE = character(1))
        
        # count genotypes of offspring, order according to gtype
        for (o in 1:length(offspring)) {
          testVec[offspring[o]] = testVec[offspring[o]] + offspring.proportion[o]
        }
        
        #store after normalizing
        tMatrix[i, j, ] = testVec / sum(testVec)
        testVec[]      = 0
      }
    }
    
    #############################################################################
    # maternal deposition
    #############################################################################
    
    # parameters: p1_MDHH, p2_MDHH, p3_MDHH, p1_MDH ,p2_MDH ,p3_MDH
    
    # saving temp variable
    tMatrix.temp = tMatrix
    
    # homozygous mothers (two H alleles) ----
    motherHH = (gtype == "HH")
    
    # probabilities of converting W into B, H, R or remaining W
    pB_MDHH = p1_MDHH
    pH_MDHH = (1 - p1_MDHH) * p2_MDHH
    pR_MDHH = (1 - p1_MDHH) * (1 - p2_MDHH) * p3_MDHH
    pW_MDHH = (1 - p1_MDHH) * (1 - p2_MDHH) * (1 - p3_MDHH)
    
    # zygotes with one W alleles
    tMatrix[motherHH, , "BB"] = tMatrix[motherHH, , "BB"] + pB_MDHH * tMatrix.temp[motherHH, , "BW"]
    tMatrix[motherHH, , "BH"] = tMatrix[motherHH, , "BH"] + pH_MDHH * tMatrix.temp[motherHH, , "BW"]
    tMatrix[motherHH, , "BR"] = tMatrix[motherHH, , "BR"] + pR_MDHH * tMatrix.temp[motherHH, , "BW"]
    tMatrix[motherHH, , "BW"] = tMatrix[motherHH, , "BW"] - (1 - pW_MDHH) *
      tMatrix.temp[motherHH, , "BW"]
    tMatrix[motherHH, , "BH"] = tMatrix[motherHH, , "BH"] + pB_MDHH * tMatrix.temp[motherHH, , "HW"]
    tMatrix[motherHH, , "HH"] = tMatrix[motherHH, , "HH"] + pH_MDHH * tMatrix.temp[motherHH, , "HW"]
    tMatrix[motherHH, , "HR"] = tMatrix[motherHH, , "HR"] + pR_MDHH * tMatrix.temp[motherHH, , "HW"]
    tMatrix[motherHH, , "HW"] = tMatrix[motherHH, , "HW"] - (1 - pW_MDHH) *
      tMatrix.temp[motherHH, , "HW"]
    tMatrix[motherHH, , "BR"] = tMatrix[motherHH, , "BR"] + pB_MDHH * tMatrix.temp[motherHH, , "RW"]
    tMatrix[motherHH, , "HR"] = tMatrix[motherHH, , "HR"] + pH_MDHH * tMatrix.temp[motherHH, , "RW"]
    tMatrix[motherHH, , "RR"] = tMatrix[motherHH, , "RR"] + pR_MDHH * tMatrix.temp[motherHH, , "RW"]
    tMatrix[motherHH, , "RW"] = tMatrix[motherHH, , "RW"] - (1 - pW_MDHH) *
      tMatrix.temp[motherHH, , "RW"]
    
    # zygotes with two W alleles
    p_MDHH_WW_BB = (factorial(2) / (factorial(2) * factorial(0) * factorial(0) *
                                      factorial(0))) * (pB_MDHH ^ 2 * pH_MDHH ^ 0 * pR_MDHH ^ 0 * pW_MDHH ^ 0)
    p_MDHH_WW_HH = (factorial(2) / (factorial(0) * factorial(2) * factorial(0) *
                                      factorial(0))) * (pB_MDHH ^ 0 * pH_MDHH ^ 2 * pR_MDHH ^ 0 * pW_MDHH ^ 0)
    p_MDHH_WW_RR = (factorial(2) / (factorial(0) * factorial(0) * factorial(2) *
                                      factorial(0))) * (pB_MDHH ^ 0 * pH_MDHH ^ 0 * pR_MDHH ^ 2 * pW_MDHH ^ 0)
    p_MDHH_WW_WW = (factorial(2) / (factorial(0) * factorial(0) * factorial(0) *
                                      factorial(2))) * (pB_MDHH ^ 0 * pH_MDHH ^ 0 * pR_MDHH ^ 0 * pW_MDHH ^ 2)
    p_MDHH_WW_BH = (factorial(2) / (factorial(1) * factorial(1) * factorial(0) *
                                      factorial(0))) * (pB_MDHH ^ 1 * pH_MDHH ^ 1 * pR_MDHH ^ 0 * pW_MDHH ^ 0)
    p_MDHH_WW_BR = (factorial(2) / (factorial(1) * factorial(0) * factorial(1) *
                                      factorial(0))) * (pB_MDHH ^ 1 * pH_MDHH ^ 0 * pR_MDHH ^ 1 * pW_MDHH ^ 0)
    p_MDHH_WW_BW = (factorial(2) / (factorial(1) * factorial(0) * factorial(0) *
                                      factorial(1))) * (pB_MDHH ^ 1 * pH_MDHH ^ 0 * pR_MDHH ^ 0 * pW_MDHH ^ 1)
    p_MDHH_WW_HR = (factorial(2) / (factorial(0) * factorial(1) * factorial(1) *
                                      factorial(0))) * (pB_MDHH ^ 0 * pH_MDHH ^ 1 * pR_MDHH ^ 1 * pW_MDHH ^ 0)
    p_MDHH_WW_HW = (factorial(2) / (factorial(0) * factorial(1) * factorial(0) *
                                      factorial(1))) * (pB_MDHH ^ 0 * pH_MDHH ^ 1 * pR_MDHH ^ 0 * pW_MDHH ^ 1)
    p_MDHH_WW_RW = (factorial(2) / (factorial(0) * factorial(0) * factorial(1) *
                                      factorial(1))) * (pB_MDHH ^ 0 * pH_MDHH ^ 0 * pR_MDHH ^ 1 * pW_MDHH ^ 1)
    
    tMatrix[motherHH, , "BB"] = tMatrix[motherHH, , "BB"] + p_MDHH_WW_BB * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "HH"] = tMatrix[motherHH, , "HH"] + p_MDHH_WW_HH * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "RR"] = tMatrix[motherHH, , "RR"] + p_MDHH_WW_RR * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "WW"] = tMatrix[motherHH, , "WW"] - (1 - p_MDHH_WW_WW) * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "BH"] = tMatrix[motherHH, , "BH"] + p_MDHH_WW_BH * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "BR"] = tMatrix[motherHH, , "BR"] + p_MDHH_WW_BR * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "BW"] = tMatrix[motherHH, , "BW"] + p_MDHH_WW_BW * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "HR"] = tMatrix[motherHH, , "HR"] + p_MDHH_WW_HR * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "HW"] = tMatrix[motherHH, , "HW"] + p_MDHH_WW_HW * tMatrix.temp[motherHH, , "WW"]
    tMatrix[motherHH, , "RW"] = tMatrix[motherHH, , "RW"] + p_MDHH_WW_RW * tMatrix.temp[motherHH, , "WW"]
    
    # hemizygous mothers (only one H allele) ----
    motherH  <- (unlist(lapply(strsplit(gtype, split = ''),
                               function(x) {
                                 any(x == "H")
                               })) & motherHH == FALSE)
    
    # probabilities of converting W into B, H, R or remaining W
    pB_MDH = p1_MDH
    pH_MDH = (1 - p1_MDH) * p2_MDH
    pR_MDH = (1 - p1_MDH) * (1 - p2_MDH) * p3_MDH
    pW_MDH = (1 - p1_MDH) * (1 - p2_MDH) * (1 - p3_MDH)
    
    # probabilities of converting W into B, H, R or remaining W
    pB_MDH = p1_MDH
    pH_MDH = (1 - p1_MDH) * p2_MDH
    pR_MDH = (1 - p1_MDH) * (1 - p2_MDH) * p3_MDH
    pW_MDH = (1 - p1_MDH) * (1 - p2_MDH) * (1 - p3_MDH)
    
    # zygotes with one W alleles
    tMatrix[motherH, , "BB"] = tMatrix[motherH, , "BB"] + pB_MDH * tMatrix.temp[motherH, , "BW"]
    tMatrix[motherH, , "BH"] = tMatrix[motherH, , "BH"] + pH_MDH * tMatrix.temp[motherH, , "BW"]
    tMatrix[motherH, , "BR"] = tMatrix[motherH, , "BR"] + pR_MDH * tMatrix.temp[motherH, , "BW"]
    tMatrix[motherH, , "BW"] = tMatrix[motherH, , "BW"] - (1 - pW_MDH) * tMatrix.temp[motherH, , "BW"]
    tMatrix[motherH, , "BH"] = tMatrix[motherH, , "BH"] + pB_MDH * tMatrix.temp[motherH, , "HW"]
    tMatrix[motherH, , "HH"] = tMatrix[motherH, , "HH"] + pH_MDH * tMatrix.temp[motherH, , "HW"]
    tMatrix[motherH, , "HR"] = tMatrix[motherH, , "HR"] + pR_MDH * tMatrix.temp[motherH, , "HW"]
    tMatrix[motherH, , "HW"] = tMatrix[motherH, , "HW"] - (1 - pW_MDH) * tMatrix.temp[motherH, , "HW"]
    tMatrix[motherH, , "BR"] = tMatrix[motherH, , "BR"] + pB_MDH * tMatrix.temp[motherH, , "RW"]
    tMatrix[motherH, , "HR"] = tMatrix[motherH, , "HR"] + pH_MDH * tMatrix.temp[motherH, , "RW"]
    tMatrix[motherH, , "RR"] = tMatrix[motherH, , "RR"] + pR_MDH * tMatrix.temp[motherH, , "RW"]
    tMatrix[motherH, , "RW"] = tMatrix[motherH, , "RW"] - (1 - pW_MDH) * tMatrix.temp[motherH, , "RW"]
    
    # zygotes with two W alleles
    p_MDH_WW_BB = (factorial(2) / (factorial(2) * factorial(0) * factorial(0) *
                                     factorial(0))) * (pB_MDH ^ 2 * pH_MDH ^ 0 * pR_MDH ^ 0 * pW_MDH ^ 0)
    p_MDH_WW_HH = (factorial(2) / (factorial(0) * factorial(2) * factorial(0) *
                                     factorial(0))) * (pB_MDH ^ 0 * pH_MDH ^ 2 * pR_MDH ^ 0 * pW_MDH ^ 0)
    p_MDH_WW_RR = (factorial(2) / (factorial(0) * factorial(0) * factorial(2) *
                                     factorial(0))) * (pB_MDH ^ 0 * pH_MDH ^ 0 * pR_MDH ^ 2 * pW_MDH ^ 0)
    p_MDH_WW_WW = (factorial(2) / (factorial(0) * factorial(0) * factorial(0) *
                                     factorial(2))) * (pB_MDH ^ 0 * pH_MDH ^ 0 * pR_MDH ^ 0 * pW_MDH ^ 2)
    p_MDH_WW_BH = (factorial(2) / (factorial(1) * factorial(1) * factorial(0) *
                                     factorial(0))) * (pB_MDH ^ 1 * pH_MDH ^ 1 * pR_MDH ^ 0 * pW_MDH ^ 0)
    p_MDH_WW_BR = (factorial(2) / (factorial(1) * factorial(0) * factorial(1) *
                                     factorial(0))) * (pB_MDH ^ 1 * pH_MDH ^ 0 * pR_MDH ^ 1 * pW_MDH ^ 0)
    p_MDH_WW_BW = (factorial(2) / (factorial(1) * factorial(0) * factorial(0) *
                                     factorial(1))) * (pB_MDH ^ 1 * pH_MDH ^ 0 * pR_MDH ^ 0 * pW_MDH ^ 1)
    p_MDH_WW_HR = (factorial(2) / (factorial(0) * factorial(1) * factorial(1) *
                                     factorial(0))) * (pB_MDH ^ 0 * pH_MDH ^ 1 * pR_MDH ^ 1 * pW_MDH ^ 0)
    p_MDH_WW_HW = (factorial(2) / (factorial(0) * factorial(1) * factorial(0) *
                                     factorial(1))) * (pB_MDH ^ 0 * pH_MDH ^ 1 * pR_MDH ^ 0 * pW_MDH ^ 1)
    p_MDH_WW_RW = (factorial(2) / (factorial(0) * factorial(0) * factorial(1) *
                                     factorial(1))) * (pB_MDH ^ 0 * pH_MDH ^ 0 * pR_MDH ^ 1 * pW_MDH ^ 1)
    
    tMatrix[motherH, , "BB"] = tMatrix[motherH, , "BB"] + p_MDH_WW_BB * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "HH"] = tMatrix[motherH, , "HH"] + p_MDH_WW_HH * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "RR"] = tMatrix[motherH, , "RR"] + p_MDH_WW_RR * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "WW"] = tMatrix[motherH, , "WW"] - (1 - p_MDH_WW_WW) * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "BH"] = tMatrix[motherH, , "BH"] + p_MDH_WW_BH * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "BR"] = tMatrix[motherH, , "BR"] + p_MDH_WW_BR * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "BW"] = tMatrix[motherH, , "BW"] + p_MDH_WW_BW * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "HR"] = tMatrix[motherH, , "HR"] + p_MDH_WW_HR * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "HW"] = tMatrix[motherH, , "HW"] + p_MDH_WW_HW * tMatrix.temp[motherH, , "WW"]
    tMatrix[motherH, , "RW"] = tMatrix[motherH, , "RW"] + p_MDH_WW_RW * tMatrix.temp[motherH, , "WW"]
    
    #############################################################################
    # Mask - fitness
    #############################################################################
    
    # initialize viability mask to account for male fertility
    viabilityMask = array(
      data = 1,
      dim = c(size, size, size),
      dimnames = list(gtype, gtype, gtype)
    )
    
    for (i in 1:length(fc)) {
      viabilityMask[, , names(fc[i])] = fc[i] * viabilityMask[, , names(fc[i])]
    }
    
    # genotype-specific modifiers
    modifiers = cubeModifiers(
      gtype,
      eta   = eta,
      phi   = phi,
      omega = omega,
      xiF   = xiF,
      xiM   = xiM,
      s     = s
    )
    
    # put everything into a labeled list to return
    return(
      list(
        ih          = tMatrix,
        tau         = viabilityMask,
        genotypesID = gtype,
        genotypesN  = size,
        wildType    = "WW",
        eta         = modifiers$eta,
        phi         = modifiers$phi,
        omega       = modifiers$omega,
        xiF         = modifiers$xiF,
        xiM         = modifiers$xiM,
        s           = modifiers$s,
        releaseType = "HH"
      )
    )
  }
