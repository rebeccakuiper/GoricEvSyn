
#' GORIC(A) evidence synthesis
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param TypeEv The type of evidence-synthesis approach: Equal-evidence approach (0) or Added-evidence approach (1).
#' In case of an equal-evidence approach, aggregating evidence from, say, 5 studies with n=100 observations is the same as obtaining evidence from 1 study (as if it was possible) with n=500 observations (like meta-analysis does).
#' In the added-evidence approach, the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param Param_studies List of S 'named' vectors with the k_s (standardized) parameter estimates of interest of Study s. Thus, there are S items in the list and each item is a 'named' vector with k_s elements: the k_s number of parameter estimates relevant for that study. In case each study has the same number of parameters (k) which denote the same (in terms of hypothesis specification), Param_studies can be an S x k 'named' matrix. Note: The names of the vectors (or the column names of the S x 'k' matrix) with estimates should be used in the hypothesis specification.
#' @param CovMx_studies List of the S covariance matrices of the (standardized) parameter estimates of interest (of size k_s x k_s). In case number of parameters are the same, it can also be a S*k_s x k_s matrix. Note: The columns (and rows) do not need to be named.
#' @param SameHypo Indicator whether the same hypotheses are used (1) or not (0) in all S studies. If SameHypo = 1, then the same estimates in Param_studies should have the same name.
#' @param NrHypos The number of theory-based hypotheses that will be evaluated within each study (is a scalar with an integer value).
#' @param Hypo_studies A vector of strings containing the NrHypos theory-based hypotheses. If SameHypo = 0, then there should be S specifications of the NrHypos theory-based hypotheses, that is S times NrHypos strings.
#' @param Safeguard Indicator of which safeguard-hypothesis should be used: "unconstrained" (default; i.e., all possible theories including the one specified), "none" (only advised when set of hyptheses cover all theories), or (only when 'NrHypos = 1') "complement" (i.e., the remaining theories).
#'
#' @return The output comprises, among others, the overall evidence for the theory-based hypotheses.
#' @importFrom restriktor goric
#' @export
#' @examples
#'
#' ### Example 1: 2 ANOVA studies (Monin and Holubar) ###
#' ### We will use the estimates of the ANOVA models
#'
#' est_1 <- c(1.88   2.54   0.02)
#' names(est_1) <- c("group1", "group2", "group3")
#' vcov_est_1 <- diag(0.2149074, 0.2149074, 0.1408014)
#' est_2 <- c(0.98, 0.02, 0.27)
#' names(est_2) <- c("gr1", "gr2", "gr3") # Use different names to show use of different set of hypotheses
#' vcov_est_2 <- diag(0.1382856, 0.1024337, 0.0987754)
#'
#' # If number of parameters differ per study (but can also be used when they are the same): make lists
#' #
#' # beta values from the analyses
#' Param_studies <- list(est_1, est_2)
#' Param_studies
#' #
#' # standard error of the beta's (from the S primary studies)
#' CovMx_studies <- list(vcov_est_1, vcov_est_2)
#' CovMx_studies
#'
#' # Set of hypotheses for each study
#' # Note: in this case we could make the names of the estimates in est_1 and est_2 the same.
#' # In general, when not same number of parameters they names will not be the same and/or the set of hypotheses, so we use that here.
#' SameHypo <- 0
#' NrHypos <- 2
#' # names(est_1) # Specify restrictions using those names
#' H11 <- 'group1 == group2; group2 > group3'  # Note: cannot use group1 == group2 > group3
#' H12 <- 'group2 > group1; group1 > group3'   # Note: cannot use group2 > group1 > group3
#' # names(est_2) # Specify restrictions using those names
#' H21 <- 'gr1 == gr2; gr2 > gr3'  # Note: cannot use gr1 == gr2 > gr3
#' H22 <- 'gr2 > gr1; gr1 > gr3'   # Note: cannot use gr2 > gr1 > gr3
#' Hypo_studies <- c(H11, H12, H21, H22)
#' #
#' # Evidence synthesis
#' TypeEv <- 1 # Added-evidence approach
#' GoricEvSyn(TypeEv, S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies)
#'
#'
#'
#' ### Example 2: 4 trust studies discussed in Kuiper et al. (2013) ###
#' ### Different statistical models are used: linear regression, probit regression, and three-level logististic regression
#'
#' S <- 4 # Number of primary studies for which the evidence for a central theory should be determined and synthesized.
#'
#' # Example 1a #
#' # If same number of parameters per study
#' #
#' NrParam <- 1 # In this example, each study has one parameter of interest (i.e., hypotheses below address only this parameter).
#' # beta values from the analyses
#' Param_studies <- matrix(c(0.09, 0.14, 1.09, 1.781), nrow = S, ncol = NrParam)
#' colnames(Param_studies) <- "beta1" # Should have names and should have the same as in the hypotheses below
#' Param_studies
#' #
#' # standard error of the beta's (from the S primary studies)
#' CovMx_studies <- matrix(c(0.029^2, 0.054^2, 0.093^2, 0.179^2), nrow = S, ncol = NrParam) # Note: no names needed
#' #
#'
#' # OR
#'
#' # Example 1b #
#' # If number of parameters differ per study (but can also be used when they are the same)
#' #
#' # beta values from the analyses
#' est_1 <- matrix(c(0.09), nrow = 1)
#' colnames(est_1) <- "beta1"
#' est_2 <- matrix(c(0.14), nrow = 1)
#' colnames(est_2) <- "beta1"
#' est_3 <- matrix(c(1.09), nrow = 1)
#' colnames(est_3) <- "beta1"
#' est_4 <- matrix(c(1.781), nrow = 1)
#' colnames(est_4) <- "beta1"
#' Param_studies <- list(est_1, est_2, est_3, est_4)
#' Param_studies
#' #
#' # standard error of the beta's (from the S primary studies)
#' vcov_est_1 <- matrix(c(0.029^2), nrow = 1)
#' vcov_est_2 <- matrix(c(0.054^2), nrow = 1)
#' vcov_est_3 <- matrix(c(0.093^2), nrow = 1)
#' vcov_est_4 <- matrix(c(0.179^2), nrow = 1)
#' CovMx_studies <- list(vcov_est_1, vcov_est_2, vcov_est_3, vcov_est_4)
#' CovMx_studies
#'
#'
#' # Set of hypotheses for each study
#' # Note: in this case the same for each study
#' SameHypo <- 1
#' NrHypos <- 3
#' H0 <- "beta1 == 0"
#' Hpos <- "beta1 > 0"
#' Hneg <- "beta1 < 0"
#' Hypo_studies <- c(H0, Hpos, Hneg)
#' # Since this covers the whole space / covers all theories, we do not need a safeguard-hypothesis
#' Safeguard <- "none"
#' #
#' # Evidence synthesis
#' TypeEv <- 1 # Added-evidence approach
#' GoricEvSyn(TypeEv, S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies, Safeguard)
#'
#'
#'
#' ### Example 3: 3 'fictional' studies; for ease, all 3 use linear regression ###
#' ### where continuous predictor variable need to be standardized and
#' ### where we will the complement of our theory as competing hypothesis ###
#'
#' S <- 3
#' ratio <- c(1,1.1,1.2)
#'
#' # Generate data1
#' n1 <- 30
#' x11 <- rnorm(n1)
#' x12 <- rnorm(n1)
#' x13 <- rnorm(n1)
#' data <- cbind(x11, x12, x13)
#' # Standardize data - since parameters for continuous variables will be compared
#' data1 <- as.data.frame(scale(data))
#' y1 <- ratio[1]*data1$x11 + ratio[2]*data1$x12 + ratio[3]*data1$x13 + rnorm(n1)
#' # Note: since there is one outcome, the outcome does not need to be standardized
#' # Fit regression model
#' fit.lm1 <- lm(y1 ~ -1 + x11 + x12 + x13, data = data1)
#'
#' # Generate data2
#' n2 <- 50
#' x21 <- rnorm(n2)
#' x22 <- rnorm(n2)
#' x23 <- rnorm(n2)
#' data <- cbind(x21, x22, x23)
#' # Standardize data - since parameters for continuous variables will be compared
#' data2 <- as.data.frame(scale(data))
#' y2 <- ratio[1]*data2$x21 + ratio[2]*data2$x22 + ratio[3]*data2$x23 + rnorm(n2)
#' # Note: since there is one outcome, the outcome does not need to be standardized
#' # Fit regression model
#' fit.lm2 <- lm(y2 ~ -1 + x21 + x22 + x23, data = data2)
#'
#' # Generate data3
#' n3 <- 100
#' x31 <- rnorm(n3)
#' x32 <- rnorm(n3)
#' x33 <- rnorm(n3)
#' data <- cbind(x31, x32, x33)
#' # Standardize data - since parameters for continuous variables will be compared
#' data3 <- as.data.frame(scale(data))
#' y3 <- ratio[1]*data3$x31 + ratio[2]*data3$x32 + ratio[3]*data3$x33 + rnorm(n3)
#' # Note: since there is one outcome, the outcome does not need to be standardized
#' # Fit regression model
#' fit.lm3 <- lm(y3 ~ -1 + x31 + x32 + x33, data = data3)
#'
#' # Extract estimates and their covariance matrix (per study)
#' est_1 <- coef(fit.lm1)
#' est_2 <- coef(fit.lm2)
#' est_3 <- coef(fit.lm3)
#' vcov_est_1 <- vcov(fit.lm1)
#' vcov_est_2 <- vcov(fit.lm2)
#' vcov_est_3 <- vcov(fit.lm3)
#' #
#' #
#' ## If same number and type of parameters per study (since they will obtain the same name)
#' #
#' NrParam <- length(est_1)
#' # Parameter estimate values from the S primary studies
#' Param_studies <- matrix(c(est_1, est_2, est_3), byrow = T, nrow = S, ncol = NrParam)
#' colnames(Param_studies) <- c("x1", "x2", "x3")
#' Param_studies
#' # standard error of the beta's (from the S primary studies)
#' CovMx_studies <- matrix(c(vcov_est_1, vcov_est_2, vcov_est_3), byrow = T, nrow = S*NrParam, ncol = NrParam) # Note: no names needed
#' CovMx_studies
#' #
#' # Set of hypotheses for each study
#' # Note: in this case the same for each study
#' SameHypo <- 1
#' # colnames(est) # Specify restrictions using those names
#' H1 <- 'x1 < x2; x2 < x3'   # Note: cannot use x1 < x2 < x3
#' # Since estimates of continuous variables are compared in our theory, we standardized the data before to obtain comparable estimates.
#' NrHypos <- 1
#' # Since we have only one theory-based hypothesis, we will use the (more powerful) complement of the hypothesis (Vanbrabant, Van Loey, Kuiper, 2019)
#' # The complement represents the remaing 11 theories, while the unconstrained reflects all 12 possible theories including H1.
#' Safeguard <- "complement"
#' # Evidence synthesis
#' TypeEv <- 1 # Added-evidence approach
#' GoricEvSyn(TypeEv, S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies = H1, Safeguard)



GoricEvSyn <- function(TypeEv, S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies, Safeguard = "unconstrained") {

  # Checks op input
  #
  if(length(TypeEv) != 1){
    print(paste("The type of evidence-synthesis approach (TypeEv) should be a scalar; more specifically, it should be 0 or 1."))
    stop()
  }
  if(TypeEv != 1 & TypeEv != 0){
    print(paste("The type of evidence-synthesis approach (TypeEv) should be 0 or 1."))
    stop()
  }
  #
  if(length(S) != 1){
    print(paste("The number of studies (S) should be a scalar; more specifically, an integer value."))
    stop()
  }
  if(S %% 1 != 0){
    print(paste("The number of studies (S) should be an integer value."))
    stop()
  }
  #
  Param_list <- NA
  if(!is.list(Param_studies)){ # Not a list
    if(length(dim(Param_studies)) != 2){ # Check if matrix
      print(paste("The argument Param_studies should be a matrix (i.e., 2-dimensional array) of size S x 'number of parameters in each study (k)' or it should be a list of S matrices."))
      stop()
    }else if(dim(Param_studies)[1] != S){ # Matrix, but not S rows
      print(paste("The number of rows in the argument Param_studies does not equal S = ", S, "."))
      stop()
    }else if(is.null(colnames(Param_studies))){ # Matrix with S rows, but no column names
      print(paste("The S x 'k' matrix Param_studies does not have column names (which should be the same as in the hypotheses)."))
      stop()
    }else{ # If S x k matrix
      Param_list <- 0
    }
  }else{ # If list
    if(length(Param_studies) != S){ # does not contain S items
      print(paste("The list Param_studies does not have S = ", S, "items."))
      stop()
    }
    for(s in 1:S){ # Check if each item...
      est <- Param_studies[[s]]
      if(!is.vector(est)){# Check if each item = vector
        if(length(dim(est)) == 2 & dim(est)[1] == 1){ # Check if matrix with one row, then make a vector
          if(is.null(colnames(est))){ # Check if matrix has names
            print(paste0("The matrix in the ", s, "th element of the list Param_studies does not have column names (which should be the same as in the hypotheses)."))
            stop()
          }else{ # make a named vector
            Param_studies[[s]] <- as.vector(est)
            names(Param_studies[[s]]) <- colnames(est)
          }
        } else{ # not a vector and not a 1 x k_smatrix
          print(paste0("The ", s, "th element of the list Param_studies does not contain a (named) vector or (named) matrix with one row."))
          stop()
        }
      }else if(is.null(names(est))){ # Check if each item = vector has names
        print(paste0("The vector in the ", s, "th element of the list Param_studies does not have names (which should be the same as in the hypotheses).."))
        stop()
      }
    }
    Param_list <- 1
  }
  #
  CovMx_list <- NA
  if(!is.list(CovMx_studies)){ # Not a list
    if(length(dim(CovMx_studies)) != 2){ # Check if matrix
      print(paste("The argument CovMx_studies should be a matrix (i.e., 2-dimensional array) of size S*k x k, with k = number of parameters in each study, or it should be a list of S (k_s x k_s) matrices."))
      stop()
    }else if(dim(CovMx_studies)[1] != S*dim(CovMx_studies)[2]){ # Matrix, but not S*k rows
      print(paste0("The number of rows in the argument CovMx_studies does not equal S * 'k' = ", S, " * ", dim(CovMx_studies)[2], " = ", S*dim(CovMx_studies)[2], "."))
      stop()
    }else{ # If S*k x k matrix
      k <- dim(CovMx_studies)[2]
      if(k > 1){
        teller_s <- 0
        for(s in 1:S){ # Check if each item...
          CovMx <- CovMx_studies[(teller_s+1):(teller_s+k),]
          if(!isSymmetric(CovMx)){ # Matrix, but not k_s rows and columns
            print(paste0("The ", s, "th covariance matrix in CovMx_studies is not a symmetric matrix."))
            stop()
          }
          teller_s <- teller_s + k
        }
      }
      CovMx_list <- 0
    }
  }else{ # If list
    if(length(CovMx_studies) != S){ # does not contain S items
      print(paste("The list CovMx_studies does not have S =", S, "items."))
      stop()
    }
    for(s in 1:S){ # Check if each item...
      CovMx <- CovMx_studies[[s]]
      if(length(dim(CovMx)) != 2){ # Check if matrix
        print(paste0("The ", s, "th element of the list CovMx_studies does not contain a matrix (i.e., 2-dimensional array of size k_s x k_s, with k_s = number of parameters in Study s)."))
        stop()
      }else if(dim(CovMx)[1] != dim(CovMx)[2]){ # Matrix, but not k_s rows and columns
        print(paste0("The ", s, "th element of the list CovMx_studies does not contain a square matrix of size k_s x k_s, with k_s = number of parameters in Study s."))
        stop()
      }else if(dim(CovMx)[1] > 1){
        if(!isSymmetric(CovMx)){ # Matrix, but not k_s rows and columns
          print(paste0("The ", s, "th element of the list CovMx_studies does not contain a symmetric matrix."))
          stop()
        }
      }
    }
    CovMx_list <- 1
  }
  #
  if(length(SameHypo) != 1){
    print(paste("The indicator for using the same set of hypotheses in all studies (SameHypo) should be a scalar; more specifically, it should be 0 or 1."))
    stop()
  }
  if(SameHypo != 1 & SameHypo != 0){
    print(paste("The indicator for using the same set of hypotheses in all studies (SameHypo) should be 0 or 1."))
    stop()
  }
  #
  if(length(NrHypos) != 1){
    print(paste("The number of hypotheses used in all studies (NrHypos) should be a scalar (integer)."))
    stop()
  }
  if(length(NrHypos) == 1){
    if(NrHypos %% 1 != 0){
      print(paste("The number of hypotheses used in all studies (NrHypos) should be an integer value."))
      stop()
    }
  }
  #
  if(SameHypo ==1){
    if(length(Hypo_studies) != NrHypos){
      print(paste("The argument (Hypo_studies) should consist of NrHypos  = ", NrHypos, " elements (which are character strings / text elements). Note: SameHypo = 1."))
      stop()
    }
  }else if(length(Hypo_studies) != S*NrHypos){
    print(paste("The argument (Hypo_studies) should consist of S * NrHypos  = ", S, " * ", NrHypos, " = ", S*NrHypos, " elements (which are character strings / text elements). Note: SameHypo = 0."))
    stop()
  }
  if(!is.character(Hypo_studies)){
    print(paste("The argument (Hypo_studies) does not contain character strings / text elements."))
    stop()
  }
  #
  if(NrHypos == 1){
    if(length(Safeguard) != 1){
      print(paste("The type of safeguard-hypothesis (Safeguard) should be one word ('unconstrained', 'none', or (if NrHypos = 1) 'complement')."))
      stop()
    }
    if(Safeguard != "unconstrained" & Safeguard != "none" & Safeguard != "complement"){
      print(paste("The type of safeguard-hypothesis (Safeguard) should be 'unconstrained', 'none', or (if NrHypos = 1) 'complement'."))
      stop()
    }
  }

  NrHypos_incl <- NrHypos + 1
  if(Safeguard == "none"){
    NrHypos_incl <- NrHypos
  }
  GORICA_m  <- array(data=NA, dim=c(S,NrHypos_incl))
  weight_m <- array(data=NA, dim=c(S,NrHypos_incl))
  LL <- matrix(NA, nrow = S, ncol = NrHypos_incl)
  PT <- matrix(NA, nrow = S, ncol = NrHypos_incl)
  rownames(LL) <- rownames(PT) <- paste0("Study", 1:S)
  overallGorica <- matrix(NA, nrow = S, ncol = NrHypos_incl)
  overallGoricaWeights <- matrix(NA, nrow = S, ncol = NrHypos_incl)
  rownames(overallGorica) <- rownames(overallGoricaWeights) <- paste0("Study", 1:S)
  if(NrHypos == 1 & Safeguard == "complement"){
    colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(overallGorica) <- colnames(overallGoricaWeights) <- c("H1", "Hc1")
    rel.weight_mu <- array(data=NA, dim=c(S,1))
  }else if(Safeguard == "none"){
    colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(overallGorica) <- colnames(overallGoricaWeights) <- c(paste0("H", 1:NrHypos))
    rel.weight_mu <- array(data=NA, dim=c(S,NrHypos_incl))
  }else{
    colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(overallGorica) <- colnames(overallGoricaWeights) <- c(paste0("H", 1:NrHypos), "Hu")
    rel.weight_mu <- array(data=NA, dim=c(S,NrHypos_incl))
  }
  rownames(GORICA_m) <- rownames(rel.weight_mu) <- rownames(weight_m) <- paste0("Study", 1:S)
  #
  #
  if(SameHypo==1){ # if same hypotheses for all studies
    for(HypoTeller in 1:NrHypos){
      assign(paste0("H", HypoTeller), Hypo_studies[HypoTeller])
    }
  }
  teller <- 0
  teller_k <- 0
  for(s in 1:S){
    #
    if(SameHypo==0){ # if NOT same hypotheses for all studies
      for(HypoTeller in 1:NrHypos){
        assign(paste0("H", HypoTeller), Hypo_studies[(teller+HypoTeller)])
      }
    }
    #
    HypoSet <- noquote(paste0("H", 1:NrHypos, collapse = ", "))
    #
    # GORICA
    if(Param_list == 0){
      est <- Param_studies[s,]
      names(est) <- colnames(Param_studies)
    }else{
      est <- Param_studies[[s]]
    }
    if(CovMx_list == 0){
      cov <- CovMx_studies[(teller_k+1):(teller_k+k),]
    }else{
      cov <- CovMx_studies[[s]]
    }
    if(length(cov) == 1){
      cov <- matrix(cov)
    }
    #
    # Run GORICA
    if(NrHypos == 1 & Safeguard == "complement"){ # vs complement
      eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                              HypoSet,
                              ", type = 'gorica', comparison = Safeguard)")))
      rel.weight_mu[s,] <- result2$relative.gw[1, NrHypos_incl]
    } else{ # vs unconstrained (default)
      eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                              HypoSet,
                              ", type = 'gorica', comparison = Safeguard)")))
      #result2 <- goric(est, VCOV = cov, HypoSet, type = "gorica", comparison = Safeguard)
      if(Safeguard == "unconstrained"){
        rel.weight_mu[s,] <- result2$relative.gw[, NrHypos_incl]
      }
    }
    LL[s,] <- result2$result[,2]
    PT[s,] <- result2$result[,3]
    GORICA_m[s,] <- result2$result[,4]
    weight_m[s,] <- result2$result[,5]
    #
    teller <- teller + NrHypos
    teller_k <- teller_k + k
  }

  sumLL <- 0
  sumPT <- 0
  if(TypeEv == 1){ # added-ev approach
    for(s in 1:S){
      sumLL <- sumLL + LL[s,]
      sumPT <- sumPT + PT[s,]
      overallGorica[s,] <- -2 * sumLL + 2 * sumPT
      overallGoricaWeights[s,] <- exp(-0.5*overallGorica[s,]) / sum(exp(-0.5*overallGorica[s,]))
    }
    EvSyn_approach <- "Added-evidence approach"
  }else{ # equal-ev approach
    for(s in 1:S){
      sumLL <- sumLL + LL[s,]
      sumPT <- sumPT + PT[s,]
      overallGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      overallGoricaWeights[s,] <- exp(-0.5*overallGorica[s,]) / sum(exp(-0.5*overallGorica[s,]))
    }
    EvSyn_approach <- "Equal-evidence approach"
  }

  Final.GORICA <- matrix(overallGorica[S,], nrow = 1)
  Final.GORICA.weights <- overallGoricaWeights[S,]
  Final.rel.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  Final.GORICA.weights <- matrix(Final.GORICA.weights, nrow = 1)
  rownames(Final.GORICA.weights) <- "Final"
  #
  colnames(Final.GORICA) <- colnames(Final.GORICA.weights) <- c(paste0("H", 1:NrHypos_incl))
  rownames(Final.rel.GORICA.weights) <- c(paste0("H", 1:NrHypos_incl))
  colnames(Final.rel.GORICA.weights) <- c(paste0("vs H", 1:NrHypos_incl))

  if(NrHypos == 1 & Safeguard == "complement"){
    colnames(rel.weight_mu) <- c("H1 vs Hc1")
    colnames(Final.GORICA) <- colnames(Final.GORICA.weights) <- c("H1", "Hc1")
    rownames(Final.rel.GORICA.weights) <- c("H1", "Hc1")
    colnames(Final.rel.GORICA.weights) <- c("vs H1", "vs Hc1")
    final <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, rel.GORICA.weight_mc = rel.weight_mu, LL_m = LL, PT_m = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  } else if(Safeguard == "none"){
    final <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, LL_m = LL, PT_m = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  }else{ # unc
    colnames(rel.weight_mu) <- c(paste0("H", 1:NrHypos, " vs Unc."), "Unc. vs Unc.")
    final <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, rel.GORICA.weight_mu = rel.weight_mu, LL_m = LL, PT_m = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  }
  return(final)

}
