
#' GORIC(A) evidence synthesis
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param TypeEv The type of evidence-synthesis approach: Equal-evidence approach (0) or Added-evidence approach (1).
#' In case of an equal-evidence approach, aggregating evidence from, say, 5 studies with n=100 observations is the same as obtaining evidence from 1 study (as if it was possible) with n=500 observations (like meta-analysis does).
#' In the added-evidence approach, the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param Param_studies List of S 'named' vectors with the (standardized) parameter estimates of interest. Thus, there are S items in the list and each item is a 'named' vector with the number of parameter estimates relevant for that study. Note: The names of the estimates should be used in the hypothesis specification.
#' @param CovMx_studies List of the S covariance matrices of the (standardized) parameter estimates of interest.
#' @param SameHypo Indicator whether the same hypotheses are used (1) or not (0). If SameHypo = 1, then the same estimates in Param_studies should have the same name.
#' @param NrHypos The number of theory-based hypotheses that will be evaluated.
#' @param Hypo_studies The NrHypos theory-based hypotheses. If SameHypo = 0, then there should be S specifications of the NrHypos theory-based hypotheses.
#' @param ComplOrUnc Indicator of which safeguard-hypothesis should be used in case 'NrHypos = 1': Complement (0; i.e., the remaining theories) or unconstrained (1; i.e., all possible theories including the one specified). If NrHypos > 1, then automatically the unconstrained will be used.
#'
#' @return The output comprises, among others, the overall evidence for the theory-based hypotheses.
#' @importFrom restriktor goric
#' @export
#' @examples
#'
#' # In progress

# TO DO denk na over hoe parameters e.d. opgeven!

GoricEvSyn <- function(TypeEv, S, Param_studies, CovMx_studies, SameHypo, NrHypos, Hypo_studies, ComplOrUnc = 1) {

  # TO DO
  # Checks op input!

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


  GORICA_m  <- array(data=NA, dim=c(S,(NrHypos+1)))
  weight_m <- array(data=NA, dim=c(S,(NrHypos+1)))
  LL <- matrix(NA, nrow = S, ncol = (NrHypos+1))
  PT <- matrix(NA, nrow = S, ncol = (NrHypos+1))
  rownames(LL) <- rownames(PT) <- paste0("Study", 1:S)
  overallGorica <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  overallGoricaWeights <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  rownames(overallGorica) <- rownames(overallGoricaWeights) <- paste0("Study", 1:S)
  if(NrHypos == 1 & ComplOrUnc == 0){
    colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(overallGorica) <- colnames(overallGoricaWeights) <- c("H1", "Hc1")
    rel.weight_mu <- array(data=NA, dim=c(S,1))
  }else{
    colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(overallGorica) <- colnames(overallGoricaWeights) <- c(paste0("H", 1:NrHypos), "Hu")
    rel.weight_mu <- array(data=NA, dim=c(S,(NrHypos+1)))
  }
  rownames(GORICA_m) <- rownames(rel.weight_mu) <- rownames(weight_m) <- paste0("Study", 1:S)
  #
  #
  param <- Param_studies[c(F,T)] # TO DO
  names <- Param_studies[c(T,F)]
  #outputList <- list()
  if(SameHypo==1){ # if same hypotheses for all studies
    for(HypoTeller in 1:NrHypos){
      assign(paste0("H", HypoTeller), Hypo_studies[HypoTeller])
    }
    HypoSet <- noquote(paste0("H", 1:NrHypos, collapse = ", "))
    #HypoSet <- paste(Hypo_studies[1:NrHypos], collapse = "; ")
    for(s in 1:S){
      # GORICA
      est <- as.vector(as.matrix(read.table(text=param[s])))
      names(est) <- as.vector(as.matrix(read.table(text=names[s])))
      njoint <- length(est)
      cov <- CovMx_studies[1:njoint,1:njoint,s]
      #
      # Run GORICA
      if(NrHypos == 1 & ComplOrUnc == 0){ # vs complement
        eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                                HypoSet,
                                ", type = 'gorica', comparison = 'complement')")))
        rel.weight_mu[s,] <- result2$relative.gw[1, (NrHypos+1)]
      } else{ # vs unconstrained (default)
        eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                                HypoSet,
                                ", type = 'gorica')")))
        #result2 <- goric(est, VCOV = cov, HypoSet, type = "gorica")
        rel.weight_mu[s,] <- result2$relative.gw[, (NrHypos+1)]
      }
      LL[s,] <- result2$result[,2]
      PT[s,] <- result2$result[,3]
      GORICA_m[s,] <- result2$result[,4]
      weight_m[s,] <- result2$result[,5]
    }
  }else{ # if NOT same Hypo's per study
    teller <- 0
    for(s in 1:S){
      for(HypoTeller in 1:NrHypos){
        assign(paste0("H", HypoTeller), Hypo_studies[(teller+HypoTeller)])
      }
      HypoSet <- noquote(paste0("H", 1:NrHypos, collapse = ", "))
      #HypoSet <- paste(Hypo_studies[(teller+1):(teller+NrHypos)], collapse = "; ")
      #
      # GORICA
      est <- as.vector(as.matrix(read.table(text=param[s])))
      names(est) <- as.vector(as.matrix(read.table(text=names[s])))
      njoint <- length(est)
      cov <- CovMx_studies[1:njoint,1:njoint,s]
      #
      # Run GORICA
      if(NrHypos == 1 & ComplOrUnc == 0){
        eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                                HypoSet,
                                ", type = 'gorica', comparison = 'complement')")))
        rel.weight_mu[s,] <- result2$relative.gw[1, (NrHypos+1)]
        colnames(rel.weight_mu) <- c("H1 vs Hc1")
      } else{
        eval(parse(text = paste("result2 <- goric(est, VCOV = cov, ",
                                HypoSet,
                                ", type = 'gorica')")))
        #result2 <- goric(est, VCOV = cov, HypoSet, type = "gorica")
        rel.weight_mu[s,] <- result2$relative.gw[, (NrHypos+1)]
      }
      LL[s,] <- result2$result[,2]
      PT[s,] <- result2$result[,3]
      GORICA_m[s,] <- result2$result[,4]
      weight_m[s,] <- result2$result[,5]
      #
      teller <- teller + NrHypos
    }
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
  rownames(Final.rel.GORICA.weights) <- c(paste0("H", 1:(NrHypos + 1)))
  colnames(Final.rel.GORICA.weights) <- c(paste0("vs H", 1:(NrHypos + 1)))

  if(NrHypos == 1 & ComplOrUnc == 0){
    colnames(rel.weight_mu) <- c("H1 vs Hc1")
    final <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, rel.GORICA.weight_mc = rel.weight_mu, LL_m = LL, PT_m = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  } else{
    colnames(rel.weight_mu) <- c(paste0("H", 1:NrHypos, " vs Unc."), "Unc. vs Unc.")
    final <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, rel.GORICA.weight_mu = rel.weight_mu, LL_m = LL, PT_m = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  }
  return(final)

}
