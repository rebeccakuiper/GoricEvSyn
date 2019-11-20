
#' GORIC(A) evidence synthesis based on log likelihood and penalty values
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param TypeEv The type of evidence-synthesis approach: Equal-evidence approach (0) or Added-evidence approach (1).
#' In case of an equal-evidence approach, aggregating evidence from, say, 5 studies with n=100 observations is the same as obtaining evidence from 1 study (as if it was possible) with n=500 observations (like meta-analysis does).
#' In the added-evidence approach, the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param LL A matrix with log likelihood values of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained).
#' @param PT A matrix with penalty values of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained).
#'
#' @return The output comprises, among others, the overall evidence for the theory-based hypotheses.
#' @export
#' @examples
#'
#' S <- 4
#' # Example based on S = 4 studies and 3 hypotheses:
#' # H0 <- "beta1 == 0"
#' # Hpos <- "beta1 > 0"
#' # Hneg <- "beta1 < 0"
#' # Note that in this set the whole space is (all theories are) covered so the unconstrained is not needed as safeguard-hypothesis
#' LL <- myLLs
#' PT <- myPTs
#'
#' # Added-evidence approach
#' TypeEv <- 1
#' GoricEvSyn_LLandPT(TypeEv, S, LL, PT)
#'
#' # Equal-evidence approach
#' TypeEv <- 0
#' GoricEvSyn_LLandPT(TypeEv, S, LL, PT)


GoricEvSyn_LLandPT <- function(TypeEv, S, LL, PT) {

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
  if(length(dim(LL)) != 2){
    print(paste0("The LL matrix should have 2 dimensions; namely, S rows and 'NrHypos+1' columns. It should not be an array with more than 2 dimensions."))
    stop()
  }
  if(dim(LL)[1] != S){
    print(paste0("The number of rows in the LL matrix (", dim(LL)[1], ") does not equal S = ", S, "."))
    stop()
  }
  if(length(dim(PT)) != 2){
    print(paste0("The PT matrix should have 2 dimensions; namely, S rows and 'NrHypos+1' columns. It should not be an array with more than 2 dimensions."))
    stop()
  }
  if(dim(PT)[1] != S){
    print(paste0("The number of rows in the PT matrix (", dim(PT)[1], ") does not equal S = ", S, "."))
    stop()
  }
  NrHypos <- dim(LL)[2] - 1
  NrHypos_PT <- dim(PT)[2] - 1
  if(NrHypos != NrHypos_PT){
    print(paste0("The number of columns in the LL matrix (", dim(LL)[2], ") does not equal the number of columns in the PT matrix (", dim(PT)[2], "). Both should equal 'NrHypos+1'."))
    stop()
  }


  overallGorica <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  overallGoricaWeights <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  colnames(overallGorica) <- colnames(overallGoricaWeights) <- colnames(LL) <- colnames(PT) <- paste0("H", 1:(NrHypos + 1))
  rownames(overallGorica) <- rownames(overallGoricaWeights) <- rownames(LL) <- rownames(PT) <- paste0("Study", 1:S)


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

  final <- list(LL = LL, PT = PT,
                  EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                  Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  return(final)

}
