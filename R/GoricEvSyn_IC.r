
#' GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#' In case IC values are used as input, the added-evidence approach is used in which the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#'
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param IC A matrix with information criteria (AIC, ORIC, GORIC, or GORICA) values of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained).
#'
#' @return The output comprises, among others, the overall evidence for the theory-based hypotheses.
#' @export
#' @examples
#'
#' S <- 4
#' IC <- myGORICs # Example based on S = 4 studies and 3 hypotheses:
#' # H0 <- "beta1 == 0"
#' # Hpos <- "beta1 > 0"
#' # Hneg <- "beta1 < 0"
#' # Note that in this set the whole space is (all theories are) covered so the unconstrained is not needed as safeguard-hypothesis
#' GoricEvSyn_IC(S, IC)


GoricEvSyn_IC <- function(S, IC) {

  # Checks op input
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
  if(length(dim(IC)) != 2){
    print(paste0("The IC matrix should have 2 dimensions; namely, S rows and 'NrHypos+1' columns. It should not be an array with more than 2 dimensions."))
    stop()
  }
  if(dim(IC)[1] != S){
    print(paste0("The number of rows in the IC matrix (", dim(IC)[1], ") does not equal S = ", S, "."))
    stop()
  }
  NrHypos <- dim(IC)[2] - 1


  overallGorica <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  overallGoricaWeights <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  colnames(overallGorica) <- colnames(overallGoricaWeights) <- colnames(IC) <- paste0("H", 1:(NrHypos + 1))
  rownames(overallGorica) <- rownames(overallGoricaWeights) <- rownames(IC) <- paste0("Study", 1:S)


  sumIC <- 0
  for(s in 1:S){
    sumIC <- sumIC + IC[s,]
    overallGorica[s,] <- sumIC
    overallGoricaWeights[s,] <- exp(-0.5*overallGorica[s,]) / sum(exp(-0.5*overallGorica[s,]))
  }
  EvSyn_approach <- "Added-evidence approach"

  Final.GORICA <- matrix(overallGorica[S,], nrow = 1)
  Final.GORICA.weights <- overallGoricaWeights[S,]
  Final.rel.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  Final.GORICA.weights <- matrix(Final.GORICA.weights, nrow = 1)
  rownames(Final.GORICA.weights) <- "Final"
  rownames(Final.rel.GORICA.weights) <- c(paste0("H", 1:(NrHypos + 1)))
  colnames(Final.rel.GORICA.weights) <- c(paste0("vs H", 1:(NrHypos + 1)))

  final <- list(IC = IC,
                EvSyn_approach = EvSyn_approach, overallGorica = overallGorica, overallGoricaWeights = overallGoricaWeights,
                Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  return(final)

}
