
#' GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#' In case IC values are used as input, the added-evidence approach is used in which the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#'
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param Weights A matrix with weigths (i.e., IC weights or posterior model probabilities) of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained). Notably, only when the set of hypotheses cover the whol space / all theories (e.g., positive versus negative effect), then you can do without a safeguard hypothesis.
#' @param PriorWeights Optional. Vector containing 'NrHypos+1' numbers that represent the prior belief for this model. By default, equal prior weights are used (i.e., 1/(NrHypos+1)). Notably, in case the prior weights do not sum to 1, it will be rescaled such that it does; which implies that relative importance can be used and not per se weights.
#' @param Name_studies Optional. Vector of S numbers or S characters to be printed at the x-axis of the plot with GORIC(A) weights. Default: Name_studies = 1:S.
#' @param Name_Hypo Optional. Vector containing 'NrHypos+1' characters which will be used for labelling the hypothesis. Default: H1, H2, ....
#' @param PrintPlot Optional. Indicator whether plot of GORIC(A) weigths should be printed (TRUE; default) or not (FALSE). The GORIC(A) weights per study are plotted and the cumulative GORIC(A) weights (where those for the last study are the final ones).
#'
#' @return The output comprises, among other things, the cumulative and final evidence for the theory-based hypotheses.
#' @export
#' @examples
#'
#' S <- 4
#' Weights <- myWeights # Example based on S = 4 studies and 3 hypotheses:
#' # H0 <- "beta1 == 0"  # this hypothesis could have been left out
#' # Hpos <- "beta1 > 0"
#' # Hneg <- "beta1 < 0"
#' # Note that in this set the whole space is (all theories are) covered so the unconstrained is not needed as safeguard-hypothesis
#' GoricEvSyn_IC(S, Weights)
#'
#' # Change labels on x-axis in GORIC(A) weigths plot and give names to hypotheses #
#' # For example, let us say that the studies come from the years 2015, 2016, 2017, 2019.
#' # Because of unequal spacing, you may want to use numbers instead of characters:
#' Name_studies <- c(2015, 2016, 2017, 2019)
#' Name_Hypo <- c("H0", "Hpos", "Hneg")
#' GoricEvSyn_IC(S, Weights, Name_studies, Name_Hypo)

GoricEvSyn_weights <- function(S, Weights, PriorWeights = NULL, Name_studies = 1:S, Name_Hypo = NULL, PrintPlot = T) {

  # Checks on input
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
  if(length(dim(Weights)) != 2){
    print(paste0("The weights matrix Weights should have 2 dimensions; namely, S rows and 'NrHypos+1' columns. It should not be an array with more than 2 dimensions."))
    stop()
  }
  if(dim(Weights)[1] != S){
    print(paste0("The number of rows in the weights matrix Weights (", dim(Weights)[1], ") does not equal S = ", S, "."))
    stop()
  }
  NrHypos <- dim(Weights)[2] - 1
  #
  if(is.null(PriorWeights)){
    PriorWeights <- rep(1/(NrHypos + 1), (NrHypos + 1))
  }
  if(length(PriorWeights) != (NrHypos+1)){
    print(paste("The argument 'PriorWeights' should consist of 'NrHypos+1' = ", (NrHypos+1), " elements."))
    stop()
  }
  if(!all(is.numeric(PriorWeights))){
    print(paste("The argument 'PriorWeights' should consist of solely numbers ('NrHypos+1 = ", (NrHypos+1), " numbers)."))
    stop()
  }
  PriorWeights <- PriorWeights/sum(PriorWeights) # To make it sum to 1 (if it not already did)
  #
  if(PrintPlot != T & PrintPlot != F){
    print(paste("The argument 'PrintPlot' should be TRUE or FALSE, not ", PrintPlot, "."))
    stop()
  }
  if(length(Name_studies) != S){
    print(paste("The argument 'Name_studies' should consist of S = ", S, " elements (either all numbers or all characters)."))
    stop()
  }
  if(!all(is.numeric(Name_studies)) & !all(is.character(Name_studies))){
    print(paste("The argument 'Name_studies' should consist of either S = ", S, " numbers or S = ", S, " characters."))
    stop()
  }
  #
  if(is.null(Name_Hypo)){
    Name_Hypo <- paste0("H", 1:(NrHypos+1))
  }
  if(length(Name_Hypo) != (NrHypos+1)){
    print(paste("The argument 'Name_Hypo' should consist of 'NrHypos+1' = ", (NrHypos+1), " elements (all characters)."))
    stop()
  }
  if(!all(is.character(Name_Hypo))){
    print(paste("The argument 'Name_Hypo' should consist of solely characters ('NrHypos+1 = ", (NrHypos+1), " characters)."))
    stop()
  }


  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeWeights) <- Name_Hypo
  rownames(CumulativeWeights) <- c(paste0("Study", 1:S), "Final")

  CumulativeWeights[1,] <- PriorWeights * Weights[1,] / sum( PriorWeights * Weights[1,] )
  for(s in 2:S){
    CumulativeWeights[s,] <- CumulativeWeights[(s-1),] * Weights[s,] / sum( CumulativeWeights[(s-1),] * Weights[s,] )
  }
  EvSyn_approach <- "Added-evidence approach (which is the only option when the input consists of Weights)"

  CumulativeWeights[(S+1),] <- CumulativeWeights[S,]
  #
  Final.weights <- CumulativeWeights[S,]
  Final.rel.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.rel.weights) <- Name_Hypo
  colnames(Final.rel.weights) <- paste0("vs ", Name_Hypo)


  # Plot
  if(PrintPlot == T){
    weight_m <- Weights
    CumulativeGoricaWeights <- CumulativeWeights
    #
    NrHypos_incl <- (NrHypos + 1)
    Legend <- c("per study", "cumulative", Name_Hypo)
    Pch <- c(16, 8, rep(NA, NrHypos_incl))
    Col <- c(1, 1, 1:NrHypos_incl)
    Lty <- c(NA, 1, rep(1,NrHypos_incl))
    dev.off() # to reset the graphics pars to defaults
    par(mar=c(par('mar')[1:3], 0)) # optional, removes extraneous right inner margin space
    plot.new()
    l <- legend(0, 0, bty='n', Legend,
                plot=FALSE, pch=Pch, lty=Lty, col=Col)
    # calculate right margin width in ndc
    w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
    par(omd=c(0, 1-w, 0, 1))
    #
    teller_col <- 1
    #plot(1:S, weight_m[,1], pch = 16, col = teller_col, xlab = "Studies", ylab = "GORIC(A) weights", ylim = c(0,1), main = "GORIC(A) weights \n per study and cumulative")
    if(all(is.numeric(Name_studies))){
      X <- Name_studies
      plot(X, weight_m[,1], pch = 16, col = teller_col, xlab = "Studies", ylab = "GORIC(A) weights", ylim = c(0,1), main = "GORIC(A) weights \n per study and cumulative", xaxt="n")
      axis(1, at=X, labels=Name_studies)
    }else{
      X <- 1:S
      plot(X, weight_m[,1], pch = 16, col = teller_col, xlab = "Studies", ylab = "GORIC(A) weights", ylim = c(0,1), main = "GORIC(A) weights \n per study and cumulative", xaxt="n")
      axis(1, at=X, labels=Name_studies)
    }
    for(i in 2:NrHypos_incl){
      teller_col <- teller_col + 1
      points(X, weight_m[,i], pch = 16, col = teller_col)
    }
    teller_col <- 0
    for(i in 1:NrHypos_incl){
      teller_col <- teller_col + 1
      points(X, CumulativeGoricaWeights[1:S,i], pch = 8, col = teller_col)
      lines(X, CumulativeGoricaWeights[1:S,i], lty = 1, lwd = 1, col = teller_col)
    }
    #
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
           Legend, pch=Pch, lty=Lty, col=Col)
  }


  # Ouput
  final <- list(weight_m = Weights,
                EvSyn_approach = EvSyn_approach,
                Cumulative.weights = CumulativeWeights,
                Final.rel.weights = Final.rel.weights)
  return(final)

}
