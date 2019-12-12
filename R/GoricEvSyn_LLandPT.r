
#' GORIC(A) evidence synthesis based on log likelihood and penalty values
#'
#' GORIC(A) evidence synthesis (GoricEvSyn) aggregates the evidence for theory-based hypotheses from multiple studies that may use diverse designs to investigate the same central theory. There is also an interactive web application on my website to perform GoricEvSyn: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param TypeEv The type of evidence-synthesis approach: Equal-evidence approach (0) or Added-evidence approach (1).
#' In case of an equal-evidence approach, aggregating evidence from, say, 5 studies with n=100 observations is the same as obtaining evidence from 1 study (as if it was possible) with n=500 observations (like meta-analysis does).
#' In the added-evidence approach, the aggregated evidence from, says, 5 studies is stronger than as if the data were combined (as if that was possible).
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param LL A matrix with log likelihood values of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained). Notably, only when the set of hypotheses cover the whol space / all theories (e.g., positive versus negative effect), then you can do without a safeguard hypothesis.
#' @param PT A matrix with penalty values of size S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained). Notably, only when the set of hypotheses cover the whol space / all theories (e.g., positive versus negative effect), then you can do without a safeguard hypothesis.
#' @param Name_studies Optional. Vector of S numbers or S characters to be printed at the x-axis of the plot with GORIC(A) weights. Default: Name_studies = 1:S.
#' @param Name_Hypo Optional. Vector containing 'NrHypos+1' characters which will be used for labelling the hypothesis. Default: H1, H2, ....
#' @param PrintPlot Optional. Indicator whether plot of GORIC(A) weigths should be printed (TRUE; default) or not (FALSE). The GORIC(A) weights per study are plotted and the cumulative GORIC(A) weights (where those for the last study are the final ones).
#'
#' @return The output comprises, among other things, the cumulative and final evidence for the theory-based hypotheses.
#' @export
#' @examples
#'
#' S <- 4
#' # Example based on S = 4 studies and 3 hypotheses:
#' # H0 <- "beta1 == 0"  # this hypothesis could have been left out
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
#'
#' # Change labels on x-axis in GORIC(A) weigths plot and give names to hypotheses #
#' # For example, let us say that the studies come from the years 2015, 2016, 2017, 2019.
#' # Because of unequal spacing, you may want to use numbers instead of characters:
#' Name_studies <- c(2015, 2016, 2017, 2019)
#' Name_Hypo <- c("H0", "Hpos", "Hneg")
#' GoricEvSyn_LLandPT(TypeEv, S, LL, PT, Name_studies, Name_Hypo)


GoricEvSyn_LLandPT <- function(TypeEv, S, LL, PT, Name_studies = 1:S, Name_Hypo = NULL, PrintPlot = T) {

  # Checks on input
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


  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  #CumulativeGorica <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  #CumulativeGoricaWeights <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  #colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(LL) <- colnames(PT) <- Name_Hypo
  #rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- rownames(LL) <- rownames(PT) <- paste0("Study", 1:S)
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(LL) <- colnames(PT) <- colnames(weight_m) <- Name_Hypo
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Study", 1:S), "Final")
  rownames(LL) <- rownames(PT) <- rownames(weight_m) <- paste0("Study", 1:S)


  sumLL <- 0
  sumPT <- 0
  IC <- -2 * LL + 2 * PT
  if(TypeEv == 1){ # added-ev approach
    for(s in 1:S){
      minIC <- min(IC[s,])
      weight_m[s,] <- exp(-0.5*(IC[s,]-minIC)) / sum(exp(-0.5*(IC[s,]-minIC)))
      #
      sumLL <- sumLL + LL[s,]
      sumPT <- sumPT + PT[s,]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT
      #CumulativeGoricaWeights[s,] <- exp(-0.5*CumulativeGorica[s,]) / sum(exp(-0.5*CumulativeGorica[s,]))
      minGoric <- min(CumulativeGorica[s,])
      CumulativeGoricaWeights[s,] <- exp(-0.5*(CumulativeGorica[s,]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s,]-minGoric)))
    }
    EvSyn_approach <- "Added-evidence approach"
  }else{ # equal-ev approach
    for(s in 1:S){
      sumLL <- sumLL + LL[s,]
      sumPT <- sumPT + PT[s,]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      #CumulativeGoricaWeights[s,] <- exp(-0.5*CumulativeGorica[s,]) / sum(exp(-0.5*CumulativeGorica[s,]))
      minGoric <- min(CumulativeGorica[s,])
      CumulativeGoricaWeights[s,] <- exp(-0.5*(CumulativeGorica[s,]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s,]-minGoric)))
    }
    EvSyn_approach <- "Equal-evidence approach"
  }

  CumulativeGorica[(S+1),] <- CumulativeGorica[S,]
  CumulativeGoricaWeights[(S+1),] <- CumulativeGoricaWeights[S,]
  #
  #Final.GORICA <- matrix(CumulativeGorica[S,], nrow = 1)
  Final.GORICA.weights <- CumulativeGoricaWeights[S,]
  Final.rel.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  #Final.GORICA.weights <- matrix(Final.GORICA.weights, nrow = 1)
  #rownames(Final.GORICA) <- "Final"
  #rownames(Final.GORICA.weights) <- "Final"
  rownames(Final.rel.GORICA.weights) <- Name_Hypo
  colnames(Final.rel.GORICA.weights) <- paste0("vs ", Name_Hypo)


  # Plot
  if(PrintPlot == T){
    NrHypos_incl <- (NrHypos + 1)
    Legend <- c("per study", "cumulative", Name_Hypo)
    Pch <- c(16, 8, rep(NA, NrHypos_incl))
    Col <- c(1, 1, 1:NrHypos_incl)
    Lty <- c(NA, 1, rep(1,NrHypos_incl))
    while (!is.null(dev.list()))  dev.off() # to reset the graphics pars to defaults
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
  #final <- list(LL_m = LL, PT_m = PT, GORICA_m = IC, GORICA.weight_m = weight_m,
  #              EvSyn_approach = EvSyn_approach, Cumulative.GORICA = CumulativeGorica, Cumulative.GORICA.weights = CumulativeGoricaWeights,
  #              Final.GORICA = Final.GORICA, Final.GORICA.weights = Final.GORICA.weights, Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  final <- list(LL_m = LL, PT_m = PT, GORICA_m = IC, GORICA.weight_m = weight_m,
                  EvSyn_approach = EvSyn_approach, Cumulative.GORICA = CumulativeGorica, Cumulative.GORICA.weights = CumulativeGoricaWeights,
                  Final.rel.GORICA.weights = Final.rel.GORICA.weights)
  return(final)

}
