#' Frequency distribution of the CDR3 length
#'
#' Calculate length of CDR3.
#' @param data "CDR3_IMGT_length" colums from IMGT/HighV-Quest table 1 coled "1_Summary.txt"
#' @return Output Output is the number of amino acids in CDR3.
#' @references references Bischof J. and Ibrahim SM. bcRep: R Package for Comprehensive Analysis of B Cell Receptor Repertoire Data. PLoS One. 11(8):e0161569, 2016. doi: 10.1371/journal.pone.0161569.
#' @examples data(IMGTtab1)
#' lenght_tab<-lengthCDR3(as.numeric(IMGTtab1$CDR3_IMGT_length))
#' @export

lengthCDR3<-function (data){

  out.list <- list()

  minCDR3 <- min(as.numeric(data), na.rm = T)
  maxCDR3 <- max(as.numeric(data), na.rm = T)

  tab.bar <- t(data.frame(apply(data.frame(seq(minCDR3, maxCDR3,
                                               1)), 1, function(x) {
                                                 length(which(data == x))
                                               })))
  rownames(tab.bar)<-c("")
  colnames(tab.bar) <- paste(" ", as.character(seq(minCDR3,
                                                   maxCDR3, 1)), sep = "")

  tab.bar <- tab.bar/length(data)*100

  return(tab.bar)
}
