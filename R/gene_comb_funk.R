#' gene/gene ratios of two gene families
#'
#' The calculation was obtained according bcRep package in R language. Every combination of two families was summarized and divided by the sum of all combinations across each group.
#' @param family1 "V_GENE_and_allele" colums from IMGT/HighV-Quest table 1 coled "1_Summary.txt"
#' @param family2 "J_GENE_and_allele" or "D_GENE_and_allele" colums from IMGT/HighV-Quest table 1 coled "1_Summary.txt"
#' @return Output is the data matrix of gene/gene ratios of two gene families.
#'
#' @references Bischof J. and Ibrahim SM. bcRep: R Package for Comprehensive Analysis of B Cell Receptor Repertoire Data. PLoS One. 11(8):e0161569, 2016. doi: 10.1371/journal.pone.0161569.
#' @examples data(IMGTtab1)
#' gane_comb<-gene_comb_funk(family1 = IMGTtab1$V_GENE_and_allele,
#'                           family2 = IMGTtab1$J_GENE_and_allele)
#' gane_comb
#'
#' @export
#'

gene_comb_funk<-function (family1,
                          family2){

  out.list <- list()

  family <- substr(strsplit(family1[which(family1 != "")][1],
                            split = " ")[[1]][2], 1, 4)

  gene<-substr(family1,1,20)
  genelist1 <- unlist(apply(data.frame(gene), 1, function(x) {
    strsplit(x, split = " |,|[.]|;|[|]|_|[*]|-|/")[[1]]
  }))

  das_wird_gebraucht <- genelist1[grep(family, genelist1)]
  das_wird_gebraucht <- substr(das_wird_gebraucht, 5, 6)
  das_wird_gebraucht<- as.numeric(das_wird_gebraucht)

  family_2 <- substr(strsplit(family2[which(family2 != "")][1],
                            split = " ")[[1]][2], 1, 4)

  gene<-substr(family2,1,20)
  genelist2 <- unlist(apply(data.frame(gene), 1, function(x) {
    strsplit(x, split = " |,|[.]|;|[|]|_|[*]|-|/")[[1]]
  }))
  das_wird_gebraucht_2 <- genelist2[grep(family_2, genelist2)]
  das_wird_gebraucht_2 <-substr(das_wird_gebraucht_2, 5, 6)
  das_wird_gebraucht_2<- as.numeric(das_wird_gebraucht_2)
  matrix<-matrix(rep(0, 7*7), ncol=7, byrow=7)
  colnames(matrix)<-paste(family,c(1:7),sep="")
  rownames(matrix)<-paste(family_2,c(1:7),sep="")


  for (i in 1:length(das_wird_gebraucht)) {

    matrix[das_wird_gebraucht_2[i],das_wird_gebraucht[i]]<-matrix[das_wird_gebraucht_2[i],das_wird_gebraucht[i]]+1

  }

  matrix<-matrix/sum(matrix)
  return(matrix)
}
