#' Average frequency of nucleotide mutations resulting in amino acid substitution in CDR3
#'
#' Calculate average frequency of nucleotide mutations resulting in amino acid substitution in CDR3.
#' @param data A column CDR3_IMGT from IMGT/HighV-Quest table 7 called "7_V-REGION-mutation-and-AA-change-table.txt"
#'
#' @return Output is the data matrix of average frequency of nucleotide mutations resulting in amino acid substitution.
#'
#' @references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#' @examples data(IMGTtab7)
#' cdr3_matrix<-aa_cdr3_dist(data=IMGTtab7$CDR3_IMGT)
#'
#' @export
#'
#'
aa_cdr3_dist<-function(data){

  # data      - CDR3_IMGT Spalte in der IMGT-Tabelle 7
  # Ergebnis  - Durchschnittliche Haeufigkeit der Aminosaeurenaustausche
  #
  data<-ifelse(data=="", NA, data)
  matrix_neu<-matrix(rep(0, 20*20), ncol=20, byrow=20)

  for(i in 1:length(data)){

    matrix<-aa_zaehlen2_CDR3(data[i])

    # nur bei 1000, 2000.. usw Daten ueber aktuelle Position ausgeben
    # wie weit die Auswertung ist
    # %% rest berechnung (modulo rechnen)12%10 = 2, 10%%10 = 0
    if (i%%1000 == 0 ){
      print(i)
    }
    if(sum(matrix)>0){
      matrix_skal<- matrix/sum(matrix)
      matrix_neu <- matrix_neu + matrix_skal
    }else{

    }

  }
  return(matrix_neu/length(data))
}
