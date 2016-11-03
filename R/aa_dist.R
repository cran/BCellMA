#' Average mutation frequency of nucleotide mutations resulting in amino acid substitution in FRs and CDRs
#'
#' @param data CDRs and FRs colums from IMGT/HighV-Quest table 7 coled "7_V-REGION-mutation-and-AA-change-table.txt"
#'
#' @return Output is the data matrix of Average mutation frequency of nucleotide mutations resulting in amino acid substitution.
#'
#' @references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#' @examples ## data(IMGTtab7)
#' Regions<-cbind(IMGTtab7$FR1_IMGT,IMGTtab7$CDR1_IMGT)
#' \dontrun{Regions_matrix<-aa_dist(data=Regions)}
#' \dontrun{Regions_matrix}
#' @export
#'
aa_dist<-function(data){

  # data      - Daten aus der IMGT-Tabelle 7, die CDRs und FRs
  # Ergebnis  - Durchschnittliche Haeufigkeit der Aminosaeurenaustausche
  #

  data[,1]<-ifelse(data[,1]=="", NA, data[,1])
  data[,2]<-ifelse(data[,2]=="", NA, data[,2])
  matrix_neu<-matrix(rep(0, 20*20), ncol=20, byrow=20)

  for(i in 1:length(data[,1])){

    matrix<-aa_zaehlen2(data[i,])


    # nur bei 1000, 2000.. usw print ausgabe machen
    # wie weit die Auswertung ist
    # %% rest Berechnung (modulo rechnen)12%10 = 2, 10%%10 = 0
    if (i%%1000 == 0 ){
      print(i)
    }
    if(sum(matrix)>0){
      matrix_skal<- matrix/sum(matrix)
      matrix_neu <- matrix_neu + matrix_skal
    }else{

    }

  }
  return(matrix_neu/length(data[,1]))
}
