#' Function to calculate the procentual difference at 6 positions to a nucleotide mutation to identify the hotspot motifs
#'
#'@param data Output from the function hotspotmutMat().
#'
#'@return Output is a list of four mutation environments A, C, T and G with the percentage difference of the nucleotides A, C, T and G at six positions to a mutation.
#'
#'@references Spencer J. and Dunn-Walters DK.  Hypermutation at A-T base pairs: the A nucleotidereplacement spectrum is affected by adjacent nucleotides and there is no reverse comple-mentary of sequences flanking muated A and T nucleotides.J Immunol, 175(8):5170 - 5177,2005.
#'@references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#'
#'@examples ##data(IMGTtab2) data(IMGTtab7)
#' germline<-germlineReconstr(IMGTtab2$V_REGION, IMGTtab7$V_REGION)
#' data<-targetingMatrix(data_tab2=IMGTtab2, data_tab_germline=germline, data_tab7=IMGTtab7)
#' targeting_motiv_data<-targeting_motiv(data)
#' targeting_motiv_data
#'@export
targeting_motiv<-function(data){

  # data       - Ergebnisse aus der Funktion hotspotmutMat()
  # Ergebnis   - eine Liste von vier Mutationsumgebungen A, C,
  #              T und G mit prozentuale Differenz der Nukleotiden
  #              A, C, T und G an sechs Positionen um eine Mutation


  A<-data$mutzur_A
  T<-data$mutzur_T
  C<-data$mutzur_C
  G<-data$mutzur_G

  Ag<-data$germ_A
  Tg<-data$germ_T
  Cg<-data$germ_C
  Gg<-data$germ_G

  for(n in 1:6){
    for(i in 1:4){

      A[i,n]<-(A[i,n])/sum(A[,n])
      T[i,n]<-T[i,n]/sum(T[,n])
      C[i,n]<-C[i,n]/sum(C[,n])
      G[i,n]<-G[i,n]/sum(G[,n])


      Ag[i,n]<-Ag[i,n]/sum(Ag[,n])
      Tg[i,n]<-Tg[i,n]/sum(Tg[,n])
      Cg[i,n]<-Cg[i,n]/sum(Cg[,n])
      Gg[i,n]<-Gg[i,n]/sum(Gg[,n])

      if(A[i,n]==0 & Ag[i,n]==0){ A[i,n]=1
      Ag[i,n]=1}

      if(T[i,n]==0 & Tg[i,n]==0){ T[i,n]=1
      Tg[i,n]=1}
      if(C[i,n]==0 & Cg[i,n]==0){ C[i,n]=1
      Cg[i,n]=1}

      if(G[i,n]==0 & Gg[i,n]==0){ G[i,n]=1
      Gg[i,n]=1}

    }
  }

  # Wichtig ! - es gibt eine Aufteilung in zwei Faelle bei der
  #             Division durch Null.
  # 1.Fall rel. Anteil an Base in der mut. Sequenz = 0 und
  #        rel. Anteil an Base in der Keimb. Sequenz. = 0
  #        Das Ergebnis der Division ist 0
  # 2.Fall rel. Anteil an Base in der mut. Sequenz != 0 und
  #        rel. Anteil an Base in der Keimb. Sequenz. = 0
  #        Das Ergebnis der Divison ist  rel. Anteil (%)
  #        an Base in der mut. Sequenz

  for(n in 1:6){
    for(i in 1:4){
      if(Ag[i,n]==0){
        A[i,n]=A[i,n]*100
      }else{A[i,n]=((A[i,n]-Ag[i,n])/Ag[i,n])*100}
      if(Tg[i,n]==0){
        T[i,n]=T[i,n]*100
      }else{T[i,n]=((T[i,n]-Tg[i,n])/Tg[i,n])*100}
      if(Cg[i,n]==0){
        C[i,n]=C[i,n]*100
      }else{C[i,n]=((C[i,n]-Cg[i,n])/Cg[i,n])*100}
      if(Gg[i,n]==0){
        G[i,n]=G[i,n]*100
      }else{G[i,n]=((G[i,n]-Gg[i,n])/Gg[i,n])*100}
    }
  }


  return(list(A=A,
              T=T,
              C=C,
              G=G))

}
