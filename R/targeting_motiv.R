#' Function to calculate the difference at three positions before and after a mutation to identify the hotspot motifs
#'
#'The change around the mutated nucleotide flanking three bases each is expressed as the difference (%) between the substitution frequency in the mutated sequence and the substitution frequency in the germline sequence.
#'@param data Output from the function hotspotmutMat().
#'
#'@return Output is a list of four mutation environments A, C, T and G with the difference (%) of the nucleotides A, C, T and G at three positions before and after a mutation.
#'
#'@references Spencer J. and Dunn-Walters DK.  Hypermutation at A-T base pairs: the A nucleotidereplacement spectrum is affected by adjacent nucleotides and there is no reverse comple-mentary of sequences flanking muated A and T nucleotides.J Immunol, 175(8):5170 - 5177,2005.
#'@references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#'
#'@examples data(IMGTtab2)
#' data(IMGTtab7)
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


    A1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    T1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    C1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    G1<-matrix(rep(0, 4*6), ncol=6, byrow=4)

    Ag1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    Tg1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    Cg1<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    Gg1<-matrix(rep(0, 4*6), ncol=6, byrow=4)

    A2<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    colnames(A2)<- c("Pos-3", "Pos-2", "Pos-1", "Pos+1", "Pos+2", "Pos+3")
    rownames(A2)<- c("A", "C", "G", "T")
    T2<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    colnames(T2)<- c("Pos-3", "Pos-2", "Pos-1", "Pos+1", "Pos+2", "Pos+3")
    rownames(T2)<- c("A", "C", "G", "T")
    C2<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    colnames(C2)<- c("Pos-3", "Pos-2", "Pos-1", "Pos+1", "Pos+2", "Pos+3")
    rownames(C2)<- c("A", "C", "G", "T")
    G2<-matrix(rep(0, 4*6), ncol=6, byrow=4)
    colnames(G2)<- c("Pos-3", "Pos-2", "Pos-1", "Pos+1", "Pos+2", "Pos+3")
    rownames(G2)<- c("A", "C", "G", "T")


    for(i in 1:4){
      for(n in 1:6){


        A1[i,n]<-(A[i,n])/sum(A[,n])
        T1[i,n]<-T[i,n]/sum(T[,n])
        C1[i,n]<-C[i,n]/sum(C[,n])
        G1[i,n]<-G[i,n]/sum(G[,n])


        Ag1[i,n]<-Ag[i,n]/sum(Ag[,n])
        Tg1[i,n]<-Tg[i,n]/sum(Tg[,n])
        Cg1[i,n]<-Cg[i,n]/sum(Cg[,n])
        Gg1[i,n]<-Gg[i,n]/sum(Gg[,n])

        if(A[i,n]==0 & Ag[i,n]==0){ A1[i,n]=1
        Ag1[i,n]=1}

        if(T[i,n]==0 & Tg[i,n]==0){ T1[i,n]=1
        Tg1[i,n]=1}
        if(C[i,n]==0 & Cg[i,n]==0){ C1[i,n]=1
        Cg1[i,n]=1}

        if(G[i,n]==0 & Gg[i,n]==0){ G1[i,n]=1
        Gg1[i,n]=1}

      }
    }
    for(n in 1:6){
      for(i in 1:4){
        if(Ag1[i,n]==0){
          A2[i,n]=A1[i,n]*100
        }else{A2[i,n]=((A1[i,n]-Ag1[i,n])/Ag1[i,n])*100}
        if(Tg1[i,n]==0){
          T2[i,n]=T1[i,n]*100
        }else{T2[i,n]=((T1[i,n]-Tg1[i,n])/Tg1[i,n])*100}
        if(Cg1[i,n]==0){
          C2[i,n]=C1[i,n]*100
        }else{C2[i,n]=((C1[i,n]-Cg1[i,n])/Cg1[i,n])*100}
        if(Gg1[i,n]==0){
          G2[i,n]=G1[i,n]*100
        }else{G2[i,n]=((G1[i,n]-Gg1[i,n])/Gg1[i,n])*100}
      }
    }


    return(list(A=A2,
                T=T2,
                C=C2,
                G=G2))

  }
