#'Function to the calculation of the number of Nucleotids at three positions before and after a mutation
#'
#'@param data_tab2 Column V_Region in th IMGT table 2
#'@param data_tab_germline Output from function germlineReconstr()
#'@param data_tab7 Column V_Region in th IMGT table 7
#'
#'@return Result is a list of 8 matrices with the numbers of nucleotides at three positions before and after a mutation from A, T, C, and G, as well as in a mutated sequence and in the germline sequence
#'
#'@references Spencer J. and Dunn-Walters DK.  Hypermutation at A-T base pairs: the A nucleotidereplacement spectrum is affected by adjacent nucleotides and there is no reverse comple-mentary of sequences flanking muated A and T nucleotides.J Immunol, 175(8):5170 - 5177,2005.
#'@references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#'
#'@examples data(IMGTtab2)
#' data(IMGTtab7)
#' germline<-germlineReconstr(IMGTtab2$V_REGION, IMGTtab7$V_REGION)
#' data<-targetingMatrix(data_tab2=IMGTtab2, data_tab_germline=germline, data_tab7=IMGTtab7)
#' data
#'@export
targetingMatrix<-function(data_tab2, data_tab_germline, data_tab7){

  # data_tab2         - IMGT-Tabelle 2
  # data_tab_germline - Keimbahngensequenzen
  # data_tab7         - IMGT-Tabelle 7
  #
  # Ergebnis ist eine Liste aus 8 Matrizen mit der Anzahlen an Nukleotiden
  # an sechs Positionen um eine Mutation A, T, C, und G
  # sowie in einer Mutierten Sequenz als auch in der
  # Keimbahngensequenz

  #Funktionsdefinition
  #Konvertierung von Buchstaben in Ascii Werte a=97 z=122
  #
  data_tab2<-ifelse(data_tab2$V_REGION=="", NA, data_tab2$V_REGION)
  data_tab7<-ifelse(data_tab7$V_REGION=="", NA, data_tab7$V_REGION)

  asc <- function(x) { strtoi(charToRaw(x),16L)}
  chr <- function(n) { rawToChar(as.raw(n)) }

  # zuerst werden die Matrizen bestehend aus Nullen definiert
  # z.B. ist matrixMutSeqA ein Matrix mit Anzahl an Nukleotiden
  # in der Umgebung um eine Mutation A
  MutSeqA<-matrix(rep(0, 6*213), ncol=6, byrow=2)
  MutSeqT<-matrix(rep(0, 6*232), ncol=6, byrow=2)
  MutSeqC<-matrix(rep(0, 6*215), ncol=6, byrow=2)
  MutSeqG<-matrix(rep(0, 6*219), ncol=6, byrow=2)

  GerSeqA<-matrix(rep(0, 6*213), ncol=6, byrow=2)
  GerSeqT<-matrix(rep(0, 6*232), ncol=6, byrow=2)
  GerSeqC<-matrix(rep(0, 6*215), ncol=6, byrow=2)
  GerSeqG<-matrix(rep(0, 6*219), ncol=6, byrow=2)

  for(n in 1:length(data_tab_germline)){


    germSeq<-strsplit(data_tab_germline[n],"" )
    germSeq<-matrix(unlist(germSeq), ncol=1, byrow=TRUE)


    mutSeq<-strsplit(data_tab2[n],"" )
    mutSeq<-matrix(unlist(mutSeq), ncol=1, byrow=TRUE)


    mutZahl<-matrix(rep(0, 1*length(mutSeq)), ncol=1, byrow=2)
    germZahl<-matrix(rep(0, 1*length(germSeq)), ncol=1, byrow=2)


    # nur bei 1000, 2000.. usw print ausgabe machen
    # wie weit die Auswertung ist
    # %% rest berechnung (modulo rechnen)12%10 = 2, 10%%10 = 0
    if (n%%1000 == 0 ){
      print(n)
    }
    # Schleife geht jede Base in der mutierten Sequenz durch
    # und wandelt diese in Zahlen um

    for(k in 1:length(mutSeq)){
      mutZahl[k]<-asc(mutSeq[k])
      germZahl[k]<-asc(germSeq[k])
    }

    mutationen <- data_tab7[n]

    temp2<- strsplit(mutationen," |,|[|]|[()]|/ " )
    temp3<-matrix(unlist(temp2), ncol=1, byrow=TRUE)


    for (i in 1:length(temp3)){


      if( nchar(temp3[i])>=3) {
        erste_buchstabe<-substr(temp3[i],1,1)
        #letzte Position herausfinden
        lastpos <- nchar(temp3[i])
        zweite_buchstabe<-substr(temp3[i],lastpos,lastpos)

        #Buchstaben in Zahlen konvertieren
        erste_wert = asc(erste_buchstabe)
        zweite_wert = asc(zweite_buchstabe)

        #Nur klein?buchstaben akzeptieren
        if (erste_wert >=97 && erste_wert <=122){

          temp4<-strsplit(temp3[i]," |>|/ " )
          temp5<-matrix(unlist(temp4), ncol=1, byrow=TRUE)

          # Position in der Sequenz bestimmen
          pos_a <- nchar(temp5[1])
          pos <-substr(temp5[1], 2, pos_a)


          #fuer buchstabe a
          if(erste_wert == 97){
            # Bbuchstaben an der Position -3 aufzehlen
            MutSeqA[mutZahl[as.numeric(pos)-3], 1]<-MutSeqA[mutZahl[as.numeric(pos)-3], 1]+1
            GerSeqA[germZahl[as.numeric(pos)-3],1]<-GerSeqA[germZahl[as.numeric(pos)-3],1]+1

            # Buchstaben an der position -2
            MutSeqA[mutZahl[as.numeric(pos)-2], 2]<-MutSeqA[mutZahl[as.numeric(pos)-2], 2]+1
            GerSeqA[germZahl[as.numeric(pos)-2],2]<-GerSeqA[germZahl[as.numeric(pos)-2],2]+1

            # Buchstaben an der stelle -1
            MutSeqA[mutZahl[as.numeric(pos)-1], 3]<-MutSeqA[mutZahl[as.numeric(pos)-1], 3]+1
            GerSeqA[germZahl[as.numeric(pos)-1],3]<-GerSeqA[germZahl[as.numeric(pos)-1],3]+1

            # Buchstaben an der Stelle +1
            MutSeqA[mutZahl[as.numeric(pos)+1], 4]<-MutSeqA[mutZahl[as.numeric(pos)+1], 4]+1
            GerSeqA[germZahl[as.numeric(pos)+1],4]<-GerSeqA[germZahl[as.numeric(pos)+1],4]+1

            # Buchstaben an der Stelle +2
            MutSeqA[mutZahl[as.numeric(pos)+2], 5]<-MutSeqA[mutZahl[as.numeric(pos)+2], 5]+1
            GerSeqA[germZahl[as.numeric(pos)+2],5]<-GerSeqA[germZahl[as.numeric(pos)+2],5]+1

            # Buchstaben an der Stelle +3
            MutSeqA[mutZahl[as.numeric(pos)+3], 6]<-MutSeqA[mutZahl[as.numeric(pos)+3], 6]+1
            GerSeqA[germZahl[as.numeric(pos)+3],6]<-GerSeqA[germZahl[as.numeric(pos)+3],6]+1

          }
          #fuer buchstabe t
          if(erste_wert == 116){
            # Buchstaben an der Position -3 aufzehlen
            MutSeqT[mutZahl[as.numeric(pos)-3], 1]<-MutSeqT[mutZahl[as.numeric(pos)-3], 1]+1
            GerSeqT[germZahl[as.numeric(pos)-3],1]<-GerSeqT[germZahl[as.numeric(pos)-3],1]+1

            # Buchstaben an der position -2

            MutSeqT[mutZahl[as.numeric(pos)-2], 2]<-MutSeqT[mutZahl[as.numeric(pos)-2], 2]+1
            GerSeqT[germZahl[as.numeric(pos)-2],2]<-GerSeqT[germZahl[as.numeric(pos)-2],2]+1

            # Buchstaben an der stelle -1

            MutSeqT[mutZahl[as.numeric(pos)-1], 3]<-MutSeqT[mutZahl[as.numeric(pos)-1], 3]+1
            GerSeqT[germZahl[as.numeric(pos)-1],3]<-GerSeqT[germZahl[as.numeric(pos)-1],3]+1

            # Buchstaben an der Stelle +1
            MutSeqT[mutZahl[as.numeric(pos)+1], 4]<-MutSeqT[mutZahl[as.numeric(pos)+1], 4]+1
            GerSeqT[germZahl[as.numeric(pos)+1],4]<-GerSeqT[germZahl[as.numeric(pos)+1],4]+1

            # Buchstaben an der Stelle +2
            MutSeqT[mutZahl[as.numeric(pos)+2], 5]<-MutSeqT[mutZahl[as.numeric(pos)+2], 5]+1
            GerSeqT[germZahl[as.numeric(pos)+2],5]<-GerSeqT[germZahl[as.numeric(pos)+2],5]+1

            # Buchstaben an der Stelle +3
            MutSeqT[mutZahl[as.numeric(pos)+3], 6]<-MutSeqT[mutZahl[as.numeric(pos)+3], 6]+1
            GerSeqT[germZahl[as.numeric(pos)+3],6]<-GerSeqT[germZahl[as.numeric(pos)+3],6]+1

          }

          #fuer buchstabe c
          if(erste_wert== 99){
            # Buchstaben an der Position -3 aufzehlen
            MutSeqC[mutZahl[as.numeric(pos)-3], 1]<-MutSeqC[mutZahl[as.numeric(pos)-3], 1]+1
            GerSeqC[germZahl[as.numeric(pos)-3],1]<-GerSeqC[germZahl[as.numeric(pos)-3],1]+1

            # Buchstaben an der position -2
            MutSeqC[mutZahl[as.numeric(pos)-2], 2]<-MutSeqC[mutZahl[as.numeric(pos)-2], 2]+1
            GerSeqC[germZahl[as.numeric(pos)-2],2]<-GerSeqC[germZahl[as.numeric(pos)-2],2]+1

            # Buchstaben an der stelle -1
            MutSeqC[mutZahl[as.numeric(pos)-1], 3]<-MutSeqC[mutZahl[as.numeric(pos)-1], 3]+1
            GerSeqC[germZahl[as.numeric(pos)-1],3]<-GerSeqC[germZahl[as.numeric(pos)-1],3]+1

            # Buchstaben an der Stelle +1
            MutSeqC[mutZahl[as.numeric(pos)+1], 4]<-MutSeqC[mutZahl[as.numeric(pos)+1], 4]+1
            GerSeqC[germZahl[as.numeric(pos)+1],4]<-GerSeqC[germZahl[as.numeric(pos)+1],4]+1

            # Buchstaben an der Stelle +2
            MutSeqC[mutZahl[as.numeric(pos)+2], 5]<-MutSeqC[mutZahl[as.numeric(pos)+2], 5]+1
            GerSeqC[germZahl[as.numeric(pos)+2],5]<-GerSeqC[germZahl[as.numeric(pos)+2],5]+1

            # Buchstaben an der Stelle +3
            MutSeqC[mutZahl[as.numeric(pos)+3], 6]<-MutSeqC[mutZahl[as.numeric(pos)+3], 6]+1
            GerSeqC[germZahl[as.numeric(pos)+3],6]<-GerSeqC[germZahl[as.numeric(pos)+3],6]+1

          }
          #fuer Buchstabe g
          if(erste_wert==103){

            # Buchstaben an der Position -3 aufzehlen
            MutSeqG[mutZahl[as.numeric(pos)-3], 1]<-MutSeqG[mutZahl[as.numeric(pos)-3], 1]+1
            GerSeqG[germZahl[as.numeric(pos)-3],1]<-GerSeqG[germZahl[as.numeric(pos)-3],1]+1

            # Buchstaben an der position -2
            MutSeqG[mutZahl[as.numeric(pos)-2], 2]<-MutSeqG[mutZahl[as.numeric(pos)-2], 2]+1
            GerSeqG[germZahl[as.numeric(pos)-2],2]<-GerSeqG[germZahl[as.numeric(pos)-2],2]+1

            # Buchstaben an der stelle -1
            MutSeqG[mutZahl[as.numeric(pos)-1], 3]<-MutSeqG[mutZahl[as.numeric(pos)-1], 3]+1
            GerSeqG[germZahl[as.numeric(pos)-1],3]<-GerSeqG[germZahl[as.numeric(pos)-1],3]+1

            # Bbuchstabe an der Stelle +1
            MutSeqG[mutZahl[as.numeric(pos)+1], 4]<-MutSeqG[mutZahl[as.numeric(pos)+1], 4]+1
            GerSeqG[germZahl[as.numeric(pos)+1],4]<-GerSeqG[germZahl[as.numeric(pos)+1],4]+1

            # Buchstaben an der Stelle +2
            MutSeqG[mutZahl[as.numeric(pos)+2], 5]<-MutSeqG[mutZahl[as.numeric(pos)+2], 5]+1
            GerSeqG[germZahl[as.numeric(pos)+2],5]<-GerSeqG[germZahl[as.numeric(pos)+2],5]+1

            # Buchstaben an der Stelle +3
            MutSeqG[mutZahl[as.numeric(pos)+3], 6]<-MutSeqG[mutZahl[as.numeric(pos)+3], 6]+1
            GerSeqG[germZahl[as.numeric(pos)+3],6]<-GerSeqG[germZahl[as.numeric(pos)+3],6]+1

          }
        }
      }
    }

  }

  matrixMutSeqA<-MutSeqA[c(97, 99, 103, 116), ]
  colnames(matrixMutSeqA)<- c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixMutSeqA)<- c("A", "C", "G", "T")
  matrixMutSeqT <- MutSeqT[c(97, 99, 103, 116), ]
  colnames(matrixMutSeqT)<- c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixMutSeqT)<- c("A", "C", "G", "T")
  matrixMutSeqC <- MutSeqC[c(97, 99, 103, 116), ]
  colnames(matrixMutSeqC)<-c("Pos-3", "Pos-2", "Pos-1",
                             "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixMutSeqC)<-c("A", "C", "G", "T")
  matrixMutSeqG<-MutSeqG[c(97, 99, 103, 116),]
  colnames(matrixMutSeqG)<-c("Pos-3", "Pos-2", "Pos-1",
                             "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixMutSeqG)<-c("A", "C", "G", "T")
  matrixGermSeqA<-GerSeqA[c(97, 99, 103, 116), ]
  colnames(matrixGermSeqA)<-c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixGermSeqA)<-c("A", "C", "G", "T")
  matrixGermSeqT<-GerSeqT[c(97, 99, 103, 116), ]
  colnames(matrixGermSeqT)<-c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixGermSeqT)<-c("A", "C", "G", "T")
  matrixGermSeqC<-GerSeqC[c(97, 99, 103, 116), ]
  colnames(matrixGermSeqC)<-c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixGermSeqC)<-c("A", "C", "G", "T")
  matrixGermSeqG<-GerSeqG[c(97, 99, 103, 116),]
  colnames(matrixGermSeqG)<-c("Pos-3", "Pos-2", "Pos-1",
                              "Pos+1", "Pos+2", "Pos+3")
  rownames(matrixGermSeqG)<-c("A", "C", "G", "T")

  return(list(mutzur_A=matrixMutSeqA, mutzur_T=matrixMutSeqT,
              mutzur_C= matrixMutSeqC, mutzur_G=matrixMutSeqG,
              germ_A=matrixGermSeqA, germ_T=matrixGermSeqT,
              germ_C=matrixGermSeqC, germ_G=matrixGermSeqG))
}
