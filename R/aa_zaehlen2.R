aa_zaehlen2<- function(data){

  # data      - Daten aus der IMGT-Tabelle 7, die CDRs und FRs
  # Ergebnis  - Anzahl der Aminosaeurenaustatusche
  #

  #Funktionsdefinition
  #Konvertierung von Buchstaben in Ascii Werte A=65 Z=90
  #

  asc <- function(x) { strtoi(charToRaw(x),16L) }
  chr <- function(n) { rawToChar(as.raw(n)) }

  matrix_neu<-matrix(rep(0, 26*26), ncol=26, byrow=26)

  for(i in 1:length(data)){

    if (is.na(data[i])==FALSE){


    temp2<- strsplit(data[i]," |,|[|]|[()]|/ " )
    temp3<-matrix(unlist(temp2), ncol=1, byrow=TRUE)

    for (i in 1:length(temp3)){


      if(nchar(temp3[i,1])>=3) {
        erste_buchstabe<-substr(temp3[i],1,1)
        #letzte Position herausfinden
        lastpos <- nchar(temp3[i])
        zweite_buchstabe<-substr(temp3[i],lastpos,lastpos)

        #Buchstaben in Zahlen konvertieren
        erste_wert = asc(erste_buchstabe)
        zweite_wert = asc(zweite_buchstabe)

        #Nur Gro?buchstaben akzeptieren
        if (erste_wert >=65 && erste_wert <=91 && zweite_wert>=65 && zweite_wert <=91){
          matrix_neu[erste_wert-64,zweite_wert-64] <- matrix_neu[erste_wert-64,
                                                                 zweite_wert-64] + 1
        }else{}
      }
      else{
        # i eins hoeher zaehlen und nichts tun
      }

    }


  }
  }
    matrix_neu <- matrix_neu[ ,c(4,5,11,8,18,7,14,17,1,19,20,
                                 16,3,22,13,25,9,12,6,23)]
    matrix_neu <- matrix_neu[ c(23,6,12,9,25,13,22,3,16,20,19,
                                1,17,14,7,18,8,11,5,4), ]
    colnames(matrix_neu)<-c("D", "E", "K", "H", "R", "G", "N", "Q","A", "S",
                            "T", "P", "C", "V", "M", "Y", "I", "L", "F", "W")
    rownames(matrix_neu)<-c("W", "F", "L", "I", "Y", "M", "V", "C", "P", "T",
                            "S", "A", "Q", "N", "G", "R", "H", "K", "E", "D")

  return(matrix_neu)
}

