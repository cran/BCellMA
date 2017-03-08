#' Function to reconstruction of germline sequencies based on IMGT outputs.
#'
#'@param data_tab2 Column V_Region in the IMGT table 2
#'@param data_tab7 Column V_Region in the IMGT table 7
#'
#'@return output is a vector of germline sequencies.
#'
#' @references Brochet X., Lefranc MP., and Giudicelli V. IMGT/V-QUEST: the highly customized andintegrated system for IG and TR standardized V-Jand V-D-J sequence analysis.NucleicAcids Res., 36(Web Server issue):W503 - W508, 2008. doi: 10.1093/nar/gkn316
#' @examples data(IMGTtab2)
#' data(IMGTtab7)
#' germline<-germlineReconstr(IMGTtab2$V_REGION, IMGTtab7$V_REGION)
#' germline
#'@export
#'
germlineReconstr<- function(data_tab2, data_tab7){

  # data_tab2 - Daten aus der IMGT-Tabelle 2
  # data_tab7 - Daten aus der IMGT-Tabelle 7

  # Ergebnis  - Keimbahngen

  # zuerst wird eine Matrix bestehend aus Nullen erzeugt
  vector_germline<-matrix(rep(0, 1*length(data_tab2)), ncol=1, byrow=2)
  data_tab7<-ifelse(data_tab7 =="", 1, data_tab7)


  # Konvertierung von Buchstaben in Ascii Werte a=97 z=122
  asc <- function(x) { strtoi(charToRaw(x),16L) }
  chr <- function(n) { rawToChar(as.raw(n)) }

  for (n in 1:length(data_tab2)) {
    S<-data_tab2[n]
    S<-strsplit(S, "")
    S<-matrix(unlist(S), ncol=1, byrow=TRUE)

    temp2<- strsplit(data_tab7[n]," |,|[|]|[()]|/ " )
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

        #Nur Gro?buchstaben akzeptieren
        if (erste_wert >=97 && erste_wert <=122){
          temp4<-strsplit(temp3[i]," |>|/ " )
          temp5<-matrix(unlist(temp4), ncol=1, byrow=TRUE)
          # Position in der Sequenz bestimmen
          position <- nchar(temp5[1])
          position_in_sequenz<-substr(temp5[1], 2,position)

          # Falls die mutierende Base in den Keimbahngen nicht der Base
          # in Mutierten Sequenz entspricht,
          # soll ein NA fuer eine spaetere Ueberpr?ffung eingebaut werden.

          b <- ifelse(S[as.numeric(position_in_sequenz),1]==zweite_buchstabe,
                      erste_buchstabe, NA )
          # die Base wird in Mutierenden Sequenzn umgeschrieben
          S[as.numeric(position_in_sequenz,1)]<-b
        }
        else{
          # i eins hoeher zaehlen und nichts tun
        }

      }
    }
    vector_germline[n]<-paste(S, collapse="")

    # nur bei 1000, 2000.. usw print ausgabe machen
    # wie weit die Auswertung ist
    if (n%%1000 == 0 ){
      print(n)
    }

  }
  return(vector_germline)
}
