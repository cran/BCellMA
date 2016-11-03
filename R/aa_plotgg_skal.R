#' @import reshape2 utils
aa_plotgg_skal<-function(data){

  # data      - Ergebnis der Funktion aa_verteilung()

  # Ergebnis  - Daten fuer Aminosaeurenaustausche Plot,
  #             Einteilung in physiochemische Eigenschaften,
  #             Einteilung in prozentuale Klassen

  data<-ifelse((data) >0 & data <=0.01, 0.01, data)
  data<-ifelse((data)>0.01 & data <=0.05, 0.05, data)
  data<-ifelse((data)>0.05 & data <=0.1, 0.1, data)
  data<-ifelse((data)>0.1 & data <=0.2, 0.2, data)
  data<-ifelse((data)>0.2 & data <=0.5, 0.5, data)
  data<-ifelse((data)>0.5, 1, data)


  longData <- melt(data)
  Klassen<-read.csv("BCellMA/R/klasse1.csv",
                    header = T, sep=";")

  longData<- cbind(longData,Klassen)

  colnames(longData)<-c("Var1", "Var2", "value",  "Eigenschaften",
                        "Klassen")


  for(i in 1:400){
    if(longData$value[i] != 0){
      longData$Klassen[i] = longData$Klassen[i]

    }
    else{ longData$Klassen[i] = NA }
  }


  for (i in 1:400) {
    if(longData$Var1[i] == longData$Var2[i]){
      longData$value[i] = "gleiche AA"
    }
    else{longData$value[i] =longData$value[i]}
  }

  for (i in 1:400) {
    if(longData$value[i] == 0){
      longData$value[i] = "0%"
    }
    else{longData$value[i] =longData$value[i]}
  }

  for (i in 1:400) {
    if(longData$value[i] == 0.01){
      longData$value[i] = "<1%"
    }
    else{longData$value[i] =longData$value[i]}
  }

  for (i in 1:400) {
    if(longData$value[i] == 0.05){
      longData$value[i] = "1-5%"
    }
    else{longData$value[i] =longData$value[i]}
  }

  for (i in 1:400) {
    if(longData$value[i] == 0.1){
      longData$value[i] = "5-10%"
    }
    else{longData$value[i] =longData$value[i]}
  }

  for (i in 1:400) {
    if(longData$value[i] == 0.2){
      longData$value[i] = "10-20%"
    }
    else{longData$value[i] =longData$value[i]}
  }
  for (i in 1:400) {
    if(longData$value[i] == 0.5){
      longData$value[i] = "20-50%"
    }
    else{longData$value[i] =longData$value[i]}
  }
  for (i in 1:400) {
    if(longData$value[i] == 1){
      longData$value[i] = "50-100%"
    }
    else{longData$value[i] =longData$value[i]}
  }
  return(longData)
}
