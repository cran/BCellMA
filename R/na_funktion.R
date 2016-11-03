na_funktion<-function(data){

  # data - Ergebnisse der germlineReconstr()
  # Ergebniss- Anzahl an NA in Keimbahngenen

  n<-rep(1, length(data))

  for(i in 1:length(data)){

    split<- strsplit(data[i]," || " )
    temp3<-matrix(unlist(split), ncol=1, byrow=TRUE)

    summme<-sum(is.na(temp3))
    n[i]<-summme

  }
  return((sum(n)))

}
