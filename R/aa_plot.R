#' Plot of Average mutation frequency of nucleotide mutations resulting in amino acid substitution
#' @param data Results from the function aa_dist or aa_cdr3_dist.
#' @param text The legend of the plot.
#' @param legend.position The position of the legend. It can be "none","left", "right", "bottom", "top".
#' @param characteristics is a data set data(Klassen) in this package
#' @return Output is the plot of average mutation frequency of nucleotide mutations resulting in amino acid substitution. The classes are divided according to IMGT.
#' @references Pommie C., Levadoux S., Sabatier R., Lefranc G., and Lefranc MP. IMGT standardizedcriteria for statistical analysis of immunoglobulin V-REGION amino acid properties.J MolRecognit, 17(1):17-32, 2004.
#' @examples ## data(IMGTtab7)
#' Regions<-cbind(IMGTtab7$FR1_IMGT,IMGTtab7$CDR1_IMGT)
#' \dontrun{allRegions_matrix<-aa_dist(Regions)}
#' ## data(Klassen)
#' \dontrun{aa_plot(allRegions_matrix, "Amino acid Distribution", "right", Klassen)}
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2 utils
#'
aa_plot<-function(data, text, legend.position, characteristics){



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

  longData<- cbind(longData,characteristics)

  colnames(longData)<-c("Var1", "Var2", "value",  "Eigenschaften", "Klassen")
  #colnames(longData)<-c("Var1", "Var2", "value")

  for(i in 1:400){
    if(longData$value[i] != 0){
     longData$Klassen[i] = longData$Klassen[i]

    }
   else{ longData$Klassen[i] = NA }
   }


  for (i in 1:400) {
    if(longData$Var1[i] == longData$Var2[i]){
      longData$value[i] = "same AA"
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

    color_palette<-  c("0%"="white",
                       "<1%"="#eff3ff",
                       "1-5%"="#bdd7e7",
                       "5-10%"="#6baed6",
                       "10-20%"="#3182bd",
                       "20-50%"="#08519c",
                       "50-100%"="#006837",
                       "same AA"="black")




    zp1 <- ggplot(longData, aes(x = longData$Var2, y = longData$Var1)) + geom_tile(aes(fill = factor(longData$value)),
                                                             colour = "black")
    zp1 <- zp1 +   scale_fill_manual(values = color_palette,
                                     name = "(Average \n mutation frequency")


    zp1 <- zp1 + geom_point(aes(shape = factor(longData$Klassen)))
    zp1 <- zp1 + scale_shape_manual(name = "classes",
                                    labels = c("very similar "," similar",
                                               "different", "very different"),
                                    values=c("1"=3,"2"=2, "3"=0, "4"=1) )

    zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
    zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
    zp1 <- zp1 + coord_equal()
    zp1 <- zp1 + theme_bw()
    zp1 <- zp1 + theme(legend.position = legend.position)
    zp1 <- zp1 + labs(list(title = text, y = "from", x="to"))
    print(zp1)

    #ggdraw(switch_axis_position(zp1, axis = 'x'))

  return(zp1)
}
