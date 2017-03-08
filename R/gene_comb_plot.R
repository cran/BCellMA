#' Plot of gene/gene ratios of two gene families
#'
#' Gene/gene ratios of two gene families can be plotted as an gene_comb_plot().
#'
#' @param data Results from the function gene_comb_funk().
#' @param text The legend of the plot.
#' @param legend_position The position of the legend. It can be "none","left", "right", "bottom", "top".
#' @param a The angle of X-axis legend.
#' @param b The distance from X-axis to legend.
#'
#' @return Output is the plot of gene/gene ratios of two gene families.
#'
#' @references Bischof J. and Ibrahim SM. bcRep: R Package for Comprehensive Analysis of B Cell Receptor Repertoire Data. PLoS One. 11(8):e0161569, 2016. doi: 10.1371/journal.pone.0161569.
#' @examples data(IMGTtab1)
#' gane_comb<-gene_comb_funk(family1 = IMGTtab1$V_GENE_and_allele,
#'                           family2 = IMGTtab1$J_GENE_and_allele)
#' gene_comb_plot(gane_comb, "Plot of IGHV and IGHJ ratio", legend_position = "right", a = 35, b = 0.5)
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2 utils
#'



gene_comb_plot<-function(data, text, legend_position, a, b){

  data <- melt(data)
  colnames(data)<-c("Var1", "Var2", "value")
  zp1 <- ggplot(data,
                aes(x = data$Var2, y = data$Var1)) + geom_tile(aes(fill = data$value), colour = "white")
  zp1 <- zp1 +   scale_fill_gradient(low = "white", high = "black", name ="%")

  zp1 <- zp1 + labs(x = "", y = "")
  zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
  zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
  zp1 <- zp1 + coord_equal()
  zp1 <- zp1 + theme_bw()
  zp1 <- zp1 + theme(legend.position = legend_position, axis.text.x = element_text(angle = a, vjust = b))
  zp1 <- zp1 + labs(list(title = text))

  print(zp1)
}
