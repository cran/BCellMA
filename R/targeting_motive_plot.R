#' Polt of the difference at three positions before and after a mutation to identify the hotspot motifs.
#'
#' The the difference (%) at three positions before and after a mutation can be plotted using the function targeting_motive_plot().
#'@param pwm Result of the function hotspot().
#'@param xaxis Accept TRUE/FALSE parameter. TRUE draw X-axis.
#'@param yaxis Accept TRUE/FALSE parameter. TRUE draw Y-axis.
#'@param xfontsize The plot size of the X-axis.
#'@param yfontsize The plot size of the Y-axis.
#'@param xlim Limit of the drowed Y-axis.
#'
#'@note The Function are based on functions from package segLogo.
#'
#'@return Output is plot the sequence logos around the mutation.
#'
#'@references Spencer J. and Dunn-Walters DK.  Hypermutation at A-T base pairs: the A nucleotidereplacement spectrum is affected by adjacent nucleotides and there is no reverse comple-mentary of sequences flanking muated A and T nucleotides.J Immunol, 175(8):5170 - 5177,2005.
#'@references Zuckerman NS., Hazanov H., Barak M., Edelman H., Hess S., Shcolnik H., Dunn-Walters D.,and Mehr R. Somatic hypermutation and antigen-driven selection of B cells are altered inautoimmune diseases.J Autoimmun, 35(4):325 - 335, 2010. doi: 10.1016/j.jaut.2010.07.004.
#' @references Bembom O. seqLogo: Sequence logos for DNA sequence alignments, Status 10.08.2016. URLhttp://www.bioconductor.org/packages/release/bioc/html/seqLogo.html.
#' @references Schneider TD. and Stephens RM. Sequence logos: a new way to display consensus sequences.Nucleic Acids Res, 18(20):6097 - 6100, 1990.
#'
#' @examples data(IMGTtab2)
#' data(IMGTtab7)
#' germline<-germlineReconstr(IMGTtab2$V_REGION, IMGTtab7$V_REGION)
#' data<-targetingMatrix(data_tab2=IMGTtab2, data_tab_germline=germline, data_tab7=IMGTtab7)
#' targeting_motiv_data<-targeting_motiv(data)
#' targeting_motive_plot(targeting_motiv_data$A, xfontsize = 15, yfontsize = 15, xlim=60 )
#'
#'@export
#'@import grid
targeting_motive_plot <- function (pwm , xaxis = TRUE, yaxis = TRUE, xfontsize = 15,
                         yfontsize = 15, xlim)
{
  # pwm -   Ergebnis aus der Funktion hotspot()
  # xaxis - kann werte TRUE/FALSE annehmen. TRUE Plotet X-Achse.
  # yaxis - kann werte TRUE/FALSE annehmen. TRUE Plotet Y-Achse.
  # xfontsize - die Plot-Groesse der X-Achse
  # yfontsize - die Plot-Groesse der Y-Achse
  # xlim - Limit der Y-Achse, der geplotet werden soll
  #
  # Ergebnis - Plot der Sequenz-logos in der Umgebung um Mutation
  #
  #
  pwm_pos<-ifelse(as.matrix(pwm) >0, as.matrix(pwm), 0)

  chars <- c("A", "C", "G", "T")
  letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
  npos <- ncol(pwm_pos)

  ylab <- "% difference"
  facs <- rep(1, npos)

  wt <- 1
  x.pos <- min(pwm_pos)
  for (j in 1:npos) {
    column <- pwm_pos[, j]
    hts <- 0.95 * column * facs[j]
    letterOrder <- order(hts)
    y.pos <- min(pwm_pos[, j])
    for (i in 1:4) {
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht > 0)
        letters <- addLetter(letters, letter, x.pos,
                             y.pos, ht, wt)
      y.pos <- y.pos + ht + 0.01
    }
    x.pos <- x.pos + wt
  }
  grid.newpage()
  bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 2)
  leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 2)
  pushViewport(plotViewport(c(bottomMargin, leftMargin, 2,
                              2)))
  pushViewport(dataViewport(0:ncol(pwm_pos), -xlim:xlim, name = "vp1"))
  grid.polygon(x = unit(letters$x, "native"),
               y = unit(letters$y,"native"), id = letters$id,
               gp = gpar(fill = letters$fill, col = "transparent"))
  if (xaxis) {
    grid.xaxis(at = seq(0.5, ncol(pwm_pos) - 0.5),
               label = colnames(pwm_pos),
               gp = gpar(fontsize = xfontsize))
    grid.text("", y = unit(-3, "lines"),
              gp = gpar(fontsize = xfontsize))
  }
  if (yaxis) {
    grid.yaxis(gp = gpar(fontsize = yfontsize))
    grid.text(ylab, x = unit(-3, "lines"), rot = 90,
              gp = gpar(fontsize = yfontsize))
  }
  #popViewport()
  #popViewport()
  #par(ask = FALSE)

  grid.lines(c(.5,.5), c(.9,.1))
  grid.lines(c(.9,.1), c(.5,.5))

  pwm_neg<-ifelse(as.matrix(pwm)<0, as.matrix(pwm),0)
  pwm_neg<-abs(ifelse(pwm_neg ==0, 0.001, pwm_neg))
  lettersPlot_minus(pwm_neg)

}
