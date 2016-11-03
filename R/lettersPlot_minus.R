#'@importFrom graphics par
lettersPlot_minus <- function (pwm)
{
  # pwm - negative Werte der Ergebnisse aus der Funktion hotspot()
  # ist nur mit der Kombination mit der Funktion letersPlot() verwendbar
  #
  # Ergebnis - Plot der Sequenz-logos in der Umgebung um Mutation
  #
  #

  chars <- c("A", "C", "G", "T")
  letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
  npos <- ncol(pwm)

  #ymin<-min(pwm)
  # ylim <- max(pwm)

  ylab <- "Probability"
  facs <- rep(1, npos)
  wt <- 1
  x.pos <- min(pwm)
  for (j in 1:npos) {
    column <- pwm[, j]
    hts <- 0.95 * column * facs[j]
    letterOrder <- order(hts, decreasing = TRUE)
    y.pos <- min(pwm[, j])
    pos1<-1+length(letters$y)
    for (i in 1:4) {
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht > 0)
        letters <- addLetter(letters, letter, x.pos,
                             y.pos, ht, wt)
      y.pos <- y.pos + ht + 0.01

    }
    pos2<-length(letters$y)
    letters$y[pos1:pos2]<-letters$y[pos1:pos2]-max(letters$y[pos1:pos2])
    x.pos <- x.pos + wt
  }


  grid.polygon(x = unit(letters$x, "native"),
               y = unit(letters$y, "native"),
               id = letters$id, gp = gpar(fill = letters$fill,
                                          col = "transparent"))
  popViewport()
  popViewport()
  par(ask = FALSE)
}
