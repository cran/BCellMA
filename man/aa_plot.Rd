% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aa_plot.R
\name{aa_plot}
\alias{aa_plot}
\title{Plot of average frequency of nucleotide mutations resulting in amino acid substitution}
\usage{
aa_plot(data, text, legend.position, characteristics)
}
\arguments{
\item{data}{Results from the function aa_dist or aa_cdr3_dist.}

\item{text}{The legend of the plot.}

\item{legend.position}{The position of the legend. It can be "none","left", "right", "bottom", "top".}

\item{characteristics}{is a data set "data(Klassen)" in this package}
}
\value{
Output is the plot of average frequency of nucleotide mutations resulting in amino acid substitution. The classes are divided according to IMGT.
}
\description{
The average frequency of nucleotide mutations resulting in amino acid substitution can be plotted as an aa_plot().
}
\examples{
data(IMGTtab7)
Regions<-cbind(IMGTtab7$FR1_IMGT,IMGTtab7$CDR1_IMGT)
allRegions_matrix<-aa_dist(Regions)
data(Klassen)
aa_plot(allRegions_matrix, "Amino acid Distribution", "right", Klassen)
}
\references{
Pommie C., Levadoux S., Sabatier R., Lefranc G., and Lefranc MP. IMGT standardizedcriteria for statistical analysis of immunoglobulin V-REGION amino acid properties.J MolRecognit, 17(1):17-32, 2004.
}
