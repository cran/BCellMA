% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lengthCDR3.R
\name{lengthCDR3}
\alias{lengthCDR3}
\title{Frequency distribution of the CDR3 length}
\usage{
lengthCDR3(data)
}
\arguments{
\item{data}{"CDR3_IMGT_length" colums from IMGT/HighV-Quest table 1 coled "1_Summary.txt"}
}
\value{
Output Output is the number of amino acids in CDR3.
}
\description{
Calculate length of CDR3.
}
\examples{
data(IMGTtab1)
lenght_tab<-lengthCDR3(as.numeric(IMGTtab1$CDR3_IMGT_length))
}
\references{
references Bischof J. and Ibrahim SM. bcRep: R Package for Comprehensive Analysis of B Cell Receptor Repertoire Data. PLoS One. 11(8):e0161569, 2016. doi: 10.1371/journal.pone.0161569.
}
