\name{ddCTcalculate}
\alias{ddCTcalculate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ddCTcalculate(geneOfInterest, sampleOfInterest = "020_N", houseKeepingGene = "GAPDH", referenceSample = "020_N", data = rawData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneOfInterest}{
%%     ~~Describe \code{geneOfInterest} here~~
}
  \item{sampleOfInterest}{
%%     ~~Describe \code{sampleOfInterest} here~~
}
  \item{houseKeepingGene}{
%%     ~~Describe \code{houseKeepingGene} here~~
}
  \item{referenceSample}{
%%     ~~Describe \code{referenceSample} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (geneOfInterest, sampleOfInterest = "020_N", houseKeepingGene = "GAPDH", 
    referenceSample = "020_N", data = rawData) 
{
    rawData$ddCt_GAP_020_N = ddCTcalculate(geneOfInterest = rawData$gene.x, 
        sampleOfInterest = rawData$origin.x, houseKeepingGene = "GAPDH", 
        referenceSample = "020_N", data = rawData)
    sampleHouse = paste(sampleOfInterest, houseKeepingGene)
    SampleGene = paste(sampleOfInterest, geneOfInterest)
    houseCp = data[sampleHouse, "meanCP"]
    geneCp = data[SampleGene, "meanCP"]
    dCt = houseCp - geneCp
    refDctRowHouse = paste(referenceSample, houseKeepingGene)
    refDctRowGene = paste(referenceSample, geneOfInterest)
    referenceSample_dCt = data[refDctRowHouse, "meanCP"] - data[refDctRowGene, 
        "meanCP"]
    ddCtNotSquared = dCt - referenceSample_dCt
    ddCt = 2^ddCtNotSquared
    return(ddCt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
