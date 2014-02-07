\name{extractReplicates}
\alias{extractReplicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractReplicates(indexes, ctData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexes}{
%%     ~~Describe \code{indexes} here~~
}
  \item{ctData}{
%%     ~~Describe \code{ctData} here~~
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
function (indexes, ctData) 
{
    indexes = c(1:384)
    CtData = data[complete.cases(data[, 3]), ]
    even = indexes[indexes\%\%2 == 0]
    odd = indexes[indexes\%\%2 == 1]
    rep1 = CtData[odd, c(1:6)]
    rep1 = na.omit(rep1)
    rep2 = CtData[even, c(1:6)]
    rep2 = na.omit(rep2)
    boundData = merge(rep1, rep2, by.x = "sample", by.y = "sample")
    usefulData = boundData[, c(1, 2, 3, 4, 6, 7, 11)]
    usefulData$meanCP = rowMeans(cbind(usefulData$Cp.x, usefulData$Cp.y))
    usefulData$stdDevCP = apply(cbind(usefulData$Cp.x, usefulData$Cp.y), 
        1, sd)
    result = list(rep1, rep2, usefulData)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
