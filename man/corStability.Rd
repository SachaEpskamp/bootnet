\name{corStability}
\alias{corStability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation stability coefficient
}
\description{
This coefficient denotes the estimated maximum number of cases that can be dropped from the data to retain, with 95\% probability, a correlation of at least 0.7 (default) between statistics based on the original network and statistics computed with less cases. This coefficient should not be below 0.25 and is preferably above 0.5. See also Epskamp, Borsboom and Fried (2016) for more details.
}
\usage{
corStability(x, cor = 0.7, statistics = c("strength", "closeness", "betweenness"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output of \code{\link{bootnet}}. Must be case-drop bootstrap.
}
  \item{cor}{
The correlation level tot est at.
}
  \item{statistics}{
The statistic(s) to test for.
}
}
\references{
Epskamp, S., Borsboom, D., & Fried, E. I. (2016). Estimating psychological networks and their accuracy: a tutorial paper. arXiv preprint, arXiv:1604.08462.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{bootnet}}
}

\examples{
\dontrun{
# BFI Extraversion data from psych package:
library("psych")
data(bfi)
bfiSub <- bfi[,1:25]

# Estimate network:
Network <- estimateNetwork(bfiSub, default = "EBICglasso")

# Bootstrap 1000 values, using 8 cores:
# Bootstrap 1000 values, using 8 cores:
Results2 <- bootnet(Network, nBoots = 1000, nCores = 8, 
              type = "case")

# Compute CS-coefficients:
corStability(Results2)
}
}
