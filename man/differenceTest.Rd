\name{differenceTest}
\alias{differenceTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrapped difference test
}
\description{
This function computes the bootstrapped difference test for edge-weights and centrality indices. A confidence interval is constructured on the difference of two values, and the test is deemed significant if zero is not in this confidence interval. See also Epskamp, Borsboom and Fried (2016) for more details.
}
\usage{
differenceTest(bootobject, x, y, measure = c("strength", "closeness", "betweenness"),
                alpha = 0.05, x2, y2, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bootobject}{
Output of \code{\link{bootnet}}. Must be nonparametric or parametric bootstrap.
}
  \item{x}{
A character string indicating the ID of a node or an edge, or a number indiacting the node or edge. For an edge, can be e.g., \code{"1--2"} or \code{"x = 1, x2 = 2"}.
}
  \item{y}{
A character string indicating the ID of a node or an edge, or a number indiacting the node or edge. For an edge, can be e.g., \code{"1--2"} or \code{"y = 1, y2 = 2"}.
}
  \item{measure}{
Measure to test. Can be \code{"strength"}, \code{"closeness"}, \code{"betweenness"}, \code{"edge"} or \code{"distance"}.
}
  \item{alpha}{
Signifiance level to test at. Note that the actual signifiance level is influenced by the number of bootstrap samples, and is returned in a message.
}
  \item{x2}{
Second node in an edge. optional.
}
  \item{y2}{
Second node in an edge. optional.
}
  \item{verbose}{
Logical, should the message indiacting actual signifiance level be printed?
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
Results1 <- bootnet(Network, nBoots = 1000, nCores = 8)

# Test for difference in strength between node "A1" and "C2":
differenceTest(Results,  "A1", "C2", "strength")

# Test for difference between edge N1--N2 and N3--N4:
differenceTest(Results,  "N1--N2", "N3--N4", "edge")

# Alternative:
differenceTest(Results,  x = "N1", x2 = "N2", y = "N3", 
                y2 = "N4", measure = "edge")
}
}
