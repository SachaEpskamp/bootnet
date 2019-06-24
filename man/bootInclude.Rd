\name{bootInclude}
\alias{bootInclude}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inclusion proportion graph
}
\description{
This function takes bootstrap results and returns a inclusion probability network (edge weights indicate how often a certain edge was included in the model). Note that the plotting method automatically uses a black-white color scheme (as edges are not signed and always positive).
}
\usage{
bootInclude(bootobject, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bootobject}{
Nonparametric bootstrap results from \code{\link{bootnet}}
}
  \item{verbose}{
Logical, should progress be reported to the console?
}
}

\value{
A \code{bootnetResult} object with the following elements:
\item{graph}{The weights matrix of the network}
\item{intercepts}{The intercepts}
\item{results}{The results of the estimation procedure}
\item{labels}{A vector with node labels}
\item{nNodes}{Number of nodes in the network}
\item{nPerson}{Number of persons in the network}
\item{input}{Input used, including the result of the default set used}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{bootnet}}, \code{\link{estimateNetwork}}
}
\examples{
\dontrun{
# BFI Extraversion data from psychTools package:
library("psychTools")
data(bfi)
# Subset of data:
bfiSub <- bfi[1:250,1:25]

# Estimate ggmModSelect networks (not stepwise to increase speed):
Network <- estimateNetwork(bfiSub], default = "ggmModSelect", corMethod = "cor",
              stepwise = FALSE)

# Bootstrap 100 values, using 8 cores (100 to incease speed, preferably 1000+):
boots <- bootnet(Network, nBoots = 100, nCores = 8)

# Threshold network:
Network_inclusion <- bootInclude(boots)

# Plot:
plot(Network_inclusion)
}
}
