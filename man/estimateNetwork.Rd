\name{estimateNetwork}
\alias{estimateNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a network structure
}
\description{
This function allows for flexible estimation of a network structure using various R packages and model frameworks. This is typically done by using one of the default sets. See details for manual specification. See also Epskamp, Borsboom and Fried (2016) for more details.
}
\usage{
estimateNetwork(data, default = c("none", "EBICglasso", "pcor", "IsingFit", "IsingLL", "huge", "adalasso"), prepFun, prepArgs = list(), estFun, estArgs, graphFun, graphArgs, intFun, intArgs, labels, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or matrix containing the raw data. Must be numeric, integer or ordered factors.
}
  \item{default}{
A string indicating the method to use. Specifying a \code{default} sets default values to \code{prepFun}, \code{prepArgs}, \code{estFun}, \code{estArgs}, \code{graphFun}, \code{graphArgs}, \code{intFun} and \code{intArgs}. Setting a \code{default} can be omitted but that does require specifying all above mentioned arguments. Current options are:
\describe{
\item{\code{"EBICglasso"}}{Gaussian Markov random field estimation using graphical LASSO and extended Bayesian information criterium to select optimal regularization parameter. Using \code{\link[qgraph]{EBICglasso}} from the qgraph package.}
\item{\code{"IsingFit"}}{Ising model estimation using LASSO regularized nodewise logistic regression and extended Bayesian information criterium to select optimal regularization parameter. Using \code{\link[IsingFit]{IsingFit}} from the IsingFit package.}
\item{\code{"pcor"}}{Partial correlation network (non-regularized Gaussian Markov random field), using \code{\link[corpcor]{cor2pcor}} from the corpcor package.}
\item{\code{adalasso}}{Uses the \code{\link[parcor]{adalasso.net}} from the parcor package.}
\item{\code{huge}}{Uses EBIC model selection of GGM networs estimated via the glasso algorithm as implemented in the huge package (as opposed to glasso and qgraph packages used in \code{default = "EBICglasso"}). Uses nonparanormal transformation in preparing the data and does not use polychoric correlations.}
}
See details section for a more detailed description.
}
  \item{prepFun}{
A function that takes as input the raw data and returns whatever the estimation function needs (as first argument). Typically this function is used to correlate or binarize the data. Defaults to \code{\link{identity}} if omitted and \code{default} is not set.
}
  \item{prepArgs}{
A list with arguments for \code{prepFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
}
  \item{estFun}{
A function that takes as input the prepared raw data and estimates a network.
}
  \item{estArgs}{
A list with arguments for \code{estArgs}
}
  \item{graphFun}{
A function that takes the result of \code{estFun} and extracts the estimated weights matrix. Defaults to \code{\link{identity}} if omitted and \code{default} is not set.
}
  \item{graphArgs}{
A list with arguments for \code{graphFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
}
  \item{intFun}{
A function that takes the result of \code{estFun} and extracts the estimated intercepts. Defaults to \code{\link{null}} if omitted and \code{default} is not set.
}
  \item{intArgs}{
A list with arguments for \code{intFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
}
  \item{labels}{
A character vector containing the node labels. If omitted the column names of the data are used.
}
  \item{verbose}{
Logical, currently only used when \code{default = "EBICglasso"} in the \code{cor_auto} function.
}
}
\details{
The user can manually set the estimation method of the PMRF in bootnet using a set of arguments to the bootnet function. First, the method of preprocessing the data must be defined via the \code{prepFun} argument, which must be assigned a function that takes a dataset as input and returns the viable input for the network estimator. The argument \code{prepArgs} can be specified a list of arguments to the prepFun function. Data preprocessing typically means correlating the data for the GGM or binarizing it for the Ising model (to this end bootnet provides a \code{binarize} function). Next, we estimate the network. To do so, we assign the \code{estFun} argument, a function that takes whatever the output of \code{prepFun} was and estimates a network. The \code{estArgs} argument can be used to assign a list of additional arguments to the function used in \code{estFun}. Finally, we need to extract the weights matrix and intercepts. Assigning functions to the \code{graphFun} and \code{intFun} arguments respectively can do this. An example of how these commands work together to estimate a network using \code{default = "EBICglasso"} is shown below:
	
\code{estimateNetwork(Data, 
	prepFun = cor_auto, 
	prepArgs =  list(verbose = FALSE), 
	estFun = qgraph::EBICglasso, 
	estArgs = list(n = nrow(bfi)), 
	graphFun = identity, 
	intFun= null)}
	
	The following defaults can be used:
\describe{
\item{\code{default = "EBICglasso"}}{
\describe{
\item{prepFun}{\code{\link[qgraph]{cor_auto}} from the qgraph package}
\item{prepArgs}{\code{list(verbose=TRUE)} for estimating a network and \code{list(verbose=FALSE)} for bootstrap samples.}
\item{estFun}{\code{\link[qgraph]{EBICglasso}} from the qgraph package}
\item{estArgs}{\code{list(n = nrow(data), returnAllResults = TRUE)}}
\item{graphFun}{\code{function(x)x[['optnet']]}}
\item{graphArgs}{\code{list()}}
\item{intFun}{\code{\link[bootnet]{null}}}
\item{intArgs}{\code{list()}}
}
}
\item{\code{default = "IsingFit"}}{
\describe{
\item{prepFun}{\code{\link[bootnet]{binarize}}}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{\link[IsingFit]{IsingFit}} from the IsingFit package}
\item{estArgs}{\code{list(plot = FALSE, progress = FALSE)}}
\item{graphFun}{\code{function(x)x[['weiadj']]}}
\item{graphArgs}{\code{list()}}
\item{intFun}{\code{function(x)x[['thresholds']]}}
\item{intArgs}{\code{list()}}
}
}

\item{\code{default = "pcor"}}{
\describe{
\item{prepFun}{\code{\link[qgraph]{cor_auto}} from the qgraph package}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{\link[corpcor]{cor2pcor}} from the corpcor package}
\item{estArgs}{\code{list()}}
\item{graphFun}{\code{\link{identity}}}
\item{graphArgs}{\code{list()}}
\item{intFun}{null}
\item{intArgs}{\code{list()}}
}
}

\item{\code{default = "adalasso"}}{
\describe{
\item{prepFun}{\code{identity}}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{\link[parcor]{adalasso.net}} from the parcor package}
\item{estArgs}{\code{list()}}
\item{graphFun}{\code{function(x)as.matrix(Matrix::forceSymmetric(x$pcor.adalasso)}}
\item{graphArgs}{\code{list()}}
\item{intFun}{null}
\item{intArgs}{\code{list()}}
}
}

\item{\code{default = "huge"}}{
\describe{
\item{prepFun}{\code{function(x)huge.npn(na.omit(as.matrix(x)),verbose = FALSE)}, using \code{\link[huge]{huge.npn}} from the huge package}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{function(x)huge::huge.select(huge::huge(x,method = "glasso",verbose=FALSE), criterion = "ebic",verbose = FALSE)}, using the huge package}
\item{estArgs}{\code{list()}}
\item{graphFun}{\code{as.matrix(qgraph::wi2net(as.matrix(x$opt.icov)))}, using the huge package}
\item{graphArgs}{\code{list()}}
\item{intFun}{null}
\item{intArgs}{\code{list()}}
}
}

}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
# BFI Extraversion data from psych package:
library("psych")
data(bfi)
bfiSub <- bfi[,1:25]

# Estimate network:
Network <- estimateNetwork(bfiSub, default = "EBICglasso")

# Some pointers:
print(Network)

# Estimated network:
plot(Network, layout = 'spring') 

\dontrun{
# Centrality indices:
library("qgraph")
centralityPlot(Network)
}
}
