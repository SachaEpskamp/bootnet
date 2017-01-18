\name{estimateNetwork}
\alias{estimateNetwork}
\alias{bootnet_EBICglasso}
\alias{bootnet_pcor}
\alias{bootnet_IsingFit}
\alias{bootnet_IsingSampler}
\alias{bootnet_adalasso}
\alias{bootnet_huge}
\alias{bootnet_mgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a network structure
}
\description{
This function allows for flexible estimation of a network structure using various R packages and model frameworks. This is typically done by using one of the default sets. See details for manual specification. See also Epskamp, Borsboom and Fried (2016) for more details.
}
\usage{
estimateNetwork(data, default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", 
  "huge","adalasso","mgm"), fun, prepFun, prepArgs, estFun, estArgs, 
  graphFun, graphArgs, intFun, intArgs, labels, verbose = TRUE,  construct = 
  c("default","function","arguments"), .dots = list(), weighted = TRUE, signed = TRUE,
  \dots, .input)
                
bootnet_EBICglasso(data, tuning = 0.5,corMethod = c("cor_auto","cov","cor","npn"), 
  missing = c("pairwise","listwise","stop"), sampleSize = c("maximum","minimim"),
  verbose = verbose, corArgs = list())
  
bootnet_pcor(data, corMethod = c("cor_auto","cov","cor","npn"), missing = c("pairwise",
  "listwise","stop"), verbose = verbose, corArgs = list())
  
bootnet_IsingFit(data, tuning = 0.25, missing = c("listwise","stop"), verbose = verbose,
  rule = c("AND","OR"), split = "median")
  
bootnet_IsingSampler(data, missing = c("listwise","stop"), verbose = verbose,
  split = "median", method = c("default","ll","pl","uni","bi"))
  
bootnet_adalasso(data,missing = c("listwise","stop"), verbose = verbose, nFolds = 10)
  
bootnet_huge(data, tuning = 0.5, missing = c("listwise","stop"), verbose = verbose,
  npn = TRUE, criterion = c("ebic","ric","stars"))
  
bootnet_mgm(data, type, lev, tuning = 0.5, missing = c("listwise","stop"), 
  verbose = verbose, criterion = c("EBIC","CV"), nFolds = 10, degree = 2,
  rule = c("AND","OR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or matrix containing the raw data. Must be numeric, integer or ordered factors.
}
  \item{default}{
A string indicating the method to use. Specifying a \code{default} sets default values to \code{prepFun}, \code{prepArgs}, \code{estFun}, \code{estArgs}, \code{graphFun}, \code{graphArgs}, \code{intFun} and \code{intArgs}. Setting a \code{default} can be omitted but that does require specifying all above mentioned arguments. Current options are:
\describe{
\item{\code{"EBICglasso"}}{Gaussian Markov random field estimation using graphical LASSO and extended Bayesian information criterion to select optimal regularization parameter. Using \code{\link[qgraph]{EBICglasso}} from the qgraph package. Calls \code{bootnet_EBICglasso}.}
\item{\code{"IsingFit"}}{Ising model estimation using LASSO regularized nodewise logistic regression and extended Bayesian information criterion to select optimal regularization parameter. Using \code{\link[IsingFit]{IsingFit}} from the IsingFit package.  Calls \code{bootnet_IsingFit}.}
\item{\code{IsingSampler}}{Calls the \code{\link[IsingSampler]{EstimateIsing}} function from the IsingSampler package. }
\item{\code{"pcor"}}{Partial correlation network (non-regularized Gaussian Markov random field), using \code{\link[corpcor]{cor2pcor}} from the corpcor package.  Calls \code{bootnet_pcor}.}
\item{\code{adalasso}}{Uses the \code{\link[parcor]{adalasso.net}} from the parcor package. Calls \code{bootnet_adalasso}.}
\item{\code{huge}}{Uses EBIC model selection of GGM networks estimated via the glasso algorithm as implemented in the huge package (as opposed to glasso and qgraph packages used in \code{default = "EBICglasso"}). Uses nonparanormal transformation in preparing the data and does not use polychoric correlations. Calls \code{bootnet_huge}.}
\item{\code{mgm}}{Estimates a Mixed graphical model by using the \code{\link[mgm]{mgmfit}} function of the mgm package. Calls \code{bootnet_mgm}.}
}
See details section for a more detailed description.
}
\item{fun}{
A custom estimation function, when no default set is used. This must be a function that takes the data as input (first argument) and returns either a weights matrix or a list containing the elements \code{"graph"} for the weights matrix, \code{"intercepts"} for the intercepts (optional) and \code{"results"} for the full estimation results (optional).
}
\item{tuning}{EBIC tuning parameter, used in 'EBICglasso', 'IsingFit', 'huge' and 'mgm' default sets.}
\item{corMethod}{Correlation method, used in 'EBICglasso' and 'pcor' default sets. \code{"cor_auto"} uses \code{\link[qgraph]{cor_auto}} for polychoric and polyserial correlations, \code{"cov"} uses the \code{\link{cov}} function for covariances, \code{"cor"} will use the \code{\link{cor}} function for correlations and \code{"npn"} will apply the nonparanormal transformation (via \code{\link[huge]{huge.npn}}) and then compute correlations.}
\item{missing}{How to handle missing data? \code{"pairwise"} for pairwise deletion, \code{"listwise"} for listwise deletion and \code{"stop"} to stop with an error.}
\item{sampleSize}{How will sample size be computed in EBICglasso default set? \code{"maximum"} takes total number of rows including rows with NA, and \code{"minimum"} takes total number of rows that contain no NA.}
\item{corArgs}{A list with arguments for the function used defined by \code{corMethod}.}
\item{rule}{The rule to use to select an edge in nodewise estimation. \code{"AND"} to only select in edge if both regression coefficients are nonzero and \code{"OR"} if only one is nonzero. Used in 'IsingFit' and 'mgm' default sets.}
\item{split}{A function or character string (\code{"median"} or \code{"mean"}) indicating how to binarize values when estimating an Ising model.}
\item{method}{The estimation method used in the IsingSampler default set (see \code{\link[IsingSampler]{EstimateIsing}}).}
\item{npn}{Logical, should nonparanormal be used in huge default set?}
\item{criterion}{The criterion used in model selection. \code{"ebic"}, \code{"ric"} or \code{"stars"} in the huge default set or \code{"EBIC"} or \code{"CV"} in the mgm default set.}
\item{nFolds}{Number of folds used in k-fold cross-validation.}
\item{type}{See \code{\link[mgm]{mgmfit}}.}
\item{lev}{See \code{\link[mgm]{mgmfit}}.}
\item{degree}{Degree of mixed model, see \code{\link[mgm]{mgmfit}}}
  \item{prepFun}{
A function that takes as input the raw data and returns whatever the estimation function needs (as first argument). Typically this function is used to correlate or binarize the data. Defaults to \code{\link{identity}} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{prepArgs}{
A list with arguments for \code{prepFun}. Defaults to \code{list()} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{estFun}{
A function that takes as input the prepared raw data and estimates a network. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{estArgs}{
A list with arguments for \code{estArgs}.  ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{graphFun}{
A function that takes the result of \code{estFun} and extracts the estimated weights matrix. Defaults to \code{\link{identity}} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{graphArgs}{
A list with arguments for \code{graphFun}. Defaults to \code{list()} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{intFun}{
A function that takes the result of \code{estFun} and extracts the estimated intercepts. Defaults to \code{\link{null}} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{intArgs}{
A list with arguments for \code{intFun}. Defaults to \code{list()} if omitted and \code{default} is not set. ARGUMENT IS DEPRECATED AND WILL BE REMOVED.
}
  \item{labels}{
A character vector containing the node labels. If omitted the column names of the data are used.
}
  \item{verbose}{
Logical, currently only used when \code{default = "EBICglasso"} in the \code{cor_auto} function.
}

\item{weighted}{Logical, should the analyzed network be weighted?}
\item{signed}{Logical, should the analyzed network be signed?}
\item{construct}{
Should the estimator function be constructed using new bootnet functionality (\code{function}) or old bootnet functionality (\code{arguments})? Will default to \code{function} unless arguments such as \code{estFun} are used.
}
\item{.dots}{
A list of arguments used in the estimation function set by a default set or by the \code{fun} argument.
}
\item{\dots}{
A list of arguments used in the estimation function set by a default set or by the \code{fun} argument.
}
\item{.input}{Used internally in the bootnet function. Do not use.}
}
\details{
NEW FUNCTIONALTIY:

The user can manually specify an estimation method by assinging a custom function to the \code{'fun'} argument. This function must take data as input and output an estimated network.


DEPRECATED OLD FUNCTIONALITY:

The user can manually set the estimation method of the PMRF in bootnet using a set of arguments to the bootnet function. First, the method of preprocessing the data must be defined via the \code{prepFun} argument, which must be assigned a function that takes a dataset as input and returns the viable input for the network estimator. The argument \code{prepArgs} can be specified a list of arguments to the \code{prepFun} function. Data preprocessing typically means correlating the data for the GGM or binarizing it for the Ising model (to this end bootnet provides a \code{binarize} function). Next, we estimate the network. To do so, we assign the \code{estFun} argument, a function that takes whatever the output of \code{prepFun} was and estimates a network. The \code{estArgs} argument can be used to assign a list of additional arguments to the function used in \code{estFun}. Finally, we need to extract the weights matrix and intercepts. Assigning functions to the \code{graphFun} and \code{intFun} arguments respectively can do this. An example of how these commands work together to estimate a network using \code{default = "EBICglasso"} is shown below:
	
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

\dontrun{
# Some pointers:
print(Network)

# Estimated network:
plot(Network, layout = 'spring') 

# Centrality indices:
library("qgraph")
centralityPlot(Network)

# BIC model selection:
Network_BIC <- estimateNetwork(bfiSub, default = "EBICglasso", tuning = 0)

# Ising model:
Network_BIC <- estimateNetwork(bfiSub, default = "IsingFit")
}
}
