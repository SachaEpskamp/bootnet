\name{estimateNetwork}
\alias{estimateNetwork}
\alias{bootnet_EBICglasso}
\alias{bootnet_pcor}
\alias{bootnet_cor}
\alias{bootnet_IsingFit}
\alias{bootnet_IsingSampler}
\alias{bootnet_adalasso}
\alias{bootnet_huge}
\alias{bootnet_mgm}
\alias{bootnet_relimp}
\alias{bootnet_TMFG}
\alias{bootnet_LoGo}
\alias{bootnet_ggmModSelect}
\alias{bootnet_graphicalVAR}
\alias{bootnet_SVAR_lavaan}
\alias{bootnet_piecewiseIsing}
\alias{bootnet_ncvRegularize}
\alias{bootnet_nodeRegresIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate a network structure
}
\description{
This function allows for flexible estimation of a network structure using various R packages and model frameworks. This is typically done by using one of the default sets. See details for manual specification. See also Epskamp, Borsboom and Fried (2016) for more details. IMPORTANT: THE ESTIMATOR FUNCTIONS (e.g., \code{fun = bootnet_pcor}) ARE NOT INTENDED TO BE USED MANUALY (see details).
}
\usage{
estimateNetwork(data, default = c("none", "EBICglasso", "pcor",
                 "IsingFit", "IsingSampler", "huge", "adalasso", "mgm",
                 "relimp", "cor", "TMFG", "ggmModSelect", "LoGo",
                 "graphicalVAR", "piecewiseIsing", "SVAR_lavaan",
                 "ncvRegularize", "nodeRegresIC"),
                 fun, labels, verbose = TRUE, .dots = list(), weighted = TRUE,
                 signed = TRUE, directed, datatype, checkNumeric = FALSE, ...,
                 .input, memorysaver = FALSE)

bootnet_EBICglasso(data, tuning = 0.5, corMethod = c("cor", "cov",
                    "cor_auto", "cor_mantar", "npn", "spearman"), missing =
                    c("pairwise", "listwise", "fiml", "stackedMI", "stop"), sampleSize
                    = c("pairwise_average", "maximum", "minimum",
                    "pairwise_maximum", "pairwise_minimum",
                    "pairwise_average_v1.5", "pairwise_maximum_v1.5",
                    "pairwise_minimum_v1.5"), verbose = TRUE, corArgs =
                    list(), refit = FALSE, principalDirection = FALSE,
                    lambda.min.ratio = 0.01, nlambda = 100, threshold =
                    FALSE, unlock = FALSE, nonPositiveDefinite = c("stop",
                    "continue"), transform = c("none", "rank",
                    "quantile"), ...)

bootnet_pcor(data, corMethod = c("cor", "cov", "cor_auto", "cor_mantar",
                    "npn", "spearman"), missing = c("pairwise", "listwise",
                    "fiml", "stackedMI", "stop"), sampleSize = c("pairwise_average",
                    "maximum", "minimum", "pairwise_maximum",
                    "pairwise_minimum", "pairwise_average_v1.5",
                    "pairwise_maximum_v1.5", "pairwise_minimum_v1.5"),
                    verbose = TRUE, corArgs = list(), threshold = 0, alpha
                    = 0.05, adjacency, principalDirection = FALSE, unlock
                    = FALSE, nonPositiveDefinite = c("stop", "continue"),
                    transform = c("none", "rank", "quantile"))

bootnet_cor(data, corMethod = c("cor", "cov", "cor_auto", "cor_mantar",
                    "npn", "spearman"), missing = c("pairwise", "listwise",
                    "fiml", "stackedMI", "stop"), sampleSize = c("pairwise_average",
                    "maximum", "minimum", "pairwise_maximum",
                    "pairwise_minimum", "pairwise_average_v1.5",
                    "pairwise_maximum_v1.5", "pairwise_minimum_v1.5"),
                    verbose = TRUE, corArgs = list(), threshold = 0, alpha
                    = 0.05, principalDirection = FALSE, unlock = FALSE,
                    nonPositiveDefinite = c("stop", "continue"), transform
                    = c("none", "rank", "quantile"))

bootnet_IsingFit(data, tuning = 0.25, missing = c("listwise", "stop"),
                    verbose = TRUE, rule = c("AND", "OR"), split =
                    "median", principalDirection = FALSE,
                    min_sum = -Inf, unlock = FALSE)

bootnet_IsingSampler(data, missing = c("listwise", "stop"), verbose = TRUE,
                    split = "median", method = c("uni", "ll", "pl", "bi"),
                    principalDirection = FALSE, unlock = FALSE, threshold
                    = FALSE, alpha = 0.01, min_sum = -Inf, rule = c("AND",
                    "OR"))

bootnet_adalasso(data, missing = c("listwise", "stop"), verbose = TRUE,
                   nFolds = 10, principalDirection = FALSE, unlock =
                   FALSE, transform = c("none", "rank", "quantile"), ...)

bootnet_huge(data, tuning = 0.5, missing = c("listwise", "stop"),
                   verbose = TRUE, npn = TRUE, criterion = c("ebic",
                   "ric", "stars"), principalDirection = FALSE,
                   lambda.min.ratio = 0.01, nlambda = 100, unlock =
                   FALSE, transform = c("none", "rank", "quantile"), ...)

bootnet_mgm(data, type, level, tuning = 0.25, missing =
                   c("listwise", "stop"), verbose = TRUE, criterion =
                   c("EBIC", "CV"), nFolds = 10, order = 2, rule =
                   c("AND", "OR"), binarySign, unlock = FALSE, transform
                   = c("none", "rank", "quantile"), ...)

bootnet_relimp(data, normalized = TRUE, type = "lmg",
                   structureDefault = c("none", "custom", "EBICglasso",
                   "pcor", "IsingFit", "IsingSampler", "huge",
                   "adalasso", "mgm", "cor", "TMFG", "ggmModSelect",
                   "LoGo", "nodeRegresIC", "ncvRegularize"),
                   missing = c("listwise", "stop"), ..., verbose
                   = TRUE, threshold = 0, unlock = FALSE, transform =
                   c("none", "rank", "quantile"))

bootnet_TMFG(data, graphType = c("cor", "pcor"), corMethod =
                   c("cor", "cov", "cor", "npn", "cor_auto", "cor_mantar"),
                   missing = c("pairwise", "listwise", "fiml", "stackedMI",
                   "stop"), verbose = TRUE, corArgs = list(), nonPositiveDefinite
                   = c("stop", "continue"), principalDirection = FALSE, unlock =
                   FALSE, transform = c("none", "rank", "quantile"), ...)

bootnet_LoGo(data, corMethod = c("cor", "cov", "cor", "npn", "cor_mantar",
                   "cor_auto"), missing = c("pairwise", "listwise",
                   "fiml", "stackedMI", "stop"), verbose = TRUE, corArgs = list(),
                   principalDirection = FALSE, unlock = FALSE, nonPositiveDefinite =
                   c("stop","continue"),transform =
                   c("none", "rank", "quantile"), ...)

bootnet_graphicalVAR(data, tuning = 0.5, verbose = TRUE, principalDirection
                   = FALSE, missing = c("listwise", "stop"), unlock =
                   FALSE, transform = c("none", "rank", "quantile"), ...)

bootnet_ggmModSelect(data, tuning = 0, corMethod = c("cor", "cov",
                    "cor_auto", "cor_mantar", "npn", "spearman"), missing =
                    c("pairwise", "listwise", "fiml", "stackedMI", "stop"),
                    sampleSize = c("pairwise_average", "maximum", "minimum",
                    "pairwise_maximum", "pairwise_minimum",
                    "pairwise_average_v1.5", "pairwise_maximum_v1.5",
                    "pairwise_minimum_v1.5"), verbose = TRUE, corArgs =
                    list(), principalDirection = FALSE, start =
                    c("glasso", "empty", "full"), stepwise = TRUE, nCores
                    = 1, unlock = FALSE, nonPositiveDefinite = c("stop",
                    "continue"), transform = c("none", "rank",
                    "quantile"), ...)

bootnet_piecewiseIsing(data, cutoff, missing = c("listwise", "stop"), verbose
                    = TRUE, IsingDefault = c("IsingSampler", "IsingFit",
                    "custom"), zeroThreshold = 1, minimalN = ncol(data) +
                    1, unlock = FALSE, ...)

bootnet_SVAR_lavaan(data, verbose = TRUE, principalDirection = FALSE,
                   missing = c("listwise", "stop"), criterion = "bic",
                   eqThreshold = 1e-04, tempWhitelist, tempBlacklist,
                   contWhitelist, contBlacklist, minimalModInd = 10,
                   unlock = FALSE, transform = c("none", "rank",
                   "quantile"), ...)

bootnet_ncvRegularize(data, penalty = c("atan", "selo", "exp", "log",
                    "sica", "scad", "mcp", "glasso"), corMethod = c("cor",
                    "cov", "cor_auto", "cor_mantar", "npn", "spearman"),
                    missing = c("pairwise", "listwise", "fiml", "stackedMI",
                    "stop"), sampleSize = c("pairwise_average", "maximum",
                    "minimum", "pairwise_maximum", "pairwise_minimum",
                    "pairwise_average_v1.5", "pairwise_maximum_v1.5",
                    "pairwise_minimum_v1.5"), verbose = TRUE, corArgs =
                    list(), principalDirection = FALSE, unlock = FALSE,
                    nonPositiveDefinite = c("stop", "continue"), transform
                    = c("none", "rank", "quantile"), ...)

bootnet_nodeRegresIC(data, rule = c("AND", "OR"), regressionSampleSize =
                    c("individual","average","max","total"), criterion =
                    c("bic","aic","aicc"), corMethod = c("cor",
                    "cov", "cor_auto", "cor_mantar", "npn", "spearman"),
                    missing = c("pairwise", "listwise", "fiml", "stackedMI",
                    "stop"), verbose = TRUE, corArgs =
                    list(), principalDirection = FALSE, unlock = FALSE,
                    nonPositiveDefinite = c("stop", "continue"), transform
                    = c("none", "rank", "quantile"), ...)
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
\item{\code{"ncvRegularize"}}{Estimates a GGM using nonconvex regularization with model selection based on an information criterion. Convex regularization via graphical LASSO is also supported. Uses \code{\link[mantar]{regularization_net}} from the mantar package. Calls \code{bootnet_ncvRegularize}.}
\item{\code{"IsingFit"}}{Ising model estimation using LASSO regularized nodewise logistic regression and extended Bayesian information criterion to select optimal regularization parameter. Using \code{\link[IsingFit]{IsingFit}} from the IsingFit package.  Calls \code{bootnet_IsingFit}.}
\item{\code{"IsingSampler"}}{Calls the \code{\link[IsingSampler]{EstimateIsing}} function from the IsingSampler package. }
\item{\code{"pcor"}}{Partial correlation network (non-regularized Gaussian Markov random field), using \code{\link[corpcor]{cor2pcor}} from the corpcor package.  Calls \code{bootnet_pcor}.}
\item{\code{"cor"}}{Correlation network.}
\item{\code{"adalasso"}}{Temporarily deprecated due to required functions not being available on CRAN: Uses the \code{\link[parcor]{adalasso.net}} from the parcor package. Calls \code{bootnet_adalasso}.}
\item{\code{"huge"}}{Uses EBIC model selection of GGM networks estimated via the glasso algorithm as implemented in the huge package (as opposed to glasso and qgraph packages used in \code{default = "EBICglasso"}). Uses nonparanormal transformation in preparing the data and does not use polychoric correlations. Calls \code{bootnet_huge}.}
\item{\code{"mgm"}}{Estimates a Mixed graphical model by using the the mgm (or mgmfit in older versions) function of the mgm package. Calls \code{bootnet_mgm}.}
\item{\code{"TMFG"}}{Estimates a Triangulated Maximally Filtered Graph, using the function \code{TMFG} of the NetworkToolbox package.  Calls \code{bootnet_TMFG}. Note that this estimates a *correlation network* by default (use the 'graphType' argument to estimate a partial correlation network instead).}
\item{\code{"LoGo"}}{Estimates a Local/Global Sparse Inverse Covariance Matrix, using the function \code{LoGo} of the NetworkToolbox package.  Calls \code{bootnet_LoGo}.}
\item{\code{"relimp"}}{Estimates a (directed) relative importance network, using the function 'calc.relimp' of the 'relaimpo' package. The 'structureDefault' argument can be used to use a different default set for estimating the structure of the graph. Calls \code{bootnet_relimp}.}
\item{\code{"ggmModSelect"}}{Estimates an unregularized GGM using the glasso algorithm and stepwise model selection, using the 'ggmModSelect' function from the qgraph package. Calls \code{bootnet_ggmModSelect}.}
\item{\code{"graphicalVAR"}}{Estimates a graphical VAR model using the graphicalVAR package. This results in two networks which can be plotted using the 'graph' argument in the plot method. Calls \code{bootnet_graphicalVAR}.}
\item{\code{"nodeRegresIC"}}{Estimates an unregularized GGM using neighborhood selection (nodewise multiple regressions) based on an information criterion. Uses \code{\link[mantar]{neighborhood_net}} from the mantar package. Calls \code{bootnet_nodeRegresIC}.}
}
See details section for a more detailed description.
}
\item{fun}{
A custom estimation function, when no default set is used. This must be a function that takes the data as input (first argument) and returns either a weights matrix or a list containing the elements \code{"graph"} for the weights matrix, \code{"intercepts"} for the intercepts (optional) and \code{"results"} for the full estimation results (optional).
}
\item{tuning}{EBIC tuning parameter, used in 'EBICglasso', 'ncvRegularize', 'IsingFit', 'huge', 'mgm', and 'ggmModSelect' default sets. Note that default values differ across methods: 'EBICglasso', 'ncvRegularize' (for convex regularization), 'huge', and 'mgm' use 0.5; 'IsingFit' uses 0.25; and 'ggmModSelect' as well as 'ncvRegularize' (for nonconvex regularization) use 0, corresponding to the basic BIC.}
\item{corMethod}{Correlation method, used in 'EBICglasso' and 'pcor' default sets. \code{"cor_auto"} uses \code{\link[qgraph]{cor_auto}} and \code{"cor_mantar"} uses \code{\link[mantar]{cor_calc}} for pearson, polychoric and polyserial correlations, \code{"cov"} uses the \code{\link{cov}} function for covariances, \code{"cor"} will use the \code{\link{cor}} function for correlations and \code{"npn"} will apply the nonparanormal transformation (via \code{huge::huge.npn}) and then compute correlations.}
\item{missing}{How to handle missing data? \code{"pairwise"} for pairwise deletion, \code{"listwise"} for listwise deletion, \code{"fiml"} for full-information maximum likelihood (called two-step EM in mantar package), \code{"stackedMI"} for stacked multiple imputation, and \code{"stop"} to stop with an error.}
\item{sampleSize}{How will sample size be computed in EBICglasso default set? The default \code{"pairwise_average"} will set the sample size to the average of sample sizes used for each individual correlation. Other options are \code{"pairwise_maximum"} (largest sample sized used for each individual correlation), \code{"pairwise_minimum"} (smallest sample sized used for each individual correlation), \code{"maximum"} (takes total number of rows including rows with NA), and \code{"minimum"} (takes total number of rows that contain no NA). The arguments "pairwise_average_v1.5", "pairwise_minimum_v1.5", and "pairwise_maximum_v1.5" can be used to mimic bootnet's behavior in version 1.5 and earlier (which also computed the sample size based on the sample sizes for the variances).}
\item{regressionSampleSize}{How will sample size be computed in nodeRegresIC default set for node-wise regression models? The default \code{"individual"} uses the number of non-missing observations for the respective dependent variable. Other options are \code{"average"} (average number of non-missing observations for all variables), \code{"max"} (maximum number of non-missing observations across all variables), and \code{"total"} (takes total number of rows including rows with NA).}
\item{corArgs}{
A list of additional arguments passed to the function specified by \code{corMethod}.
For \code{"cor_mantar"} (\code{\link[mantar]{cor_calc}}), particularly relevant arguments are \code{ordered}, which allows column-wise specification of whether variables are treated as ordered, and, when data contain missing values, \code{nimp} (number of imputations) and \code{imp_method} (imputation method used).
For \code{"cor_auto"} (\code{\link[qgraph]{cor_auto}}), important arguments include \code{detectOrdinal}, which enables automatic detection of ordinal variables, and \code{ordinalLevelMax}, which defines the maximum number of levels used to classify a variable as ordinal.
All arguments have defaults in their respective functions.
}
\item{threshold}{
Thresholding to use in partial correlation networks. Can be a fixed number to threshold all absolute edges below this value, \code{'locfdr'} for local FDR, or any option corresponding to adjustments in \code{\link[psych]{corr.p}} (\code{'none'}, \code{'sig'}, \code{'holm'}, \code{'hochberg'}, \code{'hommel'}, \code{'bonferroni'}, \code{'BH'}, \code{'BY'} or \code{'fdr'})). Can also be used for \code{default = "IsingSampler"} but can only be set to a logical enabling or disabling significance thresholding.
}
\item{refit}{Logical used in EBICglasso default set: should the estimated model structure be refitted without LASSO regularization?}
\item{rule}{The rule to use to select an edge in nodewise estimation. \code{"AND"} to only select in edge if both regression coefficients are nonzero and \code{"OR"} if only one is nonzero. Used in 'IsingFit' and 'mgm' default sets.}
\item{split}{A function or character string (\code{"median"} or \code{"mean"}) indicating how to binarize values when estimating an Ising model.}
\item{method}{The estimation method used in the IsingSampler default set (see \code{\link[IsingSampler]{EstimateIsing}}).}
\item{npn}{Logical, should nonparanormal be used in huge default set?}
\item{criterion}{The criterion used in model selection. \code{"ebic"}, \code{"ric"} or \code{"stars"} in the huge default set, \code{"bic"}, \code{"aic"} or \code{"aicc"} in the nodeRegresIC default set or \code{"EBIC"} or \code{"CV"} in the mgm default set.}
\item{nFolds}{Number of folds used in k-fold cross-validation.}
\item{type}{For mgm, see mgm or mgmfit; for relative importance networks, see \code{\link[relaimpo]{calc.relimp}}}
\item{level}{See  \code{\link[mgm]{mgm}}. Automatically set if not assigned.}
\item{order}{Order up until including which interactions are included in the model. See  \code{\link[mgm]{mgm}}. Automatically set if not assigned.}
\item{binarySign}{See  \code{\link[mgm]{mgm}}. Automatically set if not assigned.}
\item{normalized}{Should normalized relative importance be used in relative importance networks?}
\item{structureDefault}{In relative importance networks, default set used to compute the graph structure. Any other arguments used (using ...) are sent to the graph estimator function as well.}
\item{graphType}{\code{"cor"} to estimate a correlation network and \code{"pcor"} to estimate a partial correlation network (GGM)}
\item{alpha}{Significance level to test at.}
\item{principalDirection}{
Rescales variables according to the sign of the first eigen-vector. This will lead to most correlations to be positive (positive manifold), leading to negative edges to be substantively interpretable.
}
\item{stepwise}{Logical indicating if 'ggmModSelect' should use stepwise estimation.}
\item{start}{See \code{\link[qgraph]{ggmModSelect}}}
  \item{labels}{
A character vector containing the node labels. If omitted the column names of the data are used.
}
  \item{verbose}{
Logical, currently only used when \code{default = "EBICglasso"} in the \code{cor_auto} function.
}

\item{weighted}{Logical, should the analyzed network be weighted?}
\item{signed}{Logical, should the analyzed network be signed?}
\item{directed}{Logical, is the analyzed network directed? Usually does not have to be set and is detected automatically.}
\item{datatype}{
\code{"normal"} if the data argument is a data frame, or \code{"graphicalVAR"} if the data argument is a data list that can be used as input to the graphicalVAR package.}
\item{checkNumeric}{
Logical: should the data be checked to be numeric?
}
\item{lambda.min.ratio}{
Minimal lambda ratio (LASSO tuning parameter)
}
\item{nlambda}{
Number of LASSO tuning parameters to test
}
\item{nCores}{
Number of cores to use in estimating networks
}
\item{.dots}{
A list of arguments used in the estimation function set by a default set or by the \code{fun} argument.
}
\item{\dots}{
A list of arguments used in the estimation function set by a default set or by the \code{fun} argument.
}
\item{.input}{Used internally in the bootnet function. Do not use.}
\item{memorysaver}{
Logical. If TRUE attempts to save memory (RAM) by removing some objects from the output. Used by bootnet by default for bootstraps.}
\item{cutoff}{
Cutoff score for sum-score to condition on when using \code{default = "piecewiseIsing"}. This is *experimental*!
}
\item{IsingDefault}{
Default set for Ising model structure estimation in piecewise Ising estimation. This is *experimental*!
}
\item{zeroThreshold}{
Used in piecewise Ising estimation. Proportion of edges needed to be exactly 0 in pieces to set edge to zero in final network. This is *experimental*!
}
\item{minimalN}{
Used in piecewise Ising estimation. Minimal sample sizes needed in piece estimation. This is *experimental*!
}
\item{eqThreshold}{
Used in SVAR_lavaan estimation (stepup SVAR estimation). This is *experimental*! Maximum difference in criterion to decide if two models are equivalent (and select one at random).
}
\item{tempWhitelist}{
Used in SVAR_lavaan estimation (step up SVAR estimation). This is *experimental*! Matrix with edges to be whitelisted in the temporal model. The matrix should contain two columns and a row for each edge. The elements should be characters indicating the variable names of each edge (from, to).
}
\item{tempBlacklist}{
Used in SVAR_lavaan estimation (step up SVAR estimation). This is *experimental*! Matrix with edges to be blacklisted in the temporal model. The matrix should contain two columns and a row for each edge. The elements should be characters indicating the variable names of each edge (from, to).
}
\item{contWhitelist}{
Used in SVAR_lavaan estimation (step up SVAR estimation). This is *experimental*! Matrix with edges to be whitelisted in the contemporaneous model. The matrix should contain two columns and a row for each edge. The elements should be characters indicating the variable names of each edge (from, to).
}
\item{contBlacklist}{
Used in SVAR_lavaan estimation (step up SVAR estimation). This is *experimental*! Matrix with edges to be blacklisted in the contemporaneous model. The matrix should contain two columns and a row for each edge. The elements should be characters indicating the variable names of each edge (from, to).
}
\item{minimalModInd}{
Minimal modification index to consider when adding parameters in SVAR search.
}
\item{adjacency}{
  An 'adjacency' matrix indicating the graph structure (zeroes indicate a missing edge).
}
\item{nonPositiveDefinite}{
Set to \code{"stop"} to stop with an error when the input matrix is not positive definite, and to \code{"continue"} (old behavior) to continue anyway.
}
\item{unlock}{
Set to \code{TRUE} to not result in a standard error. This is to prevent using the inner functions seperatly (see details below).
}
\item{transform}{Should data be transformed before estimate the network? \code{"rank"} will call \code{\link{rank_transformation}} and \code{"quantile"}" will call \code{\link{quantile_transformation}}.}
\item{penalty}{Penalty to use in the GGMncv methods.}
\item{min_sum}{ The minimum sum score that is artifically possible in the dataset. Defaults to -Inf. Set this only if you know a lower sum score is not possible in the data, for example due to selection bias.}
}
\details{

The user can manually specify an estimation method by assigning a custom function to the \code{'fun'} argument. This function must take data as input and output an estimated network. Fhe functions \code{bootnet_} correspond to the functions used when using default sets. E.g, \code{default = "pcor"} sets \code{fun = bootnet_pcor}. As the \code{...} leads to any argument to \code{estimateNetwork} to be passed to the estimator function, the arguments described above in these estimator functions can be used for the appropriate default method. For example, if \code{default = "pcor"}, the arguments of \code{fun = bootnet_pcor} can be used in \code{estimateNetwork}. IMPORTANT NOTE: DO NOT USE THE ESTIMATOR FUNCTIONS (e.g., \code{fun = bootnet_pcor}) YOURSELF, THEY ARE ONLY INCLUDED HERE TO SHOW WHICH ARGUMENTS YOU CAN USE IN  \code{estimateNetwork}.


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
# BFI Extraversion data from psychTools package:
library("psychTools")
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
