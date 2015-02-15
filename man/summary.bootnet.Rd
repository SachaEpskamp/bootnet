\name{summary.bootnet}
\alias{summary.bootnet}
\title{
Summarize bootnet results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{bootnet}(x, types = c("edge", "intercept", "strength", "closeness", "betweenness", "distance"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{bootnet} object
}
  \item{types}{
The types of statistics to include in the summary table
}
}
\value{
A \code{tbl_df} (data frame) containing summarized statistics.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{

}
