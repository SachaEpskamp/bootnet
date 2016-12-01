\name{binarize}
\alias{binarize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binarizes a dataset
}
\description{
This function will transform data into binary data (0,1). If the data is already binary, this function does nothing.
}
\usage{
binarize(x, split = "median", na.rm = TRUE, removeNArows = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data frame or matrix
}
  \item{split}{
Either a function to split on (as character or as function) or a vector. e.g., \code{split = "mean"} will split every variable on the mean of that variable, \code{split=2} will make every value above 2 a 1 and every value below 2 a 0 and a vector of the same length as each variable in the dataset will use those elements to split.
}
  \item{na.rm}{
The \code{na.rm} argument used in the split function.
}
  \item{removeNArows}{
Logical, should rows with \code{NA} be removed?
}
\item{verbose}{Output progress to the console?}
}

\value{
A binarized data frame
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
