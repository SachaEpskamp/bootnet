\name{ggmGenerator}
\alias{ggmGenerator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a function that simulates data from the Gaussian graphical model (GGM)
}
\description{
Generates data given a partial correlation network. Data can be made ordinal by using a threshold model with equally spaced thresholds.
}
\usage{
ggmGenerator(ordinal = FALSE, nLevels = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordinal}{
Logical, should ordinal data be generated?
}
  \item{nLevels}{
Number of levels used in ordinal data.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}