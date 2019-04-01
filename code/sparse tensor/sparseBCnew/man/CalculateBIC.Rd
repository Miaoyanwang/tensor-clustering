\name{CalculateBIC}
\alias{CalculateBIC}
\title{ BIC criterion
}
\description{
Calculate BIC
}
\usage{
CalculateBIC(x, bires, method="L0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix; samples are rows and columns are features.  Cannot contain missing values.
}
  \item{bires}{
The number of row clusters, i.e., the number of clusters for the observations.
}
  \item{method}{
  Two options: "L0" and "L1" which refers to the corresponding penalty.
  }
}
\value{
\item{BIC}{BIC value}
}