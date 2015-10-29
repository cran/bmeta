\name{funnel.plot}
\alias{funnel.plot}
\title{
Funnel plot to examine publication bias
}
\description{
Function to examine publication bias. For fixed-effects models, estimates from no-pooling
effects model are used as study-specific estimates and for random-effects models, 
study-specific estimates are obtained from the posterior distribution.   
}
\usage{
funnel.plot(x,xlab=NULL,ylab=NULL,title=NULL)
}
\arguments{
  \item{x}{
a \code{bmeta} object with results of the model
}
  \item{xlab}{
title of x-axis. If unspecified, the function sets an appropriate lable by default.
}
  \item{ylab}{
title of x-axis. If unspecified, the function sets an appropriate lable by default.
}
  \item{title}{
title of the plot if specified
}
}
\author{
Tao Ding
Gianluca Baio
}
\examples{
### Read and format the data (binary)
data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))

### List data for binary outcome 
data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 

### Select random-effects meta-analysis with t-distribution prior for binary
### data
x <- bmeta(data.list, outcome="bin", model="std.dt", type="ran")

### using output from bmeta to produce funnel plot 
funnel.plot(x)

### using output from bmeta and specify title of the plot 
funnel.plot(x,title="funnel plot")
}
\keyword{Funnel plot}
