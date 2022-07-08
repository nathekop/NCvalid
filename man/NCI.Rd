\name{NCI}
\alias{NCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
New correlation-based cluster validity index
}
\description{
Compute the NCI cluster validity index for the number of clusters from user specified \code{kmin} to \code{kmax} obtained from either K-means or hierarchical clustering based on the recent paper by Wiroonsri(2021).
}
\usage{
NCI(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100,NCstart=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the data set matrix where each column is a variable to be used for cluster analysis and each row is a data point.}
  \item{kmax}{the maximum number of clusters to be considered}
  \item{kmin}{the minimum number of clusters to be considered}
  \item{method}{a character string indicating which clustering method to be used ("kmeans", "hclust_complete", "hclust_average", or "hclust_single")
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed ("pearson", "kendall", "spearman")}
  \item{nstart}{for method="kmeans", enter the maximum number of initial random sets to be performed}
  \item{NCstart}{logical, if true, the NC correlation at k=1 is defined as an adjusted sd of the distances between all data points and their mean. Otherwise, the NC correlation at k=1 is defined as 0.}
}
\details{
The NC correlation computes the correlation between an actual distance between a pair of data points and a centroid distance of clusters that the two points locate in. NCI1 and NCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the NC improvement from \code{k-1} clusters to \code{k} clusters over the entire room for improvement. The second ratio is the NC improvement from \code{k} clusters to \code{k+1} clusters over the entire room for improvement. NCI is a combination of NCI1 and NCI2.
}
\value{
\item{NCI}{the NCI index from \code{kmin} to \code{kmax} defined in the recent paper Wiroonsri(2021)}
}
\references{
N. Wiroonsri. Clustering performance analysis using new correlation based cluster validity indices. \url{https://arxiv.org/abs/2109.11172}, 2021.
}
\author{
Nathakhun Wiroonsri <\email{nathakhun.wir@kmutt.ac.th}> \url{https://sites.google.com/view/nwiroonsri/}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## kmeans on iris dataset

x <- iris
x <- scale(x[,1:4])

# Compute NCI index of a kmeans result for k from 2 to 10

nc = NCI(x, kmax=10, nstart=100)

# Plot NCI

plot(nc,ylab = "NCI", xlab = "number of clusters",type='b')

##----------------------------------------------------------------

## Hierarchical clustering on generated dataset

a = runif(550,min = 0,max = 2*pi)
r = runif(550,min = 0, max = 1)

x = r*cos(a)
y = r*sin(a)

dat = cbind(x,y)

dat[1:150,1] = dat[1:150,1]-5
dat[151:300,1] = dat[151:300,1]+5
dat[301:350,1] = dat[301:350,1]-3
dat[351:400,1] = dat[351:400,1]+3

plot(dat)

# Compute NCI index of a hierarchical clustering (complete) result for k from 3 to 10

nc = NCI(dat, kmax=10, kmin=3, method = 'hclust_complete')

# Plot NCI

plot(nc, ylab = "NCI", xlab = "number of clusters",type='b')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
Footer
© 2022 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
Footer
© 2022 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
