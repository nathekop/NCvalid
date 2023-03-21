\name{Ovalid}
\alias{Ovalid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
New correlation-based cluster validity indices and other well-known cluster validity indices
}
\description{
Compute the cluster validity indices from user specified \code{kmin} to \code{kmax} used in Wiroonsri(2021) including Calinski-Harabasz(CH), Chou-Su-Lai measure(CSL), Dunn's index(DI), Davies-Bouldin's measure(DB), Davies-Bouldin* measure (DBs), Generalized Dunn 33,43,53 index(GD33,GD43,GD53), Point biserial correlation(PB), Silhouette coefficient(SC), Score function(SF), NC correlation, NCI1 and NCI2. The clustering are obtained from either K-means or hierarchical clustering.
}
\usage{
Ovalid(x,kmax,kmin=2,method='kmeans',corr='pearson',nstart=100,indexlist='all',NCstart=TRUE)
}
\arguments{
  \item{x}{the data set matrix where each column is a variable to be used for cluster analysis and each row is a data point.}
  \item{kmax}{the maximum number of clusters to be considered}
  \item{kmin}{the minimum number of clusters to be considered}
  \item{method}{a character string indicating which clustering method to be used ("kmeans", "hclust_complete", "hclust_average", or "hclust_single")
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed ("pearson", "kendall", "spearman")}
  \item{nstart}{for method="kmeans", enter the maximum number of initial random sets to be performed}
  \item{indexlist}{a character string indicating which cluster validity indices to be computed ("all", "CH", "CSL", "DB", "DBs", "DI", "PB", "SC", "SF", "NC", "NCI", "NCI1", "NCI2" or selected a set of indices as an array of any length) Note: GD33-53 are not working as they require clusterCrit which has been removed from CRAN.}
  \item{NCstart}{logical, if true, the NC correlation at k=1 is defined as an adjusted sd of the distances between all data points and their mean. Otherwise, the NC correlation at k=1 is defined as 0.}
}
\details{
The well-known cluster validity indices used in Wiroonsri(2021) include Calinski-Harabasz(CH), Chou-Su-Lai measure(CSL), Dunn's index(DI), Davies-Bouldin's measure(DB), Davies-Bouldin* measure (DBs), Generalized Dunn 33,43,53 index(GD33,GD43,GD53), Point biserial correlation(PB), Silhouette coefficient(SC) and Score function(SF).
  
The NC correlation computes the correlation between an actual distance between a pair of data points and a centroid distance of clusters that the two points locate in. NCI1 and NCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the NC improvement from \code{k-1} clusters to \code{k} clusters over the entire room for improvement. The second ratio is the NC improvement from \code{k} clusters to \code{k+1} clusters over the entire room for improvement. NCI is a combination of NCI1 and NCI2.
}
\value{
\item{NC}{the NC correlations from \code{kmin-1} to \code{kmax+1} defined in the recent paper Wiroonsri(2021)}
\item{NCI}{the NCI index from \code{kmin} to \code{kmax} defined in the recent paper Wiroonsri(2021)}
\item{NCI1}{the NCI1 index from \code{kmin} to \code{kmax} defined in the recent paper Wiroonsri(2021)}
\item{NCI2}{the NCI2 index from \code{kmin} to \code{kmax} defined in the recent paper Wiroonsri(2021)}
\item{CH}{Calinski-Harabasz index from \code{kmin} to \code{kmax}}
\item{CSL}{Chou-Su-Lai measure from \code{kmin} to \code{kmax}}
\item{DB}{Davies-Bouldin's measure from \code{kmin} to \code{kmax}}
\item{DBs}{Davies-Bouldin* measure from \code{kmin} to \code{kmax}}
\item{DI}{Dunn's index from \code{kmin} to \code{kmax}}
\item{PB}{Point biserial correlation from \code{kmin} to \code{kmax}}
\item{PBM}{PBM index from \code{kmin} to \code{kmax}}
\item{SF}{Score function from \code{kmin} to \code{kmax}}
\item{SC}{Silhouette coefficient from \code{kmin} to \code{kmax}}
\item{STR}{STR index from \code{kmin} to \code{kmax}}
}
\references{
 J. C. Bezdek, N. R. Pal, Some new indexes of cluster validity. \emph{IEEE Transactions on Systems, Man, and Cybernetics}, Part B, 28, 301-315 (1998).

G. Brock, V. Pihur, S. Datta, S. Datta, clValid: An R Package for Cluster Validation. \emph{Journal of Statistical Software}, 25(4), 1-22. (2008). \url{https://www.jstatsoft.org/v25/i04/}

 T. Calinski, J. Harabasz, A dendrite method for cluster analysis. \emph{Communications in Statistics}, 3, 1-27 (1974).

C. H. Chou, M. C. Su, E. Lai, A new cluster validity measure and its application to image compression. \emph{Pattern Anal Applic}, 7, 205-220 (2004).

D. L. Davies, D. W. Bouldin, A cluster separation measure, \emph{IEEE Trans Pattern Anal Machine Intell}, 1, 224-227 (1979).

B. Desgraupes, clusterCrit: Clustering Indices. R package version 1.2.8. (2018). \url{https://CRAN.R-project.org/package=clusterCrit}

J. C. Dunn, A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters. \emph{J Cybern}, 3(3), 32-57 (1973).

L. Kaufman, P. J. Rousseeuw, \emph{Finding Groups in Data: An Introduction to Cluster Analysis}. John Wiley & Sons (1990).

M. Kim, R. S. Ramakrishna, New indices for cluster validity assessment. \emph{Pattern Recognition Letters}, 26, 2353-2363 (2005).

M. Maechler, P. Rousseeuw, A. Struyf, M. Hubert, K. Hornik, cluster: Cluster Analysis Basics and Extensions. R package version 2.1.2. (2021).

G. W. Miligan, An examination of the effect of six types of error perturbation on fifteen clustering algorithms. \emph{Psychometrika}, 45, 325-342 (1980).
  
M. K. Pakhira, S. Bandyopadhyay and U. Maulik, Validity index for crisp and fuzzy clusters. \emph{Pattern Recogn} 37(3):487–501 (2004).

P. J. Rousseeuw, Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. \emph{Journal of Computational and Applied Mathematics}, 20, 53-65 (1987).

S. Saitta, B. Raphael, I. Smith, A bounded index for cluster validity. \emph{In Perner, P.: Machine Learning and Data Mining in Pattern Recognition, Lecture Notes in Computer Science}, 4571, Springer (2007).

A. Starczewski, A new validity index for crisp clusters. \emph{Pattern Anal Applic} 20, 687–700 (2017).   

M. Walesiak, A. Dudek, The Choice of Variable Normalization Method in Cluster Analysis, \emph{In
Soliman KS (ed.), Education Excellence and Innovation Management: A 2025 Vision to Sustain
Economic Development During Global Challenges}, 325-340. ISBN 978-0-9998551-4-1. (2020)

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

# Compute all of the indices of a kmeans result for k from 2 to 10

nc = Ovalid(x, kmax=10, nstart=100, indexlist = 'all')

# Plot the indices

ilist = c('NC', 'NCI1', 'NCI2', 'CH', 'CSL', 'DB', 'DBs', 'DI', 'GD33', 'GD43', 'GD53', 'PB', 'SF', 'SC')
par(mar = c(4, 4, 0.5, 0.5))
par(mfrow=c(3,5))
for (i in 1:14){
  plot(data.frame(nc[i]), ylab = ilist[i], xlab = "number of clusters",type='b')
}

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

# Compute some of the indices of a hierarchical clustering (average) result for k from 3 to 10

nc = Ovalid(dat, kmax=10, kmin=3, method = 'hclust_average', indexlist = c('CH','DI','NCI1','SC'))

# Plot the indices

par(mar = c(4, 4, 0.5, 0.5))
par(mfrow=c(2,2))
plot(nc$CH,ylab = "CH", xlab = "number of clusters",type='b')
plot(nc$NCI1,ylab = "NCI1", xlab = "number of clusters",type='b')
plot(nc$DI,ylab = "DI", xlab = "number of clusters",type='b')
plot(nc$SC,ylab = "SC", xlab = "number of clusters",type='b')
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
