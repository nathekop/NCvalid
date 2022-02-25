# NCvalid
New correlation based cluster validity indices

## Description

NCvalid is a package used for analysing crisp clustering results such as k-means and hierarchical clustering. 
It contains a main function (NCvalid) that computes NC correlation, NCI1 and NCI2 cluster validity indices in a user specified range of the number of clusters.
The correlation-based indices are recently introduced in Wiroonsri(2021).  Our indices constantly yield several peaks at different numbers of clusters 
which allow the users to rank their options. This benefits the users who work in the area that the number of clusters is not necessarily fixed. 
Furthermore, a global peak can accurately detect the number of cluster in the case that it is fixed. 
The package also provides a function (Ovalid) that computes other well-known validity incides used in Wiroonsri(2021).

## Installation

The lastest development version of the package can be loaded directly from GitHub using the devtools package:

```bash
install.packages("devtools")
devtools::install_github("nwiroonsri/NCvalid")
```

## Simple Example

```r
#Evaluate a k-means result on iris dataset

x <- iris
x <- scale(x[,1:4])

#Compute the indices of a kmeans result for k from 2 to 10

nc = NCvalid(x, kmax=10)

#Plot NC correlation, NCI1, NCI2

par(mar = c(4, 4, 0.5, 0.5))
par(mfrow=c(1,3))
plot(nc$NC,ylab = "NC", xlab = "number of clusters",type='b')
plot(nc$NCI1,ylab = "NCI1", xlab = "number of clusters",type='b')
plot(nc$NCI2,ylab = "NCI2", xlab = "number of clusters",type='b')
```
## More Example

```r
#Evaluate a hierarchical clustering (complete linkage) result on a simulated dataset
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
graphics.off()

# Compute 6 cluster validity indices of a hierarchical clustering (average) result for k from 2 to 10

ilist =  c('NCI1','NCI2','CH','DI','SC','GD33')
nc = Ovalid(dat, kmax=10, kmin=2, method = 'hclust_average', indexlist = ilist)

# Plot the indices in the list

par(mar = c(4, 4, 0.5, 0.5))
par(mfrow=c(3,2))
for (i in 1:6){
  plot(data.frame(nc[i]), ylab = ilist[i], xlab = "number of clusters",type='b')
}
```

## References

N. Wiroonsri (2021). Clustering performance analysis using new correlation based cluster validity indices, [arXiv:2109.11172](https://arxiv.org/abs/2109.11172).

## License

The NCvalid package as a whole is distributed under [GPL(>=3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
