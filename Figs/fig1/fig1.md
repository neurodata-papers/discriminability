# discrimability
Greg Kiar  
January 13, 2017  

# Discriminability in a Nutshell

We are going to generate data from two classes (circles and squares), and demonstrate the
strength of discriminability through evaluating different sampling or processing strategies.
We construct our two class distributions, sample points from them, reconstruct the object,
compute the discriminability of the processing. Discriminability is: computing distance between
all pairs of objects, rank the distances, compute the histogram of same-class ranks, and compute
the mean rank of the distribution.

Before getting started, we will define a few helpful functions and variables which will streamline
this process.


```r
library(abind)
library(rmarkdown)
library(knitr)
library(RColorBrewer)
library(gplots)
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```




```r
lweight <- 4
len <- 500
nsamples <- 40
nsims <- 3
noise1 <- 0.15
tsize <- 4
darkcol <- '#28282e'
lightcol <- '#98929e'

circle <- function(){
  theta <- seq(from = 0, to = 2*pi - pi/len, by = 2*pi/len)
  x <- cos(theta)
  y <- sin(theta)
  circ <- list(x, y)
}

square <- function(){
  range <- seq(from = -1, to = 1, by = 2/(len/4))
  static <- rep(1, each = len/4)
  x <- c(range, static, -range, -static) 
  y <- c(static, -range, -static, range)
  sqr <- list(x, y)
}

shape_plot <- function(xs, ys, typ){
  size <- lweight
  if (typ == 'p'){
    size <- size/3
  }
  plot(xs, ys, type=typ, axes=FALSE, xlab='', ylab='', pch=20, asp=1, lwd=size, col=darkcol)
  if (typ == 'l'){
    polygon(x=xs, y=ys, col=lightcol)
  }
}

sample1 <- function(data, n){
  pts <- sort(floor(runif(n, min=1, max=len)))
  xs <- data[[1]]
  x <- xs[pts] + noise1*runif(n, min=-1, max=1)
  x <- c(x, x[1])
  ys <- data[[2]]
  y <- ys[pts] + noise1*runif(n, min=-1, max=1)
  y <- c(y, y[1])
  samp <- list(x, y)
}

distances <- function(data){
  diff <- matrix(data=NA, nrow=2*nsims, ncol=2*nsims)
  for (i in 1:(2*nsims)){
    for (j in 1:(2*nsims)){
      diff[i,j] <- norm(data[,,i] - data[,,j], '2')
    }
  }
  cols <- colorRampPalette(c("green", "yellow", "red"))(n = 299)
  heatmap.2(diff, cellnote=round(diff,1), Colv=NA, Rowv=NA, trace="none",
            density.info="none", notecol="black", col=cols, cexRow=tsize,
            cexCol=tsize, notecex=tsize, labRow="", labCol="", key=FALSE,
            lwid=c(.5,6), lhei=c(.5,6), margins=c(4.5,4.5))

}

simulate_shape <- function(fn){
  layout(matrix(c(1,2,3, 1,4,5, 1,6,7), 3, 3, byrow = TRUE))
  par(xpd=NA)
  shape <- fn()
  shape_plot(shape[[1]], shape[[2]], 'l')
  
  arrows(x0=1.2, x1=2.2, y0=1, y1=1.5, length=0.1, lwd = lweight)
  arrows(x0=1.3, x1=2.2, y0=0, y1=0, length=0.1, lwd = lweight)
  arrows(x0=1.2, x1=2.2, y0=-1, y1=-1.5, length=0.1, lwd = lweight)
  text(x=1.5, y=1.5, expression('s'[1]^'a'), cex=tsize)
  text(x=1.7, y=0.3, expression('s'[2]^'a'), cex=tsize)
  text(x=1.5, y=-1.5, expression('s'[3]^'a'), cex=tsize)

  arrows(x0=3.5, x1=4.8, y0=1.5, y1=1.5, length=0.1, lwd=lweight)
  arrows(x0=3.5, x1=4.8, y0=0, y1=0, length=0.1, lwd=lweight)
  arrows(x0=3.5, x1=4.8, y0=-1.5, y1=-1.5, length=0.1, lwd=lweight)

  text(x=4.15, y=1.8, bquote(p[.(1)]^'a'), cex=tsize)
  text(x=4.15, y=0.3, bquote(p[.(2)]^'a'), cex=tsize)
  text(x=4.15, y=-1.2, bquote(p[.(3)]^'a'), cex=tsize)
  
  # label <- rep(as.character(substitute(fn)), nsims)
  label <- c()
  data <- array(data=NA, c(2, nsamples, nsims))
  dim(data)
  for (i in 1:nsims){
    samp <- sample1(shape, nsamples)
    shape_plot(samp[[1]], samp[[2]], 'p')
    shape_plot(samp[[1]], samp[[2]], 'l')
    data[1, ,i] <- samp[[1]][1:nsamples]
    data[2, ,i] <- samp[[2]][1:nsamples]
    label <- c(label, as.character(substitute(fn)))
  }
  values <- list(label, data)
}
```

## Sampling shapes

We will start by sampling data from the circle distribution and reconstructing the shapes
from each set of samples.


```r
# pdf("circle_plot.pdf", width=20, height=10)
vals <- simulate_shape(circle)
```

<figure><img src="./Figures/circle_sampling-1.png"><figcaption></figcaption></figure>

```r
# dev.off()
labels <- vals[[1]]
data <- vals[[2]]
```

The same process is then repeated for the square distribution.

```r
vals <- simulate_shape(square)
```

<figure><img src="./Figures/square_sampling-1.png"><figcaption></figcaption></figure>

```r
labels <- c(labels, vals[[1]])
data <- abind(data, vals[[2]])
```

Now we have a data matrix which contains `x` and `y` sampled values, and a list of labels
for each sample identifying whether it belongs to the "circle" or "square" distribution.
We can then compute the distance between observations in each class.

```r
print(labels)
```

```
# [1] "circle" "circle" "circle" "square" "square" "square"
```

```r
print(dim(data))
```

```
# [1]  2 40  6
```

```r
diff <- distances(data)
```

<figure><img src="./Figures/heatmap-1.png"><figcaption></figcaption></figure>