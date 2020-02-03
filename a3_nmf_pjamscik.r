
source("util.R")

# load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"

# the data is sorted by newsgroup; this is not exploited/known to NMF
# throghout the notebook I truncate values >0.00001 for improved visiblity
# I also decided not to plot the legend
ggplotm(pmin(P,0.000001), format="", show.axis=FALSE, mid="black",show.legend=FALSE)

# values for box boundaries are are derived in understandin_sorting_order.ipynb
P[1:200+1,383]<-1
P[100:300,609]<-1
P[200:400,718]<-1
P[100,1:609]<-1
P[200, 383:718]<-1
P[300, 609:800]<-1
ggplotm(pmin(P,0.000001), format="", show.axis=FALSE, mid="black",show.legend=FALSE)

############
## 1. NMF ##
############

## load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"

## run NMF using GKL
r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)

## also compute SVD
P.svd <- svd(P)
rownames(P.svd$v) <- colnames(P)


## print results - NMF
# non-deterministic: row order may change
# also top-10 values may change 
with(lr.gkl,
     for (k in 1:nrow(R)) {
         print(rev(sort(R[k,]))[1:10])
         cat(strrep('-',80), "\n")
     })



# 1c) compare top-10 terms of SVD and NMF
# SVD is deterministic: we will always obtain the same row/topic order and the same terms
# SVD s hard to interpret: describing objects by absence of attributes is counterintuitive for humans
# by looking at high values > 0 we can (arguably) derive row 3 = crypt (but note that 'encrypt' also among top-10 if row 4); row 2 = med

with(P.svd,
     for (k in 1:r) {
         y=order(abs(v[,k]))
         print(rev(v[y,k])[1:10])
         cat(strrep('-',80), "\n")
     })

## look at the reconstruction -  NMF

## load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"

## run NMF using GKL
r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)


Phat <- lr.gkl$L %*% lr.gkl$R
summary(as.vector(Phat))

# important to do this after the decomposition
Phat[0:200,383+1]<-1
Phat[100:300,609+1]<-1
Phat[200:400,718+1]<-1
Phat[100,1:609]<-1
Phat[200, 383:718]<-1
Phat[300, 609:800]<-1

ggplotm(pmin(Phat, 0.00001), format="", show.axis=FALSE, mid="black", show.legend=FALSE)


## look at the decomposition (here: Rtilde and V)
## NMF - Rtilde 
#ggplotm(lr.gkl$R, format="", show.axis=FALSE, mid="black") # enlarge plot or subselect columns to see something
ggplotm(lr.gkl$R[,40:60], format="", rotate.labels=TRUE, mid="black") # first 30 columns

# find out ground truth groups by manually looping through the data;
# alternative: reverse-engineering of sorting logic (see helper.ipynb that contains corresponding python file)
# first block: cipher : 250
# second block: medicin,  e.g. labels 250:450
# second last block: space, e.g. labels 600:650
# last block: religion, e.g. labels: 750:800
ggplotm(t(P.svd$v[250:300,1:4]), format="", rotate.labels=TRUE, mid="black") # first 30 columns


## recall the original data matrix
## the data is sorted by newsgroup; this is not exploited/known to NMF
## here we truncate values >0.000001 for improved visiblity
P[0:200,383+1]<-1
P[100:300,609+1]<-1
P[200:400,718+1]<-1
P[100,1:609]<-1
P[200, 383:718]<-1
P[300, 609:800]<-1
ggplotm(pmin(P,0.0001), format="", show.axis=FALSE, mid="black")


# look at the reconstruction - SVD
# observe that e.g. for topic 2 we reconstruct negative values, which is clearly a mistake
# one should never observe negative values
Phat <- P.svd$u[,1:4] %*% diag(P.svd$d[1:4]) %*% t(P.svd$v[,1:4])
Phat[0:200+1,383+1]<-1
Phat[100:300+1,609+1]<-1
Phat[200:400,718]<-1
Phat[100,1:609+1]<-1
Phat[200, 383:718]<-1
Phat[300, 609:800]<-1
#summary(as.vector(Phat))
ggplotm(pmax(pmin(Phat,0.00001),-0.00001), format="", show.axis=FALSE, mid="black", show.legend=FALSE)

# SVD - V
# have a look at the column representing drug, good illustration
#ggplotm(t(P.svd$v[,1:4]), format="", show.axis=FALSE, mid="black")
ggplotm(t(P.svd$v[50:60,1:4]), format="", rotate.labels=TRUE, mid="black") # first 30 columns

## run NMF using GKL with rank = 2
r <- 2
lr.gkl.r2 <- lee01.gkl(P, r, reps=5)

## print results
with(lr.gkl.r2,
     for (k in 1:nrow(R)) {
         print(rev(sort(R[k,]))[1:10])
         cat(strrep('-',80), "\n")
     })

## look at the reconstruction
Phat <- lr.gkl.r2$L %*% lr.gkl.r2$R
Phat[0:200+1,383+1]<- 1
Phat[100:300+1,609+1]<-  1
Phat[200:400,718]<- 1
Phat[100,1:609+1]<- 1
Phat[200, 383:718]<- 1
Phat[300, 609:800]<- 1
#summary(as.vector(Phat.r2))
ggplotm(pmin(Phat, 0.00001), format="", show.axis=FALSE, mid="black", show.legend=FALSE)

## run NMF using GKL with rank = 8
r <- 8
lr.gkl.r8 <- lee01.gkl(P, r, reps=5)

## print results
with(lr.gkl.r8,
     for (k in 1:nrow(R)) {
         print(rev(sort(R[k,]))[1:10])
         cat(strrep('-',80), "\n")
     })

# look at the reconstruction
Phat <- lr.gkl.r8$L %*% lr.gkl.r8$R
Phat[0:200+1,383+1]<- 1
Phat[100:300+1,609+1]<-  1
Phat[200:400,718]<- 1
Phat[100,1:609+1]<- 1
Phat[200, 383:718]<- 1
Phat[300, 609:800]<- 1
#summary(as.vector(Phat.r8))
ggplotm(pmin(Phat, 0.00001), format="", show.axis=FALSE, mid="black", show.legend = FALSE)

## Gaussian NMF
# Note: if we run this code cell several times, we may obtain different factorizations
# n
r <- 4
lr.gnmf <- lee01.gnmf(P, r, reps=5)

## YOUR CODE HERE
## print results
with(lr.gnmf,
     for (k in 1:nrow(R)) {
         print(rev(sort(R[k,]))[1:10])
         cat(strrep('-',80), "\n")
     })

## look at the reconstruction
Phat <- lr.gnmf$L %*% lr.gnmf$R
Phat[0:200+1,383+1]<- 1
Phat[100:300+1,609+1]<-  1
Phat[200:400,718]<- 1
Phat[100,1:609+1]<- 1
Phat[200, 383:718]<- 1
Phat[300, 609:800]<- 1
#summary(as.vector(Phat.gnmf))
ggplotm(pmin(Phat.gnmf, 0.00001), format="", show.axis=FALSE, mid="black")

#############
## 2. PLSA ##
#############

## computing the 3-matrix decompositions
r <- 4
lr.gkl <- lee01.gkl(P, r, reps=10)
lsr.gkl <- nmf.lsr(lr.gkl) ## result as: L %*% S %*% R#
slr.gkl <- nmf.slr(lr.gkl) ## result as: S %*% L %*% R

## YOUR CODE HERE

# Explore S
# S is not normalized
lsr.gkl$S
sum(lsr.gkl$S)
lsr.gkl$S/ sum(lsr.gkl$S)

# Explore S
# S is not normalized
lsr.gkl$S
sum(lsr.gkl$S)
lsr.gkl$S/ sum(lsr.gkl$S)

# Explore R
# Rows of R sum to one
sum(lsr.gkl$R[1,])
sum(lsr.gkl$R[2,])

with(lsr.gkl,
     for (k in 1:nrow(R)) {
         print(rev(sort(R[k,]))[1:10])
         cat(strrep('-',80), "\n")
     })

# Explore L
# columns of L sum to one
#lsr.gkl$L
dim(lsr.gkl$L)
sum(lsr.gkl$L[,1])
sum(lsr.gkl$L[,4]) 


# Explore L
# Rows sum to one
dim(slr.gkl$L)
sum(slr.gkl$L[1,])
sum(slr.gkl$L[400,])






# Explore R
# Rows sum to one
dim(slr.gkl$R)
sum(slr.gkl$R[3,])

 


# Explore S
with(slr.gkl,{
    O <- S %*% L %*% R
    S[ S < 0.0007045] <- 0
    C <- S %*% L %*% R
    print(norm(C - O, "F"))
    print(norm(P - O, "F"))
    print(norm(P - C, "F"))
})
sum(slr.gkl$S)
summary(slr.gkl$S)
slr.gkl$S
slr.gkl$S/sum(slr.gkl$S)

## true labels (DO NOT USE for clustering)
## 1=sci.crypt
## 2=sci.med
## 3=sci.space
## 4=soc.religion.christian
labels <- rep(1:4, each=100)

cluster <- kmeans(P, 4, nstart=100)$cluster

## to compute the confusion matrix between a clustering and the true labels, we
## first relabel every cluster so that cluster ids match labels to the extent
## possible. (Always do this.)
cm <- function(cluster.ids) {
    cluster.ids.matched <- match.labels(cluster.ids, labels)
    u = union(cluster.ids.matched, labels)
    t = table(factor(cluster.ids.matched, u), factor(labels, u))
    confusionMatrix(t)
}

## example clustering (k-means with k=4)
cluster <- kmeans(P, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"] 
cm(cluster)               

P.svd <- svd(P)
u.4 <- P.svd$u[,1:4]
d.4 <- diag(P.svd$d[1:4])

cluster <- kmeans(u.4 %*% d.4, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)
cluster <- kmeans(lr.gkl$L, 4, nstart=100)$cluster

cm(cluster)$overall["Accuracy"]
cm(cluster)

r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)
lsr.gkl <- nmf.lsr(lr.gkl) 

cluster <- kmeans(lsr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)
slr.gkl <- nmf.slr(lr.gkl) 

cluster <- kmeans(slr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# let's look into the rows of L again 
# and observe that probably the high confidence is responsible for the high accuracy
slr.gkl$L

# Cluster on L of SLR
# Acc.: 95.75% (non-deterministic)
source("util.R")

# load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news_filtered.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"


r <- 4
lr.gkl <- lee01.gkl((P>0), r, reps=5) # simply change P to (P>0)
slr.gkl <- nmf.slr(lr.gkl) 

cluster <- kmeans(slr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# Cluster on W of WH
# Acc: 42.75%; for class 3/4 we improve accuracy
r <- 4
lr.gkl <- lee01.gkl((P>0), r, reps=5)
cluster <- kmeans(lr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# cluster on U of SVD
# again, we impprove only on topic 4
P.svd <- svd((P>0))
u.4 <- P.svd$u[,1:4]
d.4 <- diag(P.svd$d[1:4])
cluster <- kmeans(u.4 %*% d.4, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# Cluster on L of SLR
# Acc 45%: Disappointing

source("util.R")

# load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news_1.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"

# the data is sorted by newsgroup; this is not exploited/known to NMF
# throghout the notebook I truncate values >0.00001 for improved visiblity
# I also decided not to plot the legend
ggplotm(pmin(P,0.000001), format="", show.axis=FALSE, mid="black",show.legend=FALSE)


r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)
slr.gkl <- nmf.slr(lr.gkl) 

cluster <- kmeans(slr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# Cluster on W of WH
# Acc:33%; for class 4 we improve accuracy
r <- 4
lr.gkl <- lee01.gkl(P, r, reps=5)
cluster <- kmeans(lr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# cluster on U of SVD
P.svd <- svd(P)
u.4 <- P.svd$u[,1:4]
d.4 <- diag(P.svd$d[1:4])
cluster <- kmeans(u.4 %*% d.4, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# Cluster on L of SLR
# Acc 44.75%: Disappointing

source("util.R")

# load the document-term matrix
Ptilde <- as.matrix( read.csv("data/news_1.csv") )
P <- Ptilde/sum(Ptilde) ## replace word frequencies by "probabilities"

r <- 4
lr.gkl <- lee01.gkl((P>0), r, reps=5)
slr.gkl <- nmf.slr(lr.gkl) 

cluster <- kmeans(slr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# Cluster on W of WH
# Acc: 47%; for class 3/4 we improve accuracy
r <- 4
lr.gkl <- lee01.gkl((P>0), r, reps=5)
cluster <- kmeans(lr.gkl$L, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)

# cluster on U of SVD
P.svd <- svd((P>0))
u.4 <- P.svd$u[,1:4]
d.4 <- diag(P.svd$d[1:4])
cluster <- kmeans(u.4 %*% d.4, 4, nstart=100)$cluster
cm(cluster)$overall["Accuracy"]
cm(cluster)
