## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(vegclust)
library(vegan)

## -----------------------------------------------------------------------------
data(wetland)
dim(wetland)

## -----------------------------------------------------------------------------
wetlandchord = decostand(wetland,"normalize")

## -----------------------------------------------------------------------------
dchord = dist(wetlandchord)

## -----------------------------------------------------------------------------
wetland.km = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KM", nstart=20)

## -----------------------------------------------------------------------------
names(wetland.km)

## -----------------------------------------------------------------------------
t(wetland.km$memb)

## -----------------------------------------------------------------------------
round(wetland.km$mobileCenters, dig=3)

## -----------------------------------------------------------------------------
wetland.kmdist = vegclustdist(x = dchord, mobileMemb=3, 
                              method="KM", nstart = 20)
names(wetland.kmdist)

## -----------------------------------------------------------------------------
wetland.kmdist$mobileCenters

## -----------------------------------------------------------------------------
t(wetland.kmdist$memb)

## -----------------------------------------------------------------------------
wetland.km$mode
wetland.kmdist$mode

## -----------------------------------------------------------------------------
round(t(wetland.km$dist2clusters), dig=2)

## -----------------------------------------------------------------------------
wetland.fcm = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=1.2, nstart=20)
round(t(wetland.fcm$memb), dig=3)

## -----------------------------------------------------------------------------
groups = defuzzify(wetland.fcm)$cluster
groups
table(groups)

## -----------------------------------------------------------------------------
groups = defuzzify(wetland.fcm, method = "cut", alpha = 0.8)$cluster
groups
table(groups, useNA = "always")

## -----------------------------------------------------------------------------
wetland.fcm2 = vegclust(x = wetlandchord, mobileCenters=3, 
                       method="FCM", m=10, nstart=20)
round(t(wetland.fcm2$memb), dig=3)

## -----------------------------------------------------------------------------
groups2 = defuzzify(wetland.fcm2, method = "cut", alpha = 0.8)$cluster
table(groups2, useNA = "always")

## -----------------------------------------------------------------------------
wetland.nc = vegclust(x = wetlandchord, mobileCenters=3,
                       method="NC", m=1.2, dnoise=0.8, nstart=20)
round(t(wetland.nc$memb), dig=2)

## -----------------------------------------------------------------------------
groups = defuzzify(wetland.nc)$cluster
groups
table(groups)

## -----------------------------------------------------------------------------
groups = defuzzify(wetland.nc, method="cut", alpha=0.8)$cluster
groups
table(groups, useNA = "always")

## -----------------------------------------------------------------------------
dist(wetland.km$mobileCenters)
dist(wetland.fcm$mobileCenters)
dist(wetland.nc$mobileCenters)

## -----------------------------------------------------------------------------
wetland.kmdd = vegclust(x = wetlandchord, mobileCenters=3, 
                      method="KMdd", nstart=20)
t(wetland.kmdd$memb)

## -----------------------------------------------------------------------------
round(wetland.kmdd$mobileCenters, dig=3)

## -----------------------------------------------------------------------------
wetland.kmdd = vegclustdist(x = dchord, mobileMemb=3, 
                      method="KMdd", nstart=20)
wetland.kmdd$mobileCenters

## -----------------------------------------------------------------------------
wetland.31 = wetlandchord[1:31,]
wetland.31 = wetland.31[,colSums(wetland.31)>0]
dim(wetland.31)
wetland.10 = wetlandchord[-(1:31),]
wetland.10 = wetland.10[,colSums(wetland.10)>0] 
dim(wetland.10)

## -----------------------------------------------------------------------------
km = kmeans(wetland.31, 2)
groups = km$cluster
groups

## -----------------------------------------------------------------------------
wetland.31.km = as.vegclust(wetland.31, groups)

## -----------------------------------------------------------------------------
wetland.31.km$method

## -----------------------------------------------------------------------------
wetland.10.km = vegclass(wetland.31.km, wetland.10)
defuzzify(wetland.10.km)$cluster

## -----------------------------------------------------------------------------
wetland.31.km.d = as.vegclust(dist(wetland.31), groups)

## -----------------------------------------------------------------------------
wetland.d.10.31 = as.data.frame(as.matrix(dchord)[32:41,1:31])

## -----------------------------------------------------------------------------
wetland.d.11.km = vegclass(wetland.31.km.d,wetland.d.10.31)
defuzzify(wetland.d.11.km)$cluster

## -----------------------------------------------------------------------------
wetland.31.nc = as.vegclust(wetland.31, groups, method="HNC", 
                            dnoise = 0.8)
wetland.10.nc = vegclass(wetland.31.nc, wetland.10)
defuzzify(wetland.10.nc)$cluster

## -----------------------------------------------------------------------------
cf = conformveg(wetland.31, wetland.10)
wetland.31.cf<- cf$x
wetland.10.cf<- cf$y
dim(wetland.31.cf)
dim(wetland.10.cf)

## -----------------------------------------------------------------------------
fixed = clustcentroid(wetland.31.cf, groups)

## -----------------------------------------------------------------------------
wetland.nc = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster

## -----------------------------------------------------------------------------
wetland.km = vegclust(wetland.10.cf, mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = "KM",
                      nstart=10)
defuzzify(wetland.km)$cluster

## -----------------------------------------------------------------------------
wetland.nc = vegclust(rbind(wetland.31.cf,wetland.10.cf), mobileCenters=1, 
                      fixedCenters = fixed, 
                      method = wetland.31.nc$method,
                      dnoise=wetland.31.nc$dnoise, nstart=10)
defuzzify(wetland.nc)$cluster

## -----------------------------------------------------------------------------
fixedDist = wetland.d.11.km$dist2clusters

## -----------------------------------------------------------------------------
wetland.km.d = vegclustdist(dist(wetland.10), mobileMemb = 1,
                            fixedDistToCenters=fixedDist, 
                            method = "KM",
                            nstart=10)
defuzzify(wetland.km.d)$cluster

## -----------------------------------------------------------------------------
fixedDist = rbind(wetland.31.km.d$dist2clusters, wetland.d.11.km$dist2clusters)

## -----------------------------------------------------------------------------
wetland.km.d = vegclustdist(dchord, mobileMemb = 1,
                            fixedDistToCenters=fixedDist, 
                            method = "KM",
                            nstart=10)
defuzzify(wetland.km.d)$cluster

## -----------------------------------------------------------------------------
groups = c(rep(1, 17), rep(2, 14), rep(3,10))

## -----------------------------------------------------------------------------
centroids = clustcentroid(wetlandchord, groups)
round(centroids, dig=3)

## -----------------------------------------------------------------------------
medoids = clustmedoid(wetlandchord, groups)
medoids

## -----------------------------------------------------------------------------
clustvar(wetlandchord, groups)

## -----------------------------------------------------------------------------
clustvar(dchord, groups)

## -----------------------------------------------------------------------------
clustvar(wetlandchord)

## -----------------------------------------------------------------------------
as.dist(as.matrix(dchord)[medoids,medoids])

## -----------------------------------------------------------------------------
dist(centroids)

## -----------------------------------------------------------------------------
interclustdist(dchord,groups)

## -----------------------------------------------------------------------------
c = clustconst(wetlandchord, memb = as.memb(groups))

## -----------------------------------------------------------------------------
d=summary(c, mode="all")

## -----------------------------------------------------------------------------
summary(c, mode="cluster", name=names(c)[1])

