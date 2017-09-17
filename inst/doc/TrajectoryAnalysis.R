### R code from vignette source 'TrajectoryAnalysis.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: TrajectoryAnalysis.Rnw:18-20
###################################################
options(width=67)
library(vegclust)


###################################################
### code chunk number 2: TrajectoryAnalysis.Rnw:24-25
###################################################
sites = c(1,1,2,2)


###################################################
### code chunk number 3: TrajectoryAnalysis.Rnw:28-32
###################################################
xy1<-rbind(c(0,0),c(0,2.5),
           c(1,0.5),c(1,3))
segmentDistances(dist(xy1),sites,distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy1),sites,distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 4: TrajectoryAnalysis.Rnw:36-40
###################################################
xy2<-rbind(c(0,0),c(0,2.5),
           c(0,1),c(1.8,2.5))
segmentDistances(dist(xy2),sites, distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy2),sites, distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 5: TrajectoryAnalysis.Rnw:44-48
###################################################
xy3<-rbind(c(0,0),c(0,2.5),
           c(-0.5,2),c(2,2))
segmentDistances(dist(xy3),sites,distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy3),sites,distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 6: TrajectoryAnalysis.Rnw:51-55
###################################################
xy4<-rbind(c(0,0),c(0,2.5),
           c(0,2.5),c(1.8,1))
segmentDistances(dist(xy4),sites, distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy4),sites, distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 7: TrajectoryAnalysis.Rnw:58-62
###################################################
xy5<-rbind(c(0,0),c(0,2.5),
           c(1,3),c(1,0.5))
segmentDistances(dist(xy5),sites, distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy5),sites, distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 8: TrajectoryAnalysis.Rnw:65-69
###################################################
xy6<-rbind(c(0,2.5),c(0,0),
           c(1,0.5),c(1,3))
segmentDistances(dist(xy6),sites, distance.type = "Hausdorff")$Dseg
segmentDistances(dist(xy6),sites, distance.type = "directed-segment")$Dseg


###################################################
### code chunk number 9: TrajectoryAnalysis.Rnw:74-75
###################################################
sites = c(1,1,1,2,2,2)


###################################################
### code chunk number 10: TrajectoryAnalysis.Rnw:78-86
###################################################
xy1<-matrix(0, nrow=6, ncol=2)
xy1[2,2]<-1
xy1[3,2]<-2
xy1[4:6,1] <- 0.5
xy1[4:6,2] <- xy1[1:3,2]
trajectoryDistances(dist(xy1),sites, distance.type = "Hausdorff")
trajectoryDistances(dist(xy1),sites, distance.type = "SPD")
trajectoryDistances(dist(xy1),sites, distance.type = "DSPD")


###################################################
### code chunk number 11: TrajectoryAnalysis.Rnw:89-94
###################################################
xy2<-xy1
xy2[6,]<-c(1,1.8)
trajectoryDistances(dist(xy2),sites, distance.type = "Hausdorff")
trajectoryDistances(dist(xy2),sites, distance.type = "SPD")
trajectoryDistances(dist(xy2),sites, distance.type = "DSPD")


###################################################
### code chunk number 12: TrajectoryAnalysis.Rnw:97-102
###################################################
xy3<-xy2
xy3[4,]<-c(1,0.2)
trajectoryDistances(dist(xy3),sites, distance.type = "Hausdorff")
trajectoryDistances(dist(xy3),sites, distance.type = "SPD")
trajectoryDistances(dist(xy3),sites, distance.type = "DSPD")


###################################################
### code chunk number 13: TrajectoryAnalysis.Rnw:105-111
###################################################
xy4<-xy2
xy4[4,]<-xy2[6,]
xy4[6,]<-xy2[4,]
trajectoryDistances(dist(xy4),sites, distance.type = "Hausdorff")
trajectoryDistances(dist(xy4),sites, distance.type = "SPD")
trajectoryDistances(dist(xy4),sites, distance.type = "DSPD")


