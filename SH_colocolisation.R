rm(list=ls())
#library(compiler)
library(EBImage)
library("gtools")
library(outliers)

setwd("/scratch/SH_colocolisation/second/20200207_152252_921/")
tif<-dir()[grep(".tiff",dir())]
frame<-readImage("WellC06_PointC06_0001_ChannelCSU 640,CSU 561,CSU 488,CSU 405_Seq0091.tiff",type = "tif",all = T)
#frame<-readImage(tif[32],type = "tif",all = T)
dim.y<-dim(frame)[1]
dim.X<-dim(frame)[2]

UB<-frame[1:dim.y,1:dim.X,1]
tax<-frame[1:dim.y,1:dim.X,2]
GFP<-frame[1:dim.y,1:dim.X,3]
dapi<-frame[1:dim.y,1:dim.X,4]
colorMode(UB)<-"Grayscale"
colorMode(tax)<-"Grayscale"
colorMode(GFP)<-"Grayscale"
colorMode(dapi)<-"Grayscale"
display(UB*20,"raster")
display(tax*20,"raster")



dapi_normal<- dapi*5
display(dapi_normal,"raster")
nmask2 = thresh(dapi_normal, 85, 85,0.003)
mk3 = makeBrush(13, shape= "diamond")
nmask2 = opening(nmask2, mk3) 
nmask2 = fillHull(nmask2) 
nseg = bwlabel(nmask2)  #binery to object
seg_pink = paintObjects(nseg,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
display(nseg,"raster")
rm(nmask2)
chackpoint<-computeFeatures.shape(nseg)
nmask<-nseg
if (length(chackpoint) < 3 ){
  next()
}
if (is.null(nrow(chackpoint))){
  next()
}
if (nrow(chackpoint) > 100){
  next()
}
nf = computeFeatures.shape(nmask)
nr = which(nf[,2] < 50) 
nseg = rmObjects(nmask, nr) 
nn = max(nseg) 
colorMode(nseg)<-Grayscale
chackpoint<-computeFeatures.shape(nseg)
if (length(chackpoint) < 3 ){
  next()
}
if (is.null(nrow(chackpoint))){
  next()
}
if (nrow(chackpoint) > 2000 ){
  next()
}
int.dapi<-computeFeatures.basic(nseg,dapi)
y<-which(scores(int.dapi[,1], type="z",prob = 0.975))
tp<-as.numeric(attr(y,"names"))
if (length(tp) < 1 ){
  next()
}
nseg<-rmObjects(nseg,tp)
chackpoint<-computeFeatures.shape(nseg)
df<-as.data.frame(chackpoint)
rm(chackpoint,nmask)
xy<-computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
if (length(xy) < 3 ){
  next()
}
if (is.null(nrow(xy))){
  next()
}
if (nrow(xy) > 500 ){
  next()
}
gsegg=nseg
display(colorLabels(gsegg),"raster")

cell_normal<- GFP*10
display(cell_normal,"raster")
smooth<-makeBrush(19,shape = "disc") 
cell_normal<-filter2(cell_normal,smooth, boundary = c("circular", "replicate"))
thr_cell<-thresh(cell_normal, w=10, h=10, offset=0.1)
colorMode(thr_cell)<-Grayscale
cmask = opening(thr_cell, kern=makeBrush(7,shape="disc"))
rm(thr_cell)
ar<-as.vector(cell_normal)
ar.sum<-as.numeric(summary(ar))
open2<-opening(cell_normal>(ar.sum[2]))
rm(ar)
colorMode(open2)<-Grayscale
open_2 = paintObjects(open2,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
combine<-cmask
combine[open2 > cmask]<-open2[open2 > cmask]
#display(combine,method="raster")
combine[gsegg > combine]<-gsegg[gsegg > combine]
#display(colorLabels(combine),method="raster")
#display(colorLabels(gsegg),method="raster")

open_2 = paintObjects(combine,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)

csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=cmask)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
#display(colorLabels(csegpink),method="raster")

csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=combine)
#rm(open_2,cmask)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
#display(colorLabels(csegpink),method="raster")
xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
if (length(xy) < 3 ){
  next()
}
if (is.null(nrow(xy))){
  next()
}
if (nrow(xy) > 2000 ){
  next()
}
cf = computeFeatures.shape(csegpink)
cf = computeFeatures.shape(csegpink)
#cr = which(cf[,2] < 200) 
#csegpink = rmObjects(csegpink, cr,reenumerate = F) 
cf = computeFeatures.shape(csegpink)
cfErea<-data.frame(cf[,1])
cfErea$num<-row.names(cfErea)
ci = which(cf[,1] > 35000) 
csegpink = rmObjects(csegpink, ci,reenumerate = F) 
rm(cf,ci,cfErea)
#' select all the cells with normal ilumination propteis
int.GFP<-computeFeatures.basic(csegpink,GFP)
y<-which(scores(int.GFP[,1], type="z",prob = 0.80))
# seg_pink = paintObjects(csegpink,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
# display(seg_pink,method="raster")
ty<-as.numeric(attr(y,"names"))
csegpink<-rmObjects(csegpink,ty,reenumerate = F)
# seg_pink = paintObjects(test,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
# display(seg_pink,method="raster")




xy.gsegg<-as.numeric(row.names(computeFeatures.moment(gsegg)[,c('m.cx','m.cy')]))
length(xy.gsegg)
xy.cseg<- as.numeric(row.names(computeFeatures.moment(csegpink)[,c('m.cx','m.cy')])) 
length(xy.gsegg)
ind.deff<-setdiff(xy.gsegg,xy.cseg)
gsegg<-rmObjects(gsegg,ind.deff,reenumerate=F)
gsegg<-reenumerate(gsegg)
csegpink<-reenumerate(csegpink)
seg_pink = paintObjects(csegpink,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
display(colorLabels(csegpink),method="raster")



stack<-stackObjects(csegpink,GFP*10,ext = c(200,200))
stack.UB<-stackObjects(csegpink,UB*10,ext = c(200,200))

ind<-NULL
for (m in 1:dim(stack)[3]){
  st.temp<-stack[,,m]
  bg<-length(which(st.temp==0))
  fg<-length(which(st.temp > 0))
  if ((bg/fg) > 10){
    next()
  } else {
    ind<-c(ind,m)
  }
} 
stack<-stack[,,ind]
stack.UB<-stack.UB[,,ind]

stack.temp<-stack[,,1]
display(stack.temp,"raster")







# #display(stack.UB*10,"raster",all=T)
#
#
stack.temp.ub<-stack.UB[,,1]
display(stack.temp.ub*2,"raster")

# roi1.red = thresh(stack.temp.ub, 100, 100,0.004)
# display(colorLabels(roi1.red),"raster")
roi2.red = thresh(stack.temp.ub, 4, 4,0.02)
display(roi2.red,"raster")
R<-ifelse(roi2.red > 0,stack.temp.ub,0)
display(R*10,"raster")

roi2 = thresh(stack.temp, 4, 4,0.15)
G<-ifelse(roi2 > 0,stack.temp,0)
display(G*10,"raster")
Gcol<-ifelse(R > 0,G,0)
Gsum<-ifelse(G > 0,G,0)
m2<-(sum(as.vector(Gcol)))/(sum(as.vector(G)))

#     
#     Rcol<-ifelse(stack.temp)
#     overlap.roi<-roi2 + roi2.red
#     overlap.G<-ifelse(overlap.roi > 0 ,stack.temp,0)
#     display(overlap.G,"raster")
# 
#     
#     Green.inten<-computeFeatures.basic(roi1,stack.temp)[1]
# 
# 
# 
# 
# 
# m1<-(computeFeatures.basic(roi2.red,stack.temp)[1])/(computeFeatures.basic(roi2,stack.temp)[1])
# m2<-(computeFeatures.basic(roi2,stack.temp)[1])/(computeFeatures.basic(roi2,stack.temp.ub)[1])
# m<-m1*m2
# 
# dev.new("red")
# display(stack,"raster",all=T)
# dev.new("red")
# display(stack.UB*5,"raster",all=T)
# dev.new("red")
# display(GFP*5,"raster",all=T)
# dev.new("red")
# display(UB*20,"raster",all=T)
