args <- commandArgs(trailingOnly = TRUE)
library(EBImage)
library("gtools")
library(outliers)
frame<-readImage(args[1],type = "tif",all = T)
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
#display(UB*20,"raster")
#display(tax*20,"raster")



dapi_normal<- dapi*5
#display(dapi_normal,"raster")
nmask2 = thresh(dapi_normal, 85, 85,0.003)
mk3 = makeBrush(13, shape= "diamond")
nmask2 = opening(nmask2, mk3) 
nmask2 = fillHull(nmask2) 
nseg = bwlabel(nmask2)  #binery to object
seg_pink = paintObjects(nseg,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
#display(nseg,"raster")
rm(nmask2)
chackpoint<-computeFeatures.shape(nseg)
nmask<-nseg
if (length(chackpoint) < 3 ){
  stop()
}
if (is.null(nrow(chackpoint))){
  stop()
}
# if (nrow(chackpoint) > 100){
#   stop()
# }
nf = computeFeatures.shape(nmask)
nr = which(nf[,2] < 50) 
nseg = rmObjects(nmask, nr) 
nn = max(nseg) 
colorMode(nseg)<-Grayscale
chackpoint<-computeFeatures.shape(nseg)
if (length(chackpoint) < 3 ){
  stop()
}
if (is.null(nrow(chackpoint))){
  stop()
}
if (nrow(chackpoint) > 2000 ){
  stop()
}
int.dapi<-computeFeatures.basic(nseg,dapi)
y<-which(scores(int.dapi[,1], type="z",prob = 0.975))
tp<-as.numeric(attr(y,"names"))
if (length(tp) < 1 ){
  stop()
}
nseg<-rmObjects(nseg,tp)
chackpoint<-computeFeatures.shape(nseg)
df<-as.data.frame(chackpoint)
rm(chackpoint,nmask)
xy<-computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
if (length(xy) < 3 ){
  stop()
}
if (is.null(nrow(xy))){
  stop()
}
if (nrow(xy) > 500 ){
  stop()
}
gsegg=nseg
#display(colorLabels(gsegg),"raster")

cell_normal<- GFP*10
#display(cell_normal,"raster")
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
##display(combine,method="raster")
combine[gsegg > combine]<-gsegg[gsegg > combine]
##display(colorLabels(combine),method="raster")
##display(colorLabels(gsegg),method="raster")

open_2 = paintObjects(combine,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)

csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=cmask)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
##display(colorLabels(csegpink),method="raster")

csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=combine)
#rm(open_2,cmask)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
##display(colorLabels(csegpink),method="raster")
xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
if (length(xy) < 3 ){
  stop()
}
if (is.null(nrow(xy))){
  stop()
}
if (nrow(xy) > 2000 ){
  stop()
}
cf = computeFeatures.shape(csegpink)
cf = computeFeatures.shape(csegpink)
#cr = which(cf[,2] < 200) 
#csegpink = rmObjects(csegpink, cr,reenumerate = F) 
cf = computeFeatures.shape(csegpink)
cfErea<-data.frame(cf[,1])
cfErea$num<-row.names(cfErea)
ci = which(cf[,1] > 100000) 
csegpink = rmObjects(csegpink, ci,reenumerate = F) 
rm(cf,ci,cfErea)
#' select all the cells with normal ilumination propteis
int.GFP<-computeFeatures.basic(csegpink,GFP)
#y<-which(scores(int.GFP[,1], type="iqr"))
y<-which(scores(int.GFP[,1], type="z",prob = 0.95))
seg_pink = paintObjects(csegpink,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
#display(seg_pink,method="raster")
ty<-as.numeric(attr(y,"names"))
csegpink<-rmObjects(csegpink,ty,reenumerate = F)
seg_pink = paintObjects(csegpink,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
#display(seg_pink,method="raster")




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
#display(colorLabels(csegpink),method="raster")

#display(GFP*10,method="raster")

# test<-UB
# test[gsegg>0]<-0
# #display(test*10,method="raster")
# summary(as.vector(test))

UB[gsegg>0]<-0
GFP[gsegg>0]<-0

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
df.mcc<-NULL
for (m in 1:dim(stack)[3]){
  stack.temp<-stack[,,m]
  #display(stack.temp,"raster")
  stack.temp.ub<-stack.UB[,,m]
  #display(stack.temp.ub*2,"raster")
  roi2.red = thresh(stack.temp.ub, 4, 4,mean(as.vector(stack.temp.ub)))
  #display(roi2.red,"raster")
  R<-ifelse(roi2.red > 0,stack.temp.ub,0)
  #display(R*10,"raster")
  roi2 = thresh(stack.temp, 4, 4,mean(as.vector(stack.temp)))
  #display(roi2,"raster")
  G<-ifelse(roi2 > 0,stack.temp,0)
  #display(G*10,"raster")
  Gcol<-ifelse(R > 0,G,0)
  Gsum<-ifelse(G > 0,G,0)
  Rcol<-ifelse(G > 0,R,0)
  Rsum<-ifelse(R > 0,R,0)
  m2<-(sum(as.vector(Gcol)))/(sum(as.vector(G)))
  m1<-(sum(as.vector(Rcol)))/(sum(as.vector(R)))
  image_name<-paste0(m,"_stack",args[1])
  meanmcc<-(m1+m2)/2
  gfp.i.stack<-mean(as.vector(stack.temp))
  df.mcc<-rbind(df.mcc,c(image_name,m1,m2,m1*m2,meanmcc,gfp.i.stack))
  rm(m2,m1,image_name,meanmcc,gfp.i.stack)
  # setwd(out)
  # 
  # tiff(paste0(m,"GFP_stack_",tif[i]))
  # #display(stack.temp,"raster")
  # dev.off()
  # tiff(paste0(m,"UB_stack_",tif[i]))
  # #display(stack.temp.ub,"raster")
  # dev.off()
  # 
}
df.mcc<-as.data.frame(df.mcc,stringsAsFactors = F)
colnames(df.mcc)<-c("image_name","m1","m2","MCC","meanMCC","GFP_intensity")
write.table(df.mcc,file = paste0(args[2]),row.names = F)





