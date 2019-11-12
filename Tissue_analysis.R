library(tictoc)
library(EBImage)
library(data.table)
library(outliers)
setwd("~/Desktop/temp/Shireen/Expirement_3/Tisue_mesurments/090919/Striatum/Processed/WT/")
tifs<-dir()[grep(".tif",dir())]
make.mask<-function(tifs){
  ll<-list()
  #THersholding
  #frame<-readImage(tifs[1],type = "tif",all = T)
  frame<-readImage(tifs,type = "tif",all = T)
  #frame<-channel(frame,"gray")
  dimm<-dim(frame)
  GFP<-frame[1:dimm[1],1:dimm[2],1]
  colorMode(GFP)<-"Grayscale"
  dapi_normal<- GFP
  #display(GFP,"raster")
  nmask2 = thresh(dapi_normal, 20, 20,0.04)
  mk3 = makeBrush(7, shape= "diamond")
  nmask2 = opening(nmask2, mk3)
  nmask2 = fillHull(nmask2)
  seg = paintObjects(nmask2,toRGB(GFP),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  # dev.off()
  #display(seg,"raster")
  nseg = bwlabel(nmask2)  #binery to object
  rm(nmask2)
  chackpoint<-computeFeatures.shape(nseg)
  nr = which(chackpoint[,2] > 150)
  nseg = rmObjects(nseg, nr)
  seg = paintObjects(nseg,toRGB(GFP),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  # dev.off()
  #display(seg,"raster")
  ll[[paste0("frame")]]<-GFP
  ll[[paste0("mask")]]<-nseg
  return(ll)
}

list.df<-list()
list.seg<-list()
for (i in 1:length(tifs)){
        WT1<-try(make.mask(tifs[i]))
        if (inherits(WT1, "try-error")){
          df.temp<-data.frame("name"=paste0(tifs[i]),"sum_area"=0,"number_of_foci"=0)
          list.df[[i]]<-df.temp
          list.seg[[paste0(tifs[i])]]<-seg
          next()
        }
        frame<-WT1$frame
        dim(frame)
        newframe<-frame
        display(newframe)
        nseg<-WT1$mask
        newmask<-nseg
        seg = paintObjects(newmask,toRGB(newframe),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
        #dev.off()
        display(seg,"raster")
        nf.ctrl.sum.area = sum(computeFeatures.shape(newmask)[,1])
                if (nf.ctrl.sum.area==0){
                  nf.ctrl.sum.area<-0.000001
                }
                nr.ctrl.mean.intensity = mean(computeFeatures.basic(newmask,newframe)[,1])
                    if (is.na(nr.ctrl.mean.intensity)){
                      nr.ctrl.mean.intensity<-0.000001
                    }
        if (nf.ctrl.sum.area==0.000001){
             df.temp<-data.frame("name"=paste0(tifs[i]),"sum_area"=nf.ctrl.sum.area,"number_of_foci"=0.000001,stringsAsFactors = F)
                }else{
          df.temp<-data.frame("name"=paste0(tifs[i]),"sum_area"=nf.ctrl.sum.area,"number_of_foci"=nrow(computeFeatures.shape(newmask)),stringsAsFactors = F)
          }
        list.df[[i]]<-df.temp
        list.seg[[paste0(tifs[i])]]<-seg
}

length(list.df)
df<-data.frame("name"=NULL,"sum_area"=  NULL,"number_of_foci"=NULL,stringsAsFactors = F)
df.temp<-NULL
for (i in 1:length(list.df)){
    df.temp<- list.df[[i]]
    df<-rbind(df,df.temp)
}

write.csv(x = df,file = "WT_area.csv")




#1352
#118902

WT1<-make.mask(tifs[1])
frame<-WT1$frame
dim(frame)
newframe<-frame
display(newframe*3,"raster")
nseg<-WT1$mask
newmask<-nseg
seg = paintObjects(newmask,toRGB(newframe*3),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
dev.off()
display(seg,"raster")
nf.ctrl.sum.area = sum(computeFeatures.shape(newmask)[,1])
if (nf.ctrl.sum.area==0){
  nf.ctrl.sum.area<-0.000001
}
nr.ctrl.mean.intensity = mean(computeFeatures.basic(newmask,newframe)[,1])
if (is.na(nr.ctrl.mean.intensity)){
  nr.ctrl.mean.intensity<-0.000001
}
nrow(computeFeatures.shape(newmask))
#1315
#109765


#Test intensity per obgect or number of object - runing window of 200 x 200


array<-seq(200,1500,1)
centerx<-sample(array,1)
centery<-sample(array,2)
x<-centerx-100
xend<-centerx+100
y<-centery-100
yend<-centery+100





WT1<-make.mask(tifs[2])
frame<-WT1$frame
dim(frame)
newframe<-frame[x:xend,y:yend]
display(newframe*3)
nseg<-WT1$mask
newmask<-nseg[x:xend,y:yend]
seg = paintObjects(newmask,toRGB(newframe*10),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
dev.off()
display(seg,"raster")
nf.ctrl.sum.area = sum(computeFeatures.shape(newmask)[,1])
if (nf.ctrl.sum.area==0){
  nf.ctrl.sum.area<-0.000001
}
nr.ctrl.mean.intensity = mean(computeFeatures.basic(newmask,newframe)[,1])
if (is.na(nr.ctrl.mean.intensity)){
  nr.ctrl.mean.intensity<-0.000001
}



tax<-make.mask(tifs[1])
frame.tax<-tax$frame
newframe.tax<-frame.tax[x:xend,y:yend]
display(newframe.tax*10)
nseg.tax<-tax$mask
newmask.tax<-nseg.tax[x:xend,y:yend]
seg = paintObjects(newmask.tax,toRGB(newframe.tax*10),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
#dev.off()
display(seg,"raster")
nf.tax.sum.area = sum(computeFeatures.shape(newmask.tax)[,1])
if (nf.tax.sum.area==0){
  nf.tax.sum.area<-0.000001
}
nr.tax.mean.intensity = mean(computeFeatures.basic(newmask.tax,newframe.tax)[,1])
if (is.na(nr.tax.mean.intensity)){
  nr.tax.mean.intensity<-0.000001
}
vec.tax<-



  nper<-137
diff_obs<-nr.tax.mean.intensity-nr.ctrl.mean.intensity
diff<-rep(NA,nper)
diff[i]<-diff_obs


for (i in 1:nper){
  shafulled_labales<-sample(groups, replace = F)
  diff[i]<-mean(data[shafulled_labales==1])-mean(data[shafulled_labales==2])
}
pv<-(length(diff[abs(diff) >= abs(diff_obs)]))/nper
list.pv[[paste0(v)]]<-cbind(v,pv)





frame<-list$frame
newframe<-frame[x:xend,y:yend]
display(newframe*10)
nseg<-list$mask
newmask<-nseg[x:xend,y:yend]
seg = paintObjects(newmask,toRGB(newframe*10),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
dev.off()
display(seg,"raster")
nf.ctrl.sum.area = sum(computeFeatures.shape(newmask)[,1])
nr.ctrl.mean.intensity = mean(computeFeatures.basic(newmask,newframe)[,1])
if (is.na(nr.ctrl.mean.intensity)){next}

nper<-137
groups <- c(rep(1, l), rep(2, l))
data<-c(temp.sortcount,rand.temp)
diff_obs<-(mean(temp.sortcount))-(mean(rand.temp))
diff<-rep(NA,nper)
for (i in 1:nper){
  shafulled_labales<-sample(groups, replace = F)
  diff[i]<-mean(data[shafulled_labales==1])-mean(data[shafulled_labales==2])
}
pv<-(length(diff[abs(diff) >= abs(diff_obs)]))/nper
list.pv[[paste0(v)]]<-cbind(v,pv)




nmask<-nseg
nf = computeFeatures.shape(nmask)
nr = which(nf[,2] < 120)
nseg = rmObjects(nmask, nr)
nn = max(nseg)
colorMode(nseg)<-Grayscale
chackpoint<-computeFeatures.shape(nseg)
int.dapi<-computeFeatures.basic(nseg,dapi)
y<-which(scores(int.dapi[,1], type="z",prob = 0.95))
tp<-as.numeric(attr(y,"names"))
nseg<-rmObjects(nseg,tp)
chackpoint<-computeFeatures.shape(nseg)
df<-as.data.frame(chackpoint)
rm(chackpoint,nmask)
xy<-computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
if (is.null(nrow(xy)) || is.null(xy)){
  next
}
df<-cbind(df,xy)
df.combine<-as.data.frame(matrix(0,nrow(xy),5))
colnames(df.combine)<-c("x","y","Area_real","Areal_roundess","ratio")
df.combine$x<-xy[,1]
df.combine$y<-xy[,2]
df.combine$Area_real<-df[,1] #area of sample
df.combine$Areal_roundess<-pi*(df[,3])^2
df.combine$ratio<-df.combine[,4]/df[,1]
nr = which(df.combine[,5] > 1.3 )
gsegg = rmObjects(nseg, nr)
rm(nseg)
#rm(nseg,df)
nr = which(df.combine[,5] < 0.8 )
#rm(df.combine)
gsegg = rmObjects(gsegg, nr)
rm(nr,df.combine)
