require('worldmap')
require('fields')
#setwd('C:/Users/samjb/Desktop/Nonparametric/Final Project')

station_data=read.csv('master-location-identifier-database-20130801.csv',header=TRUE)
aus_data_na=station_data[station_data$country3=="AUS",]
aus_data=aus_data_na[is.na(aus_data_na$lon_prp)+is.na(aus_data_na$lon_prp)==0,]

aus_means=read.table('meanann.txt')

xgrid=seq(110,155,length=1361)
ygrid=seq(-40,-10,length=1681)
latlongs=matrix(0,1361,1681)

n=length(aus_data$country3)
stations=data.frame(lat=rep(0,n),long=rep(0,n),temperature=rep(0,n))
for(ii in 1:n)
{
  stations$lat[ii]=aus_data$lat_prp[ii]
  stations$long[ii]=aus_data$lon_prp[ii]
}

xlims=c(111.2,154.8)
ylims=c(-40.4,-8.6)

stations=stations[stations$lat>ylims[1] & stations$lat<ylims[2],]
stations=stations[stations$long>xlims[1] & stations$long<xlims[2],]

yindices=seq(0,1361-170,by=1)
xindices=seq(1,1681,by=1)

xgrid=seq(xlims[1],xlims[2],length=length(xindices))
ygrid=seq(ylims[1],ylims[2],length=length(yindices))
aus_matrix=as.matrix(aus_means)
N=dim(aus_means)[1]
M=dim(aus_means)[2]
temp_matrix=aus_matrix+matrix( rnorm(N*M,mean=0,sd=1), N, M) 
sub_mat_bad=temp_matrix[yindices,xindices]
sub_mat <- t(apply(sub_mat_bad, 2, rev))
sub_mat[sub_mat<(-9000)] <- NA

image.plot(xgrid,ygrid,sub_mat)
points(lat~long,data=stations)

n=length(stations$lat)
for(ii in 1:n)
{
  mylat=stations$lat[ii]
  mylon=stations$long[ii]
  
  indexLatY=min(which(ygrid > mylat))
  indexLonX=min(which(xgrid > mylon))
  
  stations$temperature[ii]=sub_mat[indexLonX,indexLatY]
}
stations=stations[!is.na(stations$temperature),]


rbPal <- colorRampPalette(c('blue','red'))
stations$colors <- rbPal(10)[as.numeric(cut(stations$temperature,breaks = 10))]

set.seed(123)
m=1000
indices=sample(1:length(stations$lat),m)
stations_sub=stations[indices,]
test_set=stations[-indices,]

layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="1000 Training Stations")
points(lat~long,data=stations_sub,col=colors,cex=.9,pch=20)

legend_image <- as.raster(matrix(rev(rbPal(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(stations$temperature),max(stations$temperature),l=5))
labs=paste(labs,'C')
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))


locmat=matrix(0,m,2)
locmat[,1]=(stations_sub$lat-min(stations_sub$lat))/max(stations_sub$lat-min(stations_sub$lat))
locmat[,2]=(stations_sub$long-min(stations_sub$long))/max(stations_sub$long-min(stations_sub$long))

coordmat=matrix(0,m,2)
coordmat[,1]=stations_sub$lat
coordmat[,2]=stations_sub$long
v1=variog(coords=coordmat,data=stations_sub$temperature)
#1 is latitude, 2 longitude

for(nu in seq(.5,3,by=.5))
{
  ptm <- proc.time()
  model1=MLESpatialProcess(locmat,stations_sub$temperature,cov.args=list(Covariance = "Matern",smoothness=nu))
  print(proc.time() - ptm)
  print(model1$pars)
}
gridsizex=50
gridsizey=50
xindices=round(seq(1,1681,length=gridsizex))
yindices=round(seq(1,1361-170,length=gridsizey))
xgrid=seq(xlims[1],xlims[2],length=gridsizex)
ygrid=seq(ylims[1],ylims[2],length=gridsizey)
sub_mat2=aus_matrix[yindices,xindices]
sub_mat2 <- t(apply(sub_mat2, 2, rev))
sub_mat2[sub_mat2==-9999] <- NA

model1=Krig(locmat,stations_sub$temperature, Covariance="Matern", theta=0.05935557, smoothness=1.5, sigma=1.05287098, rho=6.30118069)  

ygridM=as.numeric(unlist(lapply(ygrid,rep,length(xgrid))))
xgridM=rep(xgrid,length(ygrid))
predgrid=matrix(0,length(xgrid)*length(ygrid),2)
xnorm=(xgridM-min(xgridM))/(max(xgridM)-min(xgridM))
ynorm=(ygridM-min(ygridM))/(max(ygridM)-min(ygridM))
predgrid[,1]=ynorm
predgrid[,2]=xnorm

fitted_values=predict(model1,predgrid)
fit_matrix=matrix(fitted_values,nrow=gridsizex,ncol=gridsizey)
fit_matrix=fit_matrix*!is.na(sub_mat2)
fit_matrix[fit_matrix==0]=NA
image.plot(xgrid,ygrid,fit_matrix,main="Gridded Fitted Values",xlab="Longitude",ylab="Latitude")

predictionErrors=predictSE(model1,predgrid)
conf_low_vec=fitted_values-predictionErrors*qnorm(.9725)
conf_low_mat=matrix(conf_low_vec,nrow=gridsizex,ncol=gridsizey)
conf_low_mat=conf_low_mat*!is.na(sub_mat2)
conf_low_mat[conf_low_mat==0]=NA
image.plot(xgrid,ygrid,conf_low_mat,main="Lower 95%-Confidence Map",xlab="Longitude",ylab="Latitude")

SE_mat=matrix(predictionErrors,nrow=gridsizex,ncol=gridsizey)
SE_mat=SE_mat*!is.na(sub_mat2)
SE_mat[SE_mat==0]=NA
image.plot(xgrid,ygrid,SE_mat,main="Standard Errors",xlab="Longitude",ylab="Latitude")


conf_up_vec=fitted_values+predictionErrors*qnorm(.9725)
conf_up_mat=matrix(conf_up_vec,nrow=gridsizex,ncol=gridsizey)
conf_up_mat=conf_up_mat*!is.na(sub_mat2)
conf_up_mat[conf_up_mat==0]=NA
image.plot(xgrid,ygrid,conf_up_mat,main="Upper 95%-Confidence Map",xlab="Longitude",ylab="Latitude")

testgrid=matrix(0,length(test_set$lat),2)
theta1=1
theta2=1.01
testgrid[,1]=(test_set$lat-min(test_set$lat))/(theta1*(max(test_set$lat)-min(test_set$lat)))
testgrid[,2]=(test_set$long-min(test_set$long))/(theta2*(max(test_set$long)-min(test_set$long)))
test_values=predict(model1,testgrid)
#test_SEs=predictSE(model1,testgrid)
#test_conf_up=test_values+test_SEs*qnorm(.9725)
#test_conf_low=test_values-test_SE*qnorm(.9725)
#ch=seq(min(test_conf_low),max(test_conf_up),length=11)
which_ch=function(x)
{
  return(min(which(x<ch))-1) 
}

####PLOT TEST RESIDUALS
greGrad <- colorRampPalette(c('white','red'))
resSquare=abs(test_values-test_set$temperature)
print(mean(resSquare^2))

testcolors = greGrad(10)[as.numeric(cut(resSquare,breaks = 10))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Absolute Value Empirical Errors (Test Set)")
points(lat~long,data=test_set,col=testcolors,cex=1.4,pch=20)
#points(lat~long,data=test_set,col='black',cex=1,pch=21)
legend_image <- as.raster(matrix(rev(greGrad(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(resSquare),max(resSquare),l=5))
labs=paste(labs,'C')
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))
mean(resSquare^2)


####PLOT TEST LOWER CONFS
#testcolors1 = rbPal(10)[apply(test_conf_low,1,which_ch)]
#layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
#plot(newmap,xlim=xlims,ylim=ylims,main="Test Lower 95% Confidence")
#points(lat~long,data=test_set,col=testcolors1,cex=.9,pch=20)
#legend_image <- as.raster(matrix(rev(rbPal(10)), ncol=1))
#plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
#labs=round(seq(min(ch),max(ch),l=5))
#labs=paste(labs,'C')
#text(x=1.5, y = seq(0,1,l=5), labels = labs)
#rasterImage(legend_image, 0, 0, 1,1)
#layout(matrix(1:1,ncol=1), width = c(1),height = c(1))

####PLOT TEST VALUES
testcolors2 = rbPal(10)[apply(test_SEs,1,which_ch)]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Predictions on Test Set")
points(lat~long,data=test_set,col=testcolors2,cex=.9,pch=20)
legend_image <- as.raster(matrix(rev(rbPal(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(ch),max(ch),l=5))
labs=paste(labs,'C')
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))

####PLOT STANDARD ERRORS
greGrad <- colorRampPalette(c('white','blue'))
testcolors = greGrad(10)[as.numeric(cut(test_SEs,breaks = 10))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Test Set Standard Errors (Theoretical)")
points(lat~long,data=test_set,col=testcolors,cex=1.4,pch=20)
#points(lat~long,data=test_set,col='black',cex=1,pch=21)
legend_image <- as.raster(matrix(rev(greGrad(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(test_SEs),max(test_SEs),l=5),2)
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))


###NONPARAMETRIC PORTION
spline=spline(stations_sub$temperature~stations_sub$lat+stations_sub$long)
image.plot(xgrid,ygrid,spline,xlab="longitude",ylab="latitude",main="Local Linear Smoothing")

meanmodel=loess(stations_sub$temperature~stations_sub$lat+stations_sub$long,span=10,degree=2)
mean_fitted=predict(meanmodel,predgrid_s)
mean_matrix=matrix(mean_fitted,nrow=gridsizex,ncol=gridsizey)
mean_matrix=mean_matrix*!is.na(sub_mat2)
mean_matrix[mean_matrix==0]=NA
image.plot(xgrid,ygrid,mean_matrix,xlab="longitude",ylab="latitude",main="Mean Function for Matern Model")


modelloess=loess(stations_sub$temperature~stations_sub$lat+stations_sub$long,span=.05,degree=1)

predgrid_s=matrix(0,length(xgrid)*length(ygrid),2)
predgrid_s[,1]=ygridM
predgrid_s[,2]=xgridM

gridPredict=predict(modelloess,predgrid_s,se=T)
fitted_loess=gridPredict$fit #all the same value
loess_matrix=matrix(fitted_loess,nrow=gridsizex,ncol=gridsizey)
loess_matrix=loess_matrix*!is.na(sub_mat2)
loess_matrix[loess_matrix==0]=NA
image.plot(xgrid,ygrid,loess_matrix,xlab="longitude",ylab="latitude",main="Local Linear Smoothing")

loessSE=gridPredict$se.fit #all the same value
loessSE_matrix=matrix(loessSE,nrow=gridsizex,ncol=gridsizey)
loessSE_matrix=loessSE_matrix*!is.na(sub_mat2)
loessSE_matrix[loessSE_matrix==0]=NA
loessSE_matrix[loessSE_matrix>.5]=.5
image.plot(xgrid,ygrid,loessSE_matrix,xlab="longitude",ylab="latitude",main="Local Linear Standard Errors")



###Test Data
testgrid=matrix(0,length(test_set$lat),2)
testgrid[,1]=test_set$lat
testgrid[,2]=test_set$long

test_set_loess=test_set
mypredict=predict(modelloess,testgrid,se=T)
test_set_loess$pred=mypredict$fit
test_set_loess$se=mypredict$se.fit
test_set_loess=test_set_loess[!is.na(test_set_loess$pred),]

greGrad <- colorRampPalette(c('blue','red'))
testcolors = greGrad(10)[as.numeric(cut(test_set_loess$pred,breaks = 10))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Test Values: Local Linear Regression")
points(lat~long,data=test_set_loess,col=testcolors,cex=1.4,pch=20)
legend_image <- as.raster(matrix(rev(greGrad(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(test_set_loess$pred),max(test_set_loess$pred),l=5))
labs=paste(labs,'C')
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))


greGrad <- colorRampPalette(c('white','red'))
resSquare=abs(test_set_loess$pred-test_set_loess$temperature)
testcolors = greGrad(10)[as.numeric(cut(resSquare,breaks = 10))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Absolute Value Test Set Residuals")
points(lat~long,data=test_set,col=testcolors,cex=1.4,pch=20)
legend_image <- as.raster(matrix(rev(greGrad(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(resSquare),max(resSquare),l=5))
labs=paste(labs,'C')
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))
mean(resSquare^2)

greGrad <- colorRampPalette(c('white','blue'))
testcolors = greGrad(10)[as.numeric(cut(test_set_loess$se,breaks = 10))]
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(newmap,xlim=xlims,ylim=ylims,main="Standard Errors (Theoretical)")
points(lat~long,data=test_set,col=testcolors,cex=1.4,pch=20)
legend_image <- as.raster(matrix(rev(greGrad(10)), ncol=1))
plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Gradient Legend')
labs=round(seq(min(test_set_loess$se),max(test_set_loess$se),l=5),2)
text(x=1.5, y = seq(0,1,l=5), labels = labs)
rasterImage(legend_image, 0, 0, 1,1)
layout(matrix(1:1,ncol=1), width = c(1),height = c(1))

modelspline=mars(locmat,stations_sub$temperature)
fitted_values=predict(modelspline,predgrid)
fit_matrix=matrix(fitted_values,nrow=gridsizex,ncol=gridsizey)
fit_matrix=fit_matrix*!is.na(sub_mat2)
fit_matrix[fit_matrix==0]=NA
image.plot(xgrid,ygrid,fit_matrix,main="MARS Spline Fit",xlab="Longitude",ylab="Latitude")


