
bikes=read.csv("C:/Users/samjb/Desktop/Nonparametric/hw1/train.csv")

proc=bikes;
proc$season=factor(proc$season)
levels(proc$season)<-c("spring","summer","fall","winter")
proc$count=log(proc$count+1)
proc$hour=factor(proc$hour)
proc$workingday=factor(proc$workingday)
proc$holiday=factor(proc$holiday)

hdata<-proc[(proc$day<=15),]
validate<-proc[(proc$day>15),]

#basemodel=lm(count~temp*atemp + windspeed + factor(season) + factor(holiday) + factor(workingday) + factor(hour),data=hdatra
basemodel=lm(count~(month+factor(workingday)+temp+humidity+year+factor(hour))+temp:year+month:humidity+workingday:temp+workingday:humidity+workingday:hour,data=hdata)
mypredict=predict(basemodel,validate)
RMSLE=sqrt(mean((validate$count-mypredict)^2))
RMSLE

hdata2=aggregate(.~daylabel,data=hdata,mean)
plot(hdata2$daylabel,hdata2$count,xlab="DayLabel",ylab="Average Hourly Count")

xs=hdata2$daylabel
Ys=hdata2$count

#hs=(1:15)*10;
#risks=as.numeric(lapply(hs,test))
#plot(hs,risks,type="b")
#sum(modelnp$residuals^2)
fitted_values=test(40);

plot(xs,Ys,xlab="daylabels with fitted mean line, bandwidth=40",ylab="log counts")
lines(xs,fitted_values,col="red",lwd=2)
modelnp=loess(count~daylabel,data=hdata2,degree=1,span=.12,control=loess.control(surface="direct"))
lines(hdata2$daylabel,modelnp$fitted,type="l",col="blue")
lines(hdata2$daylabel,fitted_values,type="l",col="red")

fitted_values=modelnp$fitted

hdata3=hdata
xs=hdata2$daylabel
for(i in (1:length(xs)))
{
  lab=xs[i];
  hdata3$count[hdata3$daylabel==lab]=hdata3$count[hdata3$daylabel==lab]-fitted_values[i];
}

model_resids=lm(count~(month+workingday+temp+humidity+hour+year)+temp:year+month:humidity+workingday:temp+workingday:humidity+workingday:hour,data=hdata3)

predict_resids=predict(basemodel3,validate)
predict_means=predict(modelnp,validate)
val_predicts=predict_resids+predict_means
RMSLE=sqrt(mean((validate$count-val_predicts)^2))
RMSLE

mygam=gam(count~(daylabel+temp+factor(hour)+humidity+year)+daylabel:humidity+daylabel:temp+daylabel:humidity+hour:humidity+humidity:year,data=hdata);
summary(mygam)
gam_predict=predict(mygam,validate)
RMSLE=sqrt(mean((validate$count-gam_predict)^2))
RMSLE


w=c(1,1,1)
w=w/(sum(w))
avg_predict=(w[1]*mypredict+w[2]*val_predicts+w[3]*gam_predict)
RMSLE=sqrt(mean((validate$count-avg_predict)^2))
RMSLE


########USE MODELS##########


test_unproc=read.csv("C:/Users/samjb/Desktop/Nonparametric/hw1/test.csv");
test_set=test_unproc;
test_set$season=factor(test_set$season)
test_set$hour=factor(test_set$hour)
levels(test_set$season)<-c("spring","summer","fall","winter")
test_set$workingday=factor(test_set$workingday)
test_set$holiday=factor(test_set$workingday)

train_unproc=read.csv("C:/Users/samjb/Desktop/Nonparametric/hw1/train.csv");
train_set=train_unproc;
train_set$season=factor(train_set$season)
train_set$hour=factor(train_set$hour)
levels(train_set$season)<-c("spring","summer","fall","winter")
train_set$count=log(train_set$count+1)
train_set$workingday=factor(train_set$workingday)
train_set$holiday=factor(train_set$workingday)

basemodel_train=lm(count~(month+workingday+temp+humidity+hour+year)+temp:year+month:humidity+workingday:temp+workingday:humidity+workingday:hour,data=train_set)
test_predict1=predict(basemodel,test_set)

train_set_agg=aggregate(.~daylabel,data=train_set,mean);
model_mean=loess(count~daylabel,data=train_set_agg,degree=1,span=.12,control=loess.control(surface="direct"))

train_set_resids=train_set
xs=train_set_agg$daylabel
for(i in (1:length(xs)))
{
  lab=xs[i];
  train_set_resids$count[train_set_resids$daylabel==lab]=train_set_resids$count[train_set_resids$daylabel==lab]-model_mean$fitted[i];
}

model_resids=lm(count~(month+workingday+temp+humidity+hour+year)+temp:year+month:humidity+workingday:temp+workingday:humidity+workingday:hour,data=train_set_resids)

test_predict_resids=predict(model_resids,test_set)
test_predict_means=predict(model_mean,test_set)
test_predict2=test_predict_resids+test_predict_means


model_gam=gam(count~(daylabel+temp+hour+humidity+year)+daylabel:humidity+daylabel:temp+daylabel:humidity+hour:humidity+humidity:year,data=train_set);
summary(model_gam)
test_predict3=predict(model_gam,test_set)

avg_predict=(test_predict1+test_predict2+test_predict3)/3

recover_1=exp(test_predict1)-1
recover_2=exp(test_predict2)-1
recover_3=exp(test_predict3)-1

recover_counts=exp(avg_predict)-1
cat(recover_counts)

setwd("C:/Users/samjb/Desktop/Nonparametric/hw1/")
sink("assn1-samjbaugh.txt",append=TRUE)
cat(t(recover_counts),sep="\n")
closeAllConnections()