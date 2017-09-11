mydata<-read.table("./data.txt", header=TRUE, row.names="Year")

altdat=c()

wint=mydata$"DJF"[-1]
spring=mydata$"MAM"[-358]
summ=mydata$"JJA"[-358]
fall=mydata$"SON"[-358]

altdat=c()
for(i in (1:358))
{
  if(i!=1)
  {
  altdat=c(altdat,mydata$"DJF"[i])
  }
  if(i!=358)
  {
  altdat=c(altdat,mydata$"MAM"[i])
  altdat=c(altdat,mydata$"JJA"[i])
  altdat=c(altdat,mydata$"SON"[i])
  }
  
}

#Idea 1

plot(diff(altdat),type="l")

plot(diff(mydata$"DJF"[-1]),type="l")
plot(diff(mydata$"MAM"[-358]),type="l")
plot(diff(mydata$"JJA"[-358]),type="l")
plot(diff(diff(mydata$"SON"[-358])),type="l")



avcumsum <- function(season){
  sum=season[1]
  avcumsum=c()
  for(i in (2:357))
  {
    sum=sum+season[i]
    avcumsum=c(avcumsum,sum/i)
  }
  return(avcumsum)
}

season=wint
plot(avcumsum(season),type="l")
fit=lm(avcumsum(season)~years)
abline(fit$coefficients[1],fit$coefficients[2])


#Chapter 1: Visualizing Data (Symmetric Moving Average)

years=c(1:358)+1658

seasons=c()

for(i in (4:1432))
{
  seasons=c(seasons,floor(years[floor(i/4)]))
}
seasons=seasons[-1]

par(mfrow=c(1,1))

plot(seasons,altdat,type="l",xlab="Year",ylab="Degrees Celcius")

v=filter(rowMeans(mydata)[-1][-357],sides=2,filter=rep(1/25,25))
plot(years[-1][-357],v,type="l",xlab="Year",ylab="Degrees Celcius",main="Symmetric Moving Average: All Seasons")
abline(v=1726)
abline(v=1975)


par(mfrow=c(2,2))

v=filter(wint,sides=2,filter=rep(1/25,25))
plot(years[-1],v,type="l",xlab="Year",ylab="Degrees Celcius",main="Symmetric Moving Average: Winter")

v=filter(spring,sides=2,filter=rep(1/25,25))
plot(years[-358],v,type="l",xlab="year",ylab="Degrees Celcius",main="Symmetric Moving Average: Spring")

v=filter(summ,sides=2,filter=rep(1/25,25))
plot(years[-358],v,type="l",xlab="year",ylab="Degrees Celcius",main="Symmetric Moving Average: Summer")

v=filter(fall,sides=2,filter=rep(1/25,25))
plot(years[-358],v,type="l",xlab="year",ylab="Degrees Celcius",main="Symmetric Moving Average: Fall")


plot.ts(years[-1][-357],rowMeans(mydata),main="Temperature Data",xlab="Year",ylab="Degrees Celcius",type="l")

par(mfrow=c(2,2))

plot.ts(years[-1],wint,main="Temperature Data (Winter)",xlab="Year",ylab="Degrees Celcius",type="l")
plot.ts(years[-357],spring,main="Temperature Data (Spring)",xlab="Year",ylab="Degrees Celcius",type="l")
plot.ts(years[-357],summ,main="Temperature Data (Summer)",xlab="Year",ylab="Degrees Celcius",type="l")
plot.ts(years[-357],fall,main="Temperature Data (Fall)",xlab="Year",ylab="Degrees Celcius",type="l")


#Part 2: Basic Models

summary(fit <- lm(temp_data~years[-1][-357]))

plot(years[-1][-357],resid(fit),type="l",xlab="Year",ylab="Degrees Celcius",main="Detrended Temperatures")


v=filter(rowMeans(mydata),sides=2,filter=rep(1/25,25))
plot(years[-1][-357],v,type="l",xlab="Year",ylab="Degrees Celcius",main="Symmetric Moving Average: All Seasons")
abline(fit)

par(mfrow=c(1,3))
acf(rowMeans(mydata),main="Data ACF")
acf(diff(rowMeans(mydata)),main="Differenced Data ACF")
acf(resid(fit),main="Detrended ACF")

plot(years[-1][-357][-356],diff(rowMeans(mydata)),type="l",xlab="Year",ylab="Degrees Celcius",main="Differenced Temperature Data")

#Kernel Smoother
plot(years,rowMeans(mydata,na.rm=TRUE),type="l",xlab="Year",ylab="Degrees Celcius",main="Temperature Data: Kernel Smoother")
lines(ksmooth(years, rowMeans(mydata), "normal", bandwidth=5), lwd=2, col=4)

#Lowess Smoother
plot(years,rowMeans(mydata,na.rm=TRUE),type="l",xlab="Year",ylab="Degrees Celcius",main="Temperature Data: Lowess Smoother")
lines(lowess(years,rowMeans(mydata),f=1/20),lwd=2,col=2)

diff1=diff(temp_data)
par(mfrow=c(1,3))
plot(diff1,type="l",main="Series diff1",ylab="Delta Degrees Celcius")
acf(diff1,lag="100")
pacf(diff1,lag="100")

diff2=diff(overall,differences=2)
par(mfrow=c(1,3))
plot(diff2,type="l")
acf(diff2,lag="100")
pacf(diff2,lag="100")

#Fitting an MA(2) model
model1=sarima(temp_data,0,1,2)


#Unit Root Testing
adf.test(diff1, k=0) # DF test
adf.test(diff1) # ADF test
pp.test(diff1)

phi0=-0.006258
phi1=-0.6610
phi2=-0.2941 

mysim=c(diff1[1],diff1[2])

for(i in (1:355))
{
  val=phi0+phi1*mysim[i-1]+phi2*mysim[i-2]
  mysim=c(mysim,val)
}


mysim=arima.sim(n=357,model=list(order=c(2,1,0),ar=c(phi1,phi2),sd=0))+phi0
plot(years,mysim,type="l")

plot(years[-1][-357],overall[-1][-357],type="l",xlab="Year",ylab="Degrees Celcius",main="Temperature Data with AR(1)")
lines(years,mysim,col=2)


#Trying to get this to work

fit=Arima(overall,order=c(2,1,0))
plot(years[-1][-357],overall[-1][-357],type="l",xlab="Year",ylab="Degrees Celcius",main="Temperature Data with AR(2) Fit")
lines(years[-1][-357],fitted(fit)[-1][-357],col="red")


#Prediction for Model 1
model1=arima(overall,order=c(0,1,2))
last_years=tail(c(years[-1],(2017:2030)),30)
last_temps=tail(c(overall[-1][-357],rep(NaN,13)),30)
plot(last_years,last_temps,type="l",ylim=c(4,12),xlab="Year",ylab="Degrees Celcius",main="Temperature Data Forecast to 2030")
forecast1=predict(model1,n.ahead=13)

adj=c(rep(NaN,16),last_temps[17],forecast1$pred)

lines(last_years,adj,col="red",type="o")

adj1=c(rep(NaN,17),forecast1$pred+forecast1$se)
adj2=c(rep(NaN,17),forecast1$pred-forecast1$se)

lines(last_years,adj1, lty="dashed", col=4)
lines(last_years,adj2, lty="dashed", col=4)

#Prediction for SEASONS
acf(diff(wint))
pacf(diff(wint))

acf(diff(spring))
pacf(diff(spring))

acf(diff(summ))
pacf(diff(summ))

acf(diff(fall))
pacf(diff(fall))

model1=arima(overall,order=c(2,1,0))
last_years=tail(c(years[-1],(2017:2030)),30)
last_temps=tail(c(overall[-1][-357],rep(NaN,13)),30)
plot(last_years,last_temps,type="l",ylim=c(4,12),xlab="Year",ylab="Degrees Celcius",main="Temperature Data Forecast to 2030")
forecast1=predict(model1,n.ahead=13)

adj=c(rep(NaN,16),last_temps[17],forecast1$pred)

lines(last_years,adj,col="red",type="o")

adj1=c(rep(NaN,17),forecast1$pred+fore$se)
adj2=c(rep(NaN,17),forecast1$pred-fore$se)

lines(last_years,adj1, lty="dashed", col=4)
lines(last_years,adj2, lty="dashed", col=4)



modelD=arima(diff1,order=c(2,0,0))
forecastD=predict(modelD,n.ahead=50)
plot(diff1,type="l",xlim=c(0,400))
lines(forecastD$pred,col="red")
lines(forecastD$pred+fore$se, lty="dashed", col=4)
lines(forecastD$pred-fore$se, lty="dashed", col=4)

#Periodogram for all data
alldata=c(t(mydata))[-1][-1431][-1430][-1429]

mvspec(alldata,log="yes")

myperiod=mvspec(altered,log="no")



par(mfrow=c(2,2))
mvspec(wint,log="no")
mvspec(spring,log="no")
mvspec(summ,log="no")
mvspec(fall,log="no")

#Nonparametric Spectral Techniques

k=kernel("daniell",5)
mvspec(temp_data,k,log="no")

#spliced analysis

middle=tail(head(overall,1925-1659),1925-1776)
myperiod=mvspec(middle,k,klog="no")

myspecs=myperiod$spec
for(i in (1:180))
{
  if(myspecs[i]>1.5)
  {
    print(i)
  }
}

par(mfrow=c(2,2))

k=kernel("daniell",9)

mvspec(wint,k,log="no")
mvspec(spring,k,log="no")
mvspec(summ,k,log="no")
mvspec(fall,k,log="no")

k=kernel("daniell",9)

plot(spec.ar(wint,log="no"))
plot(spec.ar(spring,log="no"))
plot(spec.ar(summ,log="no"))
plot(spec.ar(fall,log="no"))


#coherency between the seasons

seasons=cbind(wint,spring,summ,fall)

par(mfrow=c(2,3))

    sr = mvspec(cbind(wint,spring), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(18+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Spring")
    abline(h = C)    

    sr = mvspec(cbind(wint,summ), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(18+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Summer")
    abline(h = C)  
    
    sr = mvspec(cbind(wint,fall), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(18+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Fall")
    abline(h = C)  
    
    sr = mvspec(cbind(summ,spring), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(18+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Summer vs. Spring")
    abline(h = C)  
    
    
    sr = mvspec(cbind(summ,fall), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(8+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Summer vs. Fall")
    abline(h = C)  
    
    
    sr = mvspec(cbind(spring,fall), kernel("daniell",9), plot=FALSE)
    sr$df # df = 35.8625
    f = qf(.999, 2, sr$df-2) # = 8.529792
    C = f/(8+f) # = 0.321517
    plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Spring vs. Fall")
    abline(h = C)  

    
  #Chapter 5 Techniques
    logtemp = log(temp_data)-mean(log(temp_data))
    temp.fd = fracdiff(logtemp, nar=0, nma=2)
    temp.fd$d # = 0.3841688
    temp.fd$stderror.dpq # = 4.589514e-06 (questionable result!!)
    p = rep(1,31)
    for (k in 1:30){ p[k+1] = (k-temp.fd$d)*p[k]/(k+1) }
    plot(1:30, p[-1], ylab=expression(pi(d)), xlab="Index", type="h")
    res.fd = diffseries(log(temp_data), temp.fd$d) # frac diff resids
    res.arima = resid(arima(log(temp_data), order=c(1,1,1))) # arima resids
    par(mfrow=c(2,1))
    acf(res.fd, 100, xlim=c(4,97), ylim=c(-.2,.2), main="Fractional Differencing ACF")
    
    
    
#FIRST EVALUATION

training_data=head(temp_data,length(temp_data)-66)
y=head(years[-1][-357],290)
    
summary(model2<-lm(training_data~y+I(y^2)+I(y^3)))
summary(quadratic_test<-lm(training_data~y+I(y^2)+I(y^3)+I(y^4)))

plot(y,training_data,type="l",xlab="Year",ylab="Degrees Celcius",main="Training Data with Cubic Trend")
lines(y,model2$fitted.values,col="blue")


plot(y,model2$residuals,xlab="Year",ylab="Degrees Celcius",main="Detrended Training Data",type="l")
abline(h=mean(model2$residuals))


Box.test(model2$residuals,type="Ljung-Box")

acf(detrended)
pacf(detrended)

detrended_model=model2

bptest(detrended_model)

detrended=model2$residuals

sarima(detrended,0,0,5)
model3=arima(detrended,order=c(0,0,5))

extra_poly=forecast(model2$fitted.values,h=2016-1949)$mean

plot(c(tail(y,99),(1950:2016)),c(tail(training_data,99),rep(NaN,67)),ylim=c(7.5,11.5),xlab="Year",ylab="Degrees Celcius",main="Predicted against Actual",type="l")
#lines((1850:2016),c(tail(model2$fitted,100),extra_poly),col="green")
y1=years[-1][-357]
summary(model<-lm(temp_data~y1+I(y1^2)+I(y1^3)))

#high:
highs=forecast(model3,h=67)$upper[,1]+extra_poly
xtemp=c(1949:2016)
lines(x=xtemp,y=c(training_data[length(training_data)],highs),col="blue",type="l",lty=3)

#median:
meds=forecast(model3,h=67)$mean+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],meds),col="red",type="l")

#low:
lows=forecast(model3,h=67)$lower[,1]+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],lows),col="blue",type="l",lty=3)


#actuals
lines(x=xtemp,y=c(training_data[length(training_data)],tail(overall,67)),col="purple",type="l",font.main=4)


#actual_mean
lines(tail(y1,67),tail(model_test$fitted.values,67),col="green",lwd=2)


mypoly <- function(x){

    return(tail(forecast(model2$fitted.values,h=x-1949)$mean,1))
  #return(-1.377e+03  +x*2.310e+00+x^2*-1.283e-03 +x^3*2.373e-07)
}
  



#SECOND EVALUATION
training_data=head(temp_data,length(temp_data)-16)
y=head(years[-1][-357],340)

summary(model2<-lm(training_data~y+I(y^2)+I(y^3)))
summary(quadratic_test<-lm(training_data~y+I(y^2)+I(y^3)+I(y^4)))

plot(y,training_data,type="l",xlab="Year",ylab="Degrees Celcius",main="Training Data with Cubic Trend")
lines(y,model2$fitted.values,col="blue")

plot(y,model2$residuals,xlab="Year",ylab="Degrees Celcius",main="Detrended Training Data",type="l")
abline(h=mean(model2$residuals))

Box.test(model2$residuals,type="Ljung-Box")
detrended=model2$residuals
par(mfrow=c(1,2))
acf(detrended)
pacf(detrended)
par(mfrow=c(1,1))
detrended_model=model2

bptest(detrended_model)

detrended=model2$residuals

sarima(detrended,3,0,0)
model3=arima(detrended,order=c(3,0,0))

extra_poly=forecast(model2$fitted.values,h=2016-1999)$mean

plot(c(tail(y,99),(2000:2016)),c(tail(training_data,99),rep(NaN,17)),ylim=c(7.5,11.5),xlab="Year",ylab="Degrees Celcius",main="Predicted against Actual 2000+",type="l")
#lines((1900:2016),c(tail(model2$fitted,100),extra_poly),col="orange")
y1=years[-1][-357]
summary(model<-lm(temp_data~y1+I(y1^2)+I(y1^3)))

#high:
highs=forecast(model3,h=17)$upper[,1]+extra_poly
xtemp=c(1999:2016)
lines(x=xtemp,y=c(training_data[length(training_data)],highs),col="blue",type="l",lty=3)

#median:
meds=forecast(model3,h=17)$mean+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],meds),col="red",type="l")

#low:
lows=forecast(model3,h=17)$lower[,1]+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],lows),col="blue",type="l",lty=3)


#actuals
lines(x=xtemp,y=c(training_data[length(training_data)],tail(overall,17)),col="purple",type="l",font.main=4)


#actual_mean
lines(tail(y1,17),tail(model_test$fitted.values,17),col="green",lwd=2)


mypoly <- function(x){
  
  return(tail(forecast(model2$fitted.values,h=x-1949)$mean,1))
  #return(-1.377e+03  +x*2.310e+00+x^2*-1.283e-03 +x^3*2.373e-07)
}






#Third evaluation

training_data=head(temp_data,length(temp_data)-116)
y=head(years[-1][-357],240)

summary(model2<-lm(training_data~y+I(y^2)+I(y^3)))
summary(quadratic_test<-lm(training_data~y+I(y^2)+I(y^3)+I(y^4)))

plot(y,training_data,type="l",xlab="Year",ylab="Degrees Celcius",main="Training Data (up to 1900) with Cubic Trend")
lines(y,model2$fitted.values,col="blue")

plot(y,model2$residuals,xlab="Year",ylab="Degrees Celcius",main="Detrended Training Data",type="l")
abline(h=mean(model2$residuals))

Box.test(model2$residuals,type="Ljung-Box")
detrended=model2$residuals
par(mfrow=c(1,2))
acf(detrended)
pacf(detrended)
par(mfrow=c(1,1))

detrended_model=model2

bptest(detrended_model)

detrended=model2$residuals

sarima(detrended,0,0,5)
model3=arima(detrended,order=c(0,0,5))

extra_poly=forecast(model2$fitted.values,h=2016-1899)$mean

plot(c(tail(y,99),(1900:2016)),c(tail(training_data,99),rep(NaN,117)),ylim=c(7.5,11.5),xlab="Year",ylab="Degrees Celcius",main="Predicted against Actual 1900+",type="l")
#lines((1900:2016),c(tail(model2$fitted,100),extra_poly),col="orange")
y1=years[-1][-357]
summary(model<-lm(temp_data~y1+I(y1^2)+I(y1^3)))

#high:
highs=forecast(model3,h=117)$upper[,1]+extra_poly
xtemp=c(1899:2016)
lines(x=xtemp,y=c(training_data[length(training_data)],highs),col="blue",type="l",lty=3)

#median:
meds=forecast(model3,h=117)$mean+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],meds),col="red",type="l")

#low:
lows=forecast(model3,h=117)$lower[,1]+extra_poly
lines(x=xtemp,y=c(training_data[length(training_data)],lows),col="blue",type="l",lty=3)


#actuals
lines(x=xtemp,y=c(training_data[length(training_data)],tail(overall,117)),col="purple",type="l",font.main=4)


#actual_mean
lines(tail(y1,117),tail(model_test$fitted.values,117),col="green",lwd=2)


mypoly <- function(x){
  
  return(tail(forecast(model2$fitted.values,h=x-1949)$mean,1))
  #return(-1.377e+03  +x*2.310e+00+x^2*-1.283e-03 +x^3*2.373e-07)
}


#Newest Periodic Stuff
y=(1660:2015)
summary(model1<-lm(temp_data~y+I(y^2)+I(y^3)))


plot(y,temp_data,type="l",ylab="Degrees Celcius",xlab="Years")
lines(y,model1$fitted.values,col="blue",lwd=2)

plot(y,model1$residuals,type="l",ylab="Degrees Celcius",xlab="Years")
abline(lwd=2,h=0,col="blue")

detrended=model1$residuals

adf.test(detrended, k=0) # DF test
adf.test(detrended) # ADF test
pp.test(detrended)

k = kernel("daniell", 6)
mvspec(detrended, k, log="no")

#season trends, season detrended, season periods, season smoothed
y=(1660:2016)
summary(modelw<-lm(wint~y))

y=(1659:2015)
summary(modelr<-lm(spring~y+I(y^2)+I(y^3)))

y=(1660:2016)
summary(models<-lm(summer~y+I(y^2)+I(y^3)))

y=(1660:2016)
summary(modelf<-lm(fall~y+I(y^2)+I(y^3)))

par(mfrow=c(2,2))

plot(y,winter,type="l",ylab="Degrees Celcius",xlab="Years",main="Winter")
lines(y,modelw$fitted.values,col="blue",lwd=2)

plot(y,spring,type="l",ylab="Degrees Celcius",xlab="Years",main="Spring")
lines(y,modelr$fitted.values,col="blue",lwd=2)

plot(y,summer,type="l",ylab="Degrees Celcius",xlab="Years",main="Summer")
lines(y,models$fitted.values,col="blue",lwd=2)

plot(y,fall,type="l",ylab="Degrees Celcius",xlab="Years",main="Fall")
lines(y,modelf$fitted.values,col="blue",lwd=2)
  
par(mfrow=c(2,2))

plot(years,modelw$residuals,type="l",ylab="Degrees Celcius",main="Winter")
abline(lwd=2,h=0,col="blue")

plot(years,modelr$residuals,type="l",ylab="Degrees Celcius",main="Spring")
abline(lwd=2,h=0,col="blue")

plot(years,models$residuals,type="l",ylab="Degrees Celcius",main="Summer")
abline(lwd=2,h=0,col="blue")

plot(years,modelf$residuals,type="l",ylab="Degrees Celcius",main="Fall")
abline(lwd=2,h=0,col="blue")

winter_detrended=modelw$residuals
spring_detrended=modelr$residuals
summer_detrended=models$residuals
fall_detrended=modelf$residuals


mvspec(winter_detrended,log="no")
mvspec(spring_detrended,log="no")
mvspec(summer_detrended,log="no")
mvspec(fall_detrended,log="no")

k=kernel()k = kernel("daniell", 6)
mvspec(winter_detrended,log="no",kernel=k)
mvspec(spring_detrended,log="no",kernel=k)
mvspec(summer_detrended,log="no",kernel=k)
mvspec(fall_detrended,log="no",kernel=k)
  




#new coherency

#coherency between the seasons

seasons=cbind(wint,spring,summ,fall)

par(mfrow=c(2,3))

sr = mvspec(cbind(winter_detrended,spring_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(18+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Spring")
abline(h = C)    

sr = mvspec(cbind(winter_detrended,summer_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(18+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Summer")
abline(h = C)  

sr = mvspec(cbind(winter_detrended,fall_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(18+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Winter vs. Fall")
abline(h = C)  

sr = mvspec(cbind(summer_detrended,spring_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(18+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Summer vs. Spring")
abline(h = C)  


sr = mvspec(cbind(summer_detrended,fall_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(8+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Summer vs. Fall")
abline(h = C)  


sr = mvspec(cbind(spring_detrended,fall_detrended), kernel("daniell",9), plot=FALSE)
sr$df # df = 35.8625
f = qf(.999, 2, sr$df-2) # = 8.529792
C = f/(8+f) # = 0.321517
plot(sr, plot.type = "coh", ci.lty = 2,main="Coh: Spring vs. Fall")
abline(h = C)  


par(mfrow=c(1,1))
(spec.ar(detrended,log="no"))
par(mfrow=c(2,2))
(spec.ar(winter_detrended))
(spec.ar(spring_detrended))
(spec.ar(summer_detrended,log="no"))
(spec.ar(fall_detrended,log="no"))