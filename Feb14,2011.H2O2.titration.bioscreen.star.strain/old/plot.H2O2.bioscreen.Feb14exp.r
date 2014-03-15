debug = 0; 

find.slope = function( xin, yin, window) {
 halfwindow = floor( window /2 )
 end = length(xin) - halfwindow
 xout = (halfwindow+1):end;
 yout = (halfwindow+1):end;
 count = 1
 for ( i in (halfwindow+1):end ) {
   pos = (i-halfwindow):(i+halfwindow)
   m = lm( yin[pos] ~ xin[pos] )
   xout[count]  = xin[i]
   yout[count] = m$coeff[2]
   count = count + 1;
 }
 return(data.frame(cbind(xout,yout)) )
}

tb = read.csv( "Feb15,2011.H2O2.titration.12.star.strain.csv");
tb = tb[!is.na(tb[,1]), ] #remove empty spots

#change time to numerica values
tb$Time = as.character(tb$Time)
bb = strsplit(tb$Time, ":")
tb$Hrs <- sapply(bb,function(x){as.numeric(x[1])+as.numeric(x[2])/60 + as.numeric(x[3])/3600 })

mycolors = c("black", "cyan", "yellow","green", "brown", "red","blue","purple","pink", "orange", "navy" )

plot( tb[,2] ~ tb$Hrs, main="orange is 6, red=7")
lines( tb[,3] ~ tb$Hrs, col="blue")
lines( tb[,4] ~ tb$Hrs, col='purple')
lines( tb[,5] ~ tb$Hrs, col='green')
lines( tb[,6] ~ tb$Hrs, col='orange')
lines( tb[,7] ~ tb$Hrs, col='red')

tbw =	read.csv("wells.Feb15.H2O2.bioscreen.csv", colClasses=c(NA,"character", NA, NA, NA));
tbw$wells2 = paste("Wells", tbw$Wells, sep='.')

window = 10;
groups = 13;
#mycolors = rainbow(9)

myH2O2 = unique(tbw[,3], na.rm=T)
myH2O2 = myH2O2[! is.na(myH2O2)] / 100 #the stock is 100X
myDilutions = unique( tbw[,5])

pdf("bioscreen.H2O2.Feb8,2011.pdf", width=8, height=8)
for ( g in 0:groups){
	#first plot wildtype
   mywell = names(tb)[g*window + 2]
   mystrain = tbw$Strain[g*window + 1]
   #plot( tb[ ,g*window + 2] ~ tb$Hrs, main=paste('group', g, mywell, mystrain))
   plot( tb[ ,g*window + 2] ~ tb$Hrs, main=paste( mywell, mystrain), xlab='Hours', ylab='OD600nm', las=2)
   
	for( i in 2:(window-2)) {
		lines( tb[ ,g*window + i+1] ~ tb$Hrs, col=mycolors[i]); # 032011 change
	}
	
	if ( g < 12) {
 	 legend( 30 , 1 ,legend=myH2O2, col=mycolors, lty=1);
 	} else {
 	 legend( 30 , 1 ,legend=myDilutions, col=mycolors, lty=1);
 	}
	
}
dev.off()






quit("yes")

wells  =	read.csv("wells.Feb9,2011.csv") #, colClasses=c("character","character"));
wells[,c(1,4)] = wells[,c(1,4)]+100
wells$wells  = paste("Well", wells[,1], sep='.')
wells$wells2 = paste("Well", wells[,4], sep='.')

tb$Time = as.character(tb$Time)
bb = strsplit(tb$Time, ":")
tb$Hrs <- sapply(bb,function(x){as.numeric(x[1])+as.numeric(x[2])/60 + as.numeric(x[3])/3600 })

base = tb[1, 2:31]
for( j in 2:31) {
  base[1,j-1] = min(tb[1:30,j])
}

#for( i in 1:length(tb[,2]) ) {  #this does not affect growth curve analysis
#  tb[i,2:31] = tb[i,2:31] - base
#}

tb = tb[, -c(14:20)] #remove empty wells

#for debug, are the derivaties calculated correctly? 
if( debug > 1) {
 tb[,2] = tb$Hrs*0.5 + 0.1  #linear model
 tb[,3] = tb$Hrs^2* + 0.1 #quadratic
}


#calculate the growth rate (derivative) in sliding window
window = 31;
halfwindow = floor( window/2 )
end = length(tb[,1]) - halfwindow 

out = data.frame(names(tb)[-c(1,length(tb[1,]))])
out$t.max.slope = NA;
out$t.max.slope.lowess = NA;

for( j in 2:(length(tb[1,])-1) ) {
#for( j in 2:5){
 res = find.slope( tb$Hrs, tb[,j], window) 
 out$t.max.slope[j-1] = res$xout[res$yout==max(res$yout)][1] #too much noise
 time  = res$xout
 slope = res$yout
 lo = lowess( res$xout, res$yout, f=1/10 )
 #res2 = data.frame(cbind(lo$x, lo$y))
 #out$t.max.slope.lowess[j-1] = res2$x[res$yout==max(res2$y)][1] #will this introduce bias
 out$t.max.slope.lowess[j-1] = lo$x[lo$y==max(lo$y)][1] #will this introduce bias
 
 ymax = max(2, max(tb[,j]));
 plot( tb[,j] ~ tb$Hrs, ylab = names(tb)[j], ylim=c(-0.1, ymax), type='l', xlim = c(0, 35) ); 
 par(new=T)
 plot( slope ~ time, col="red", type='l', xlim=c(0,35), ylim=c( min(slope, na.rm=T)-0.05, max(slope,na.rm=T)+0.05 ), axes=F, ylab='', xlab='' )
 lines( lo$y ~ lo$x, col="blue", lwd=2)  
}

rownames( out ) = out[,1]

#for standard curves
wells2 = intersect( wells$wells2, out[,1])
standard.tb = out[wells2, ]
standard.tb$cell.dilutions = wells$cell.dilutions[ match(wells2, wells$wells2)  ]
standard.tb
plot( standard.tb$ t.max.slope  ~ standard.tb$cell.dilutions, log='x') # not very good
plot( standard.tb$ t.max.slope.lowess ~ standard.tb$cell.dilutions, log='x') #Good straight line!
standard.tb[1,] = NA;
standard.tb
m = lm( t.max.slope.lowess ~ log10(cell.dilutions), data=standard.tb)
summary(m)

#for rapamycin treatment
rap.tb = out[1:12, ]
rap.tb$rapamycin = wells$Rapamycin.ng.mL.[match(rap.tb[,1], wells[,1])]
rap.tb
plot( t.max.slope.lowess ~ rapamycin, data=rap.tb)
plot( t.max.slope.lowess ~ rapamycin, data=rap.tb, log='x')

quit("yes")

strains = unique(wells[,1])


mycols = rainbow( length(wells[,1]))

j =1; 
plot( tb[, wells$well[j]] ~ tb$hours, type='l', col="white", ylim = c(0,1.8), xlab="hours", ylab="OD600"
,xlim=c(0,50)
#, log='x' 
)
for( j in 1:length(wells[,1])) {
 lines( tb[, wells$well[j]] ~ tb$hours, col=mycols[j]  );
}
legend( 40, 1.5, wells[,1], col=mycols, lwd=1)

for ( mystrain in strains) {
  wells2 = wells[ wells[,1]==mystrain, ]
  mycols = c("red","green","blue")

  plot( tb[, wells2$well[1]] ~ tb$hours, type='l', col="white", ylim = c(0,1.6), 
      xlab="hours", ylab="OD600"
      #, log='x'
      ,xlim=c(0,20) 
      )
  for( j in 1:length(wells2[,1])) {
    lines( tb[, wells2$well[j]] ~ tb$hours, col=mycols[j]  );
  }
  # legend( 40, 1.5, wells[,1], col=mycols, lwd=1)
  title( mystrain )		
}
