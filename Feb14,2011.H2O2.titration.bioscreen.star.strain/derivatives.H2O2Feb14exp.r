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

tbw =	read.csv("wells.Feb15.H2O2.bioscreen.csv", colClasses=c(NA,"character", NA, NA, NA));
tbw$wells2 = paste("Wells", tbw$Wells, sep='.')

window = 10;
groups = 13;
mycolors = c("black", "cyan", "yellow","green", "brown", "red","blue","purple","pink", "orange", "navy" )

myH2O2 = unique(tbw[,3], na.rm=T)
myH2O2 = myH2O2[! is.na(myH2O2)] / 100
myDilutions = unique( tbw[,5])

tbslope = tb; 
tbw$t.max.slope = NA; #time of max slope

myxlim = c(0, max(tb$Hrs))
myylim = c(0, 1.8)
pdf("bioscreen.H2O2.Feb8,2011.logY.pdf", width=8, height=8)
for ( g in 0:groups){
	#first plot wildtype
   mywell = names(tb)[g*window + 2]
   mystrain = tbw$Strain[g*window + 1]
   #plot( tb[ ,g*window + 2] ~ tb$Hrs, main=paste('group', g, mywell, mystrain))
   plot( tb[ ,g*window + 2] ~ tb$Hrs, main=paste( mywell, mystrain), type='l',xlab='Hours', ylab='OD600nm', las=2, xlim=myxlim, log='y' )
      
   res = find.slope( tb$Hrs, tb[ ,g*window + 2], window )
   tbslope[ ,  g*window + 2 ] = c( res$yout, rep(NA, length(tb[,1]) - length(res$yout)))   
   tbslope$time = c( res$xout, rep(NA, length(tb[,1]) - length(res$yout))) 
   
   lo = lowess( res$xout, res$yout, f=1/10 )
	tbw$t.max.slope[g*window + 1] = lo$x[lo$y==max(lo$y, na.rm=T)][1]
	
	for( i in 2:(window-1)) {
		lines( tb[ ,g*window + i + 1] ~ tb$Hrs, col=mycolors[i]);
 
      	res = find.slope( tb$Hrs, tb[ ,g*window + i + 1], window )
      	tbslope[ , g*window + i + 1] = c( res$yout, rep(NA, length(tb[,1]) - length(res$yout)))   
   		lo = lowess( res$xout, res$yout, f=1/10 )
		tbw$t.max.slope[g*window + i] = lo$x[lo$y==max(lo$y, na.rm=T)][1]
 	}
	
	if ( g < 12) {
 	 legend( 30 , 1 ,legend=myH2O2, col=mycolors, lty=1);
 	} else {
 	 legend( 30 , 1 ,legend=myDilutions, col=mycolors, lty=1);
 	}
	
	my.ylim = c(0, max(tbslope[,2]*2.0 ))
	
	plot( tbslope[ ,g*window + 2] ~ tbslope$time, type='l',main=paste( mywell, mystrain, 'slope'), xlab='Hours', ylab='slope', las=2,xlim=myxlim)
	for( i in 2:(window-1)) {
		lines( tbslope[ ,g*window + i + 1] ~ tbslope$time, col=mycolors[i]);
  	}

   myy = max( tbslope[ ,g*window + 2], na.rm=T )
	if ( g < 12) {
 	 legend( 30 , myy ,legend=myH2O2, col=mycolors, lty=1);
 	} else {
 	 legend( 30 , myy ,legend=myDilutions, col=mycolors, lty=1);
 	}

   sub = tbw[ (g*window+1) : (g*window+9)  , ]
	if ( g < 12) {
 	 	sub$H2O2 = sub[,3]/100
   	sub$H2O2[sub$H2O2==0] = 0.001
 		plot( sub$t.max.slope ~ sub$H2O2, log='x', type='b', main=paste(mystrain), xlab='H2O2', ylab='t.max.slope' )
 	} else {
 		plot( sub$t.max.slope ~ sub$Dilutions, log='x', type='b', main=paste(mystrain), xlab='Dilutions', ylab='t.max.slope' )
 	}			
}
dev.off()

write.csv(tbslope, "_slope.H2O2Feb14Exp.032111.csv");
write.csv(tbw, "_out.summary.H2O2Feb14Exp.032111.csv");

quit("yes")
