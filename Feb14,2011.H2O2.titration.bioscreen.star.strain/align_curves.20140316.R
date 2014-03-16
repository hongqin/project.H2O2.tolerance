#20140316 align growth curves by minimizing sum of errors 

rm(list=ls())
debug = 9; 

#
pairwise_sum_of_errors = function( X1, X2, Start1,End1,Start2, End2 ) {
  if( (End1 - Start1) == (End2 - Start2) ) {
    return( sum( ( X1[Start1:End1]-X2[Start2:End2] )^2))
  } else {
    print("Error: X1 and X2 should have the same lengths")
    return( NA); 
  }
}

#testing 
x1 = 1:100;  x2 = 5:104
pairwise_sum_of_errors( x1, x2, 5, 100, 1,96)
sum((x1[5:100]-x2[1:96])^2)
pairwise_sum_of_errors( x1, x2, 5, 100, 1,100)

#align two vectors with equal lengths
#Let's assume start2 lags behind start1 
#start2 = start1 + delta # 2nd start lags behind 1st start by delta
pairwise_cost_of_shifted_vectors = function( delta, X1, X2,  low.threshold=0.01, debug=10 ){
  delta = floor( (delta+0.5) ) #delta must be an integer
  X1[X1<low.threshold]=0;
  X2[X2<low.threshold]=0;
  Start1 = 1;
  Start2 = Start1 + delta;  #left truncation at delta
  End2 = length(X2)
  End1 = length(X1) - delta; #right truncation at delta
  if(debug){ print(paste('Start1', Start1, "Start2", Start2, "End1", End1, "End2", End2))
  }
  if( (End1 - Start1) == (End2 - Start2) ) {
    return( sum( ( X1[Start1:End1]-X2[Start2:End2] )^2))
  } else {
    print("Error: X1 and X2 should have the same lengths")
    return( NA); 
  }
}

#testing 1
x1 = 5:104; x2 = 1:100;  
pairwise_sum_of_errors( x1, x2, 1,96, 5,100 )
pairwise_cost_of_shifted_vectors( 4, x1, x2)

delta = 4
res1 = optim(c(delta), fn=pairwise_cost_of_shifted_vectors, X1=x1, X2=x2, method="Brent", lower=c(1), upper=c(100))
#return 3 not 4? 

#testing 2
x1 = 50:149; x2 = 1:100;  
pairwise_sum_of_errors( x1, x2, 1,50, 50,99 )
pairwise_cost_of_shifted_vectors( 49, x1, x2)

delta = 50
res2 = optim(c(delta), fn=pairwise_cost_of_shifted_vectors, X1=x1, X2=x2, method="Brent", lower=c(1), upper=c(100))
#return 50 not 49?

### real data
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
mycolors = c("black", "cyan", "yellow","green", "brown", 
             "red","blue","purple","pink", "orange", "navy" )

myH2O2 = unique(tbw[,3], na.rm=T)
myH2O2 = myH2O2[! is.na(myH2O2)] / 100
myDilutions = unique( tbw[,5])

#tbslope = tb; 
#tbw$t.max.slope = NA; #time of max slope
tb2$delta = NA

myxlim = c(0, max(tb$Hrs))
myylim = c(0, 1.8)
g= 1
#pdf("bioscreen.H2O2.Feb8,2011.logY.pdf", width=8, height=8)
#for ( g in 0:groups){
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
