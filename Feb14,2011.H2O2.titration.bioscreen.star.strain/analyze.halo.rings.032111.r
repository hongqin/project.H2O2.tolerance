
rm(list=ls())
tb0 = read.csv( "H2O2.Halo.star.strains.by.DP.ADC.csv" );
h30 = hist(tb0$DiameterOuter[tb0$H2O2Concentration==0.3], br=12)
h15 = hist(tb0$DiameterOuter[tb0$H2O2Concentration==0.15], br= h30$breaks)
   
 #generate the comparison table
 tb <-  data.frame( rbind(h30$counts, h15$counts) )  ;
 names( tb ) <- h30$mids;
 row.names(tb) <- c( "30%", "15%" )
 tb
 
 pdf( "hist.haloRings.barplot.bw.pdf", width=7, height=7 );
 barplot( as.matrix(tb), beside=T, col=c('black','gray'), 
	ylab="Frequency", xlab="Diamters(cm)",
    	legend= c( "30%", "15%" ),
 );
 title(main="" )
 dev.off()
 
quit("no");
