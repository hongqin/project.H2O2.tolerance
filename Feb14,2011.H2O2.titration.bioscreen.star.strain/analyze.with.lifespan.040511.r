lifespan = read.delim("062705.rls.cls.tab")
str(lifespan)

halo = read.csv("H2O2.Halo.star.strains.by.DP.ADC.csv")
halo = halo[, 1:8]
halo
halo30 = halo[halo$H2O2Concentration==0.3, ]

pos = match( halo30$Strain, lifespan$strain)
halo30 = cbind( halo30, lifespan[pos,])
rownames(halo30) = as.character(halo30[,1])
sub = halo30[,c(4,5,10:13)]
sub

#remove M14, M22
sub["M14", 1:2] = c(NA, NA)
sub["M22", 1:2] = c(NA, NA)

pairs(sub)

summary(lm( sub$DiameterOuter ~ sub$rls ))
summary(lm( sub$DiameterOuter ~ sub$cls ))
summary(lm( sub$DiameterOuter ~ sub$a ))
summary(lm( sub$DiameterOuter ~ sub$b ))
summary(lm( sub$DiameterOuter ~ sub$rls + sub$cls ))
summary(lm( sub$DiameterOuter ~ sub$a + sub$b + sub$cls ))

sub$longCLS = 0;
sub$longCLS[sub$cls> 6] = 1;
summary(lm( sub$DiameterOut ~ sub$longCLS))
longCLS = sub[sub$longCLS ==1, ]
shortCLS = sub[sub$longCLS ==0, ]
#wilcox.test( longCLS$DiameterOuter, shortCLS$DiameterOuter)
t.test( longCLS$DiameterOuter, shortCLS$DiameterOuter, alter="less")
t.test( longCLS$DiameterInner, shortCLS$DiameterInner, alter="less")

summary(lm( sub$DiameterInner ~ sub$rls ))
summary(lm( sub$DiameterInner ~ sub$cls ))
summary(lm( sub$DiameterInner ~ sub$a ))
summary(lm( sub$DiameterInner ~ sub$b ))

plot( DiameterOuter ~ cls, data=sub, pch=16 )
text(12, 2.2, "p=0.127 t-test")

bc = read.csv("Doubling.H2O2.concentration.040511.csv")
pos2 = match(bc$strain, halo30$strain)
bc = cbind( bc, halo30[pos2, ])

plot( bc$DoubleConc ~ bc$DiameterOuter, log='y' )
summary(lm( bc$DoubleConc ~ bc$DiameterOuter  ))
summary(lm( log( bc$DoubleConc) ~ bc$DiameterOuter  ))

summary(lm( bc$DoubleConc ~ bc$rls  ))
summary(lm( bc$DoubleConc ~ bc$cls  ))
summary(lm( log( bc$DoubleConc ) ~ bc$cls  ))

bc2 = bc[,c(1,7,8,12:16 )]
bc2
#remove M1-2
bc2[1,1] = NA
pairs( bc2 )

summary(lm( bc2$DoubleConc ~ bc2$rls  ))
summary(lm( bc2$DoubleConc ~ bc2$cls  ))
#summary(lm( log( bc2$DoubleConc ) ~ bc2$cls  ))

plot( DoubleConc ~ cls, data=bc2, pch=16 )
text(12, 1.75, "p=0.102")


