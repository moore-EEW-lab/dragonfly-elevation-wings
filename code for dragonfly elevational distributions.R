setwd('/Users/michaelmoore/desktop/Working Directory')

library(car)
library(ggplot2)
library(lmerTest)
library(phytools)
library(nlme)
library(plyr)
library(MuMIn)
library(visreg)


drags.dat <- read.csv('elevation.dat.drags.csv')
rownames(drags.dat) <- drags.dat$binom
na.drags <- read.tree('phylo.drags.elev.tre')


#### mean elevation of highest grid cell
mod00a <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
mod00b <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
mod00c <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
mod00d <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(mod00a, mod00b, mod00c, mod00d)
       # df     AICc
# mod00a  6 971.4514
# mod00b  7 730.9328
# mod00c  7 730.9328
# mod00d  6 747.3871

summary(mod00b)
# Generalized least squares fit by REML
  # Model: log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC     BIC    logLik
  # 730.5519 756.408 -358.2759

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5175356 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)           1.584491 1.7019619  0.930979  0.3526
# log(body.size)        0.310054 0.2628460  1.179604  0.2391
# log(rel.wing)         5.913270 2.3390027  2.528116  0.0120
# log(range.size)       0.349136 0.0348243 10.025643  0.0000
# log(bio1_mean + 100) -0.860483 0.1494672 -5.757005  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.468                       
# log(rel.wing)         0.252  0.004                
# log(range.size)      -0.619 -0.139 -0.076         
# log(bio1_mean + 100) -0.582 -0.135 -0.134   0.313 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -2.91298725 -0.56359896 -0.01971429  0.55110349  1.80670477 

# Residual standard error: 1.058088 
# Degrees of freedom: 302 total; 297 residual

confint(mod00b)
                          # 2.5 %     97.5 %
# (Intercept)          -1.7512935  4.9202745
# log(body.size)       -0.2051144  0.8252231
# log(rel.wing)         1.3289090 10.4976312
# log(range.size)       0.2808812  0.4173898
# log(bio1_mean + 100) -1.1534335 -0.5675330

### elevational breadth - add 0.1 to make all response values positive
r00a <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
r00b <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
r00c <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
r00d <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(r00a, r00b, r00c, r00d)
     # df     AICc
# r00a  6 1302.182
# r00b  7 1039.040
# r00c  7 1039.040
# r00d  6 1039.706

summary(r00d)

# Generalized least squares fit by REML
  # Model: log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 1039.421 1061.583 -513.7104

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          -9.922173 2.4617346 -4.030562  0.0001
# log(body.size)        0.091335 0.2996671  0.304787  0.7607
# log(rel.wing)         7.818613 2.6118204  2.993549  0.0030
# log(range.size)       0.751498 0.0576144 13.043577  0.0000
# log(bio1_mean + 100) -0.582295 0.2230271 -2.610870  0.0095

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.326                       
# log(rel.wing)         0.343  0.308                
# log(range.size)      -0.772 -0.105 -0.386         
# log(bio1_mean + 100) -0.668 -0.100 -0.278   0.318 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.59253157 -0.43807609  0.06115794  0.58757180  1.86477648 

# Residual standard error: 1.324674 
# Degrees of freedom: 302 total; 297 residual

confint(r00d)
                           # 2.5 %     97.5 %
# (Intercept)          -14.7470844 -5.0972620
# log(body.size)        -0.4960019  0.6786714
# log(rel.wing)          2.6995395 12.9376873
# log(range.size)        0.6385759  0.8644203
# log(bio1_mean + 100)  -1.0194195 -0.1451696

### mean elevation of the lowest occupied grid cell - add 0.1 to make all response values positive
min00a <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
min00b <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
min00c <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
min00d <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(min00a, min00b, min00c, min00d)
       # df      AICc
# min00a  6 1181.9400
# min00b  7  976.8986
# min00c  7  976.8986
# min00d  6  979.6529

summary(min00b)
# Generalized least squares fit by REML
  # Model: log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 976.5177 1002.374 -481.2588

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5479602 

# Coefficients:
                         # Value Std.Error    t-value p-value
# (Intercept)          23.088388  2.587703   8.922347  0.0000
# log(body.size)        0.421935  0.400614   1.053221  0.2931
# log(rel.wing)         3.562826  3.547770   1.004244  0.3161
# log(range.size)      -0.571347  0.052621 -10.857720  0.0000
# log(bio1_mean + 100) -1.121007  0.226625  -4.946527  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.469                       
# log(rel.wing)         0.250  0.002                
# log(range.size)      -0.613 -0.140 -0.073         
# log(bio1_mean + 100) -0.578 -0.136 -0.132   0.312 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.7973431 -0.5566213 -0.1297642  0.4308258  2.0877056 

# Residual standard error: 1.645537 
# Degrees of freedom: 302 total; 297 residual

confint(min00b)
                          # 2.5 %     97.5 %
# (Intercept)          18.0165828 28.1601941
# log(body.size)       -0.3632542  1.2071246
# log(rel.wing)        -3.3906768 10.5163278
# log(range.size)      -0.6744832 -0.4682115
# log(bio1_mean + 100) -1.5651844 -0.6768302

## maximum elevation of highest occupied grid cell
max00a <- gls(log(max.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
max00b <- gls(log(max.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
max00c <- gls(log(max.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
max00d <- gls(log(max.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(max00a, max00b, max00c, max00d)
       # df     AICc
# max00a  6 904.3276
# max00b  7 672.7869
# max00c  7 672.7869
# max00d  6 694.2917

summary(max00b)

# Generalized least squares fit by REML
  # Model: log(max.elev) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
      # AIC      BIC   logLik
  # 672.406 698.2621 -329.203

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5290211 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)           4.157371 1.5460806  2.688974  0.0076
# log(body.size)        0.157943 0.2390058  0.660832  0.5092
# log(rel.wing)         4.745377 2.1229868  2.235236  0.0261
# log(range.size)       0.290615 0.0315624  9.207626  0.0000
# log(bio1_mean + 100) -0.824399 0.1356411 -6.077800  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.468                       
# log(rel.wing)         0.251  0.003                
# log(range.size)      -0.617 -0.139 -0.075         
# log(bio1_mean + 100) -0.581 -0.135 -0.133   0.313 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.7848834 -0.4778526 -0.1006900  0.4707425  1.5698318 

# Residual standard error: 0.9692417 
# Degrees of freedom: 302 total; 297 residual

confint(max00b)
                          # 2.5 %     97.5 %
# (Intercept)           1.1271084  7.1876331
# log(body.size)       -0.3105001  0.6263856
# log(rel.wing)         0.5843989  8.9063544
# log(range.size)       0.2287536  0.3524760
# log(bio1_mean + 100) -1.0902509 -0.5585477

## minimum elevation of lowest occupied grid cell - add 172 to make all response values positive
minny00a <- gls(log(min.elev + 172) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
minny00b <- gls(log(min.elev + 172) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
minny00c <- gls(log(min.elev + 172) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.75, form=~binom, phy = na.drags, fixed = FALSE))
minny00d <- gls(log(min.elev + 172) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(minny00a, minny00b, minny00c, minny00d) # pagel's lambda models aren't converging well. Fit nearly as well as star phylo, so go with star phylo
         # df      AICc
# minny00a  6 1214.2517
# minny00b  7  938.5150
# minny00c  7  937.4404
# minny00d  6  942.2159

summary(minny00d)
# Generalized least squares fit by REML
  # Model: log(min.elev + 172) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 941.9311 964.0935 -464.9656

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          19.974328 2.0891184  9.561128  0.0000
# log(body.size)        0.001563 0.2543085  0.006146  0.9951
# log(rel.wing)        -3.349333 2.2164867 -1.511100  0.1318
# log(range.size)      -0.478834 0.0488937 -9.793377  0.0000
# log(bio1_mean + 100) -0.553655 0.1892689 -2.925229  0.0037

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.326                       
# log(rel.wing)         0.343  0.308                
# log(range.size)      -0.772 -0.105 -0.386         
# log(bio1_mean + 100) -0.668 -0.100 -0.278   0.318 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -4.3310923 -0.2015900  0.2657207  0.5426477  1.2229318 

# Residual standard error: 1.124167 
# Degrees of freedom: 302 total; 297 residual

confint(minny00d)
                          # 2.5 %     97.5 %
# (Intercept)          15.8797312 24.0689248
# log(body.size)       -0.4968725  0.4999983
# log(rel.wing)        -7.6935672  0.9949008
# log(range.size)      -0.5746644 -0.3830046
# log(bio1_mean + 100) -0.9246153 -0.1826947

## minimum elevation of the highest grid cell (i.e. the highest that they for sure live at) - add 2 to make all response values positive
max.min00a <- gls(log(max.min.elev + 2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
max.min00b <- gls(log(max.min.elev + 2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
max.min00c <- gls(log(max.min.elev + 2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
max.min00d <- gls(log(max.min.elev + 2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))
AICc(max.min00a, max.min00b, max.min00c, max.min00d)
           # df      AICc
# max.min00a  6 1224.9230
# max.min00b  7  962.2820
# max.min00c  7  962.2820
# max.min00d  6  964.7921

summary(max.min00b) 
# Generalized least squares fit by REML
  # Model: log(max.min.elev + 2) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
      # AIC      BIC    logLik
  # 961.901 987.7571 -473.9505

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.2200327 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          -4.195007  2.396818 -1.750240  0.0811
# log(body.size)        0.282506  0.351925  0.802744  0.4228
# log(rel.wing)         7.497236  3.280521  2.285379  0.0230
# log(range.size)       0.593579  0.051603 11.502879  0.0000
# log(bio1_mean + 100) -1.090433  0.213991 -5.095686  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.443                       
# log(rel.wing)         0.273  0.042                
# log(range.size)      -0.674 -0.124 -0.125         
# log(bio1_mean + 100) -0.619 -0.119 -0.162   0.321 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -3.41636783 -0.46534777 -0.01526039  0.65513334  2.33561340 

# Residual standard error: 1.275748 
# Degrees of freedom: 302 total; 297 residual

confint(max.min00b)
                          # 2.5 %     97.5 %
# (Intercept)          -8.8926844  0.5026707
# log(body.size)       -0.4072550  0.9722667
# log(rel.wing)         1.0675324 13.9269394
# log(range.size)       0.4924398  0.6947185
# log(bio1_mean + 100) -1.5098486 -0.6710176

### maximum elevation of the lowest grid cell (i.e. lowest that we are certain they live at) - add 1 to make all response values positive
min.max00a <- gls(log(min.max.elev + 1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
min.max00b <- gls(log(min.max.elev + 1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
min.max00c <- gls(log(min.max.elev + 1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
min.max00d <- gls(log(min.max.elev + 1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))
AICc(min.max00a, min.max00b, min.max00c, min.max00d)
           # df      AICc
# min.max00a  6 1194.7680
# min.max00b  7  970.1573
# min.max00c  7  970.1573
# min.max00d  6  977.1353

summary(min.max00b)
# Generalized least squares fit by REML
  # Model: log(min.max.elev + 1) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 969.7764 995.6325 -477.8882

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5576269 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          22.136458  2.562601  8.638277  0.0000
# log(body.size)        0.489491  0.396995  1.232990  0.2186
# log(rel.wing)         2.264530  3.510284  0.645113  0.5194
# log(range.size)      -0.493546  0.052004 -9.490524  0.0000
# log(bio1_mean + 100) -1.171496  0.224213 -5.224931  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.469                       
# log(rel.wing)         0.250  0.001                
# log(range.size)      -0.611 -0.141 -0.072         
# log(bio1_mean + 100) -0.576 -0.137 -0.131   0.312 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.2631594 -0.5522861 -0.1984277  0.4899144  1.7805132 

# Residual standard error: 1.641878 
# Degrees of freedom: 302 total; 297 residual

confint(min.max00b)
                          # 2.5 %     97.5 %
# (Intercept)          17.1138522 27.1590637
# log(body.size)       -0.2886052  1.2675866
# log(rel.wing)        -4.6155008  9.1445609
# log(range.size)      -0.5954726 -0.3916202
# log(bio1_mean + 100) -1.6109450 -0.7320472

##### mean elevation of highest occupied grid cell +  residual wing size
r.mod05a <- gls(log(max.mean.elev) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
r.mod05b <- gls(log(max.mean.elev) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
r.mod05c <- gls(log(max.mean.elev) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
r.mod05d <- gls(log(max.mean.elev) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(r.mod05a, r.mod05b, r.mod05c, r.mod05d)
         # df     AICc
# r.mod05a  6 978.3234
# r.mod05b  7 737.8048
# r.mod05c  7 737.8048
# r.mod05d  6 754.2591

summary(r.mod05b)
# Generalized least squares fit by REML
  # Model: log(max.mean.elev) ~ log(body.size) + resid.wing + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC    BIC    logLik
  # 737.4239 763.28 -361.7119

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5175356 

# Coefficients:
                          # Value Std.Error   t-value p-value
# (Intercept)           1.7735006 1.7223327  1.029708  0.3040
# log(body.size)        0.0994906 0.2753715  0.361296  0.7181
# resid.wing            0.1903666 0.0752998  2.528116  0.0120
# log(range.size)       0.3491355 0.0348243 10.025643  0.0000
# log(bio1_mean + 100) -0.8604832 0.1494672 -5.757005  0.0000

 # Correlation: 
                     # (Intr) lg(b.) rsd.wn lg(r.)
# log(body.size)       -0.529                     
# resid.wing            0.293 -0.298              
# log(range.size)      -0.615 -0.110 -0.076       
# log(bio1_mean + 100) -0.581 -0.088 -0.134  0.313

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -2.91298725 -0.56359896 -0.01971429  0.55110349  1.80670477 

# Residual standard error: 1.058088 
# Degrees of freedom: 302 total; 297 residual

confint(r.mod05b)
                           # 2.5 %     97.5 %
# (Intercept)          -1.60220948  5.1492107
# log(body.size)       -0.44022752  0.6392088
# resid.wing            0.04278172  0.3379514
# log(range.size)       0.28088125  0.4173898
# log(bio1_mean + 100) -1.15343346 -0.5675330

## elevational breadth with residual wing size - add 0.1 to make all response values positive
r.r00a <- gls(log(range.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
r.r00b <- gls(log(range.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
r.r00c <- gls(log(range.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
r.r00d <- gls(log(range.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(r.r00a, r.r00b, r.r00c, r.r00d)
       # df     AICc
# r.r00a  6 1309.054
# r.r00b  7 1045.912
# r.r00c  7 1045.912
# r.r00d  6 1046.578

summary(r.r00d)
# Generalized least squares fit by REML
  # Model: log(range.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 1046.293 1068.455 -517.1464

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          -9.672261 2.4915706 -3.881994  0.0001
# log(body.size)       -0.187076 0.2851352 -0.656095  0.5123
# resid.wing            0.251706 0.0840826  2.993549  0.0030
# log(range.size)       0.751498 0.0576144 13.043577  0.0000
# log(bio1_mean + 100) -0.582295 0.2230271 -2.610870  0.0095

 # Correlation: 
                     # (Intr) lg(b.) rsd.wn lg(r.)
# log(body.size)       -0.449                     
# resid.wing            0.372 -0.003              
# log(range.size)      -0.776  0.016 -0.386       
# log(bio1_mean + 100) -0.670 -0.014 -0.278  0.318

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.59253157 -0.43807609  0.06115794  0.58757180  1.86477648 

# Residual standard error: 1.324674 
# Degrees of freedom: 302 total; 297 residual

confint(r.r00d)
                            # 2.5 %     97.5 %
# (Intercept)          -14.55564982 -4.7888727
# log(body.size)        -0.74593039  0.3717790
# resid.wing             0.08690659  0.4165044
# log(range.size)        0.63857590  0.8644203
# log(bio1_mean + 100)  -1.01941954 -0.1451696

## mean elevation of lowest occupied grid cell + resid wing size - add 0.1 to make all response values positive
r.min00a <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
r.min00b <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
r.min00c <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
r.min00d <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(r.min00a, r.min00b, r.min00c, r.min00d)
         # df      AICc
# r.min00a  6 1188.8120
# r.min00b  7  983.7706
# r.min00c  7  983.7706
# r.min00d  6  986.5249

summary(r.min00b)
# Generalized least squares fit by REML
  # Model: log(min.mean.elev + 0.1) ~ log(body.size) + resid.wing + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 983.3897 1009.246 -484.6948

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.5479602 

# Coefficients:
                         # Value Std.Error    t-value p-value
# (Intercept)          23.202270 2.6183996   8.861241  0.0000
# log(body.size)        0.295068 0.4198505   0.702792  0.4827
# resid.wing            0.114698 0.1142138   1.004244  0.3161
# log(range.size)      -0.571347 0.0526213 -10.857720  0.0000
# log(bio1_mean + 100) -1.121007 0.2266251  -4.946527  0.0000

 # Correlation: 
                     # (Intr) lg(b.) rsd.wn lg(r.)
# log(body.size)       -0.529                     
# resid.wing            0.291 -0.299              
# log(range.size)      -0.609 -0.112 -0.073       
# log(bio1_mean + 100) -0.577 -0.091 -0.132  0.312

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.7973431 -0.5566213 -0.1297642  0.4308258  2.0877056 

# Residual standard error: 1.645537 
# Degrees of freedom: 302 total; 297 residual

confint(r.min00b)
                          # 2.5 %     97.5 %
# (Intercept)          18.0703008 28.3342385
# log(body.size)       -0.5278241  1.1179595
# resid.wing           -0.1091564  0.3385533
# log(range.size)      -0.6744832 -0.4682115
# log(bio1_mean + 100) -1.5651844 -0.6768301

############# predictors of high elevation life - finer resolution
f.mod00a <- gls(log(max.mean.elev.finer) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
f.mod00b <- gls(log(max.mean.elev.finer) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
f.mod00c <- gls(log(max.mean.elev.finer) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
f.mod00d <- gls(log(max.mean.elev.finer) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(f.mod00a, f.mod00b, f.mod00c, f.mod00d)
         # df     AICc
# f.mod00a  6 955.8267
# f.mod00b  7 716.6999
# f.mod00c  7 716.6999
# f.mod00d  6 726.9044

summary(f.mod00b)
# Generalized least squares fit by REML
  # Model: log(max.mean.elev.finer) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio1_mean + 100) 
  # Data: drags.dat 
      # AIC      BIC    logLik
  # 716.319 742.1751 -351.1595

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.3611738 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)           1.810259 1.6226413  1.115625  0.2655
# log(body.size)        0.193266 0.2456787  0.786662  0.4321
# log(rel.wing)         6.079748 2.2404753  2.713598  0.0070
# log(range.size)       0.371807 0.0341318 10.893269  0.0000
# log(bio1_mean + 100) -0.924974 0.1439874 -6.423993  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.458                       
# log(rel.wing)         0.261  0.021                
# log(range.size)      -0.647 -0.132 -0.096         
# log(bio1_mean + 100) -0.602 -0.126 -0.146   0.318 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -2.90903567 -0.59209844 -0.03065673  0.54603681  2.24496276 

# Residual standard error: 0.9173663 
# Degrees of freedom: 302 total; 297 residual

confint(f.mod00b)
                          # 2.5 %     97.5 %
# (Intercept)          -1.3700592  4.9905778
# log(body.size)       -0.2882552  0.6747877
# log(rel.wing)         1.6884976 10.4709993
# log(range.size)       0.3049098  0.4387040
# log(bio1_mean + 100) -1.2071839 -0.6427638

###### elevational breadth using finer resolution - add 0.1 to make all response values positive
f.r00a <- gls(log(range.elev.finer + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
f.r00b <- gls(log(range.elev.finer + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
f.r00c <- gls(log(range.elev.finer + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
f.r00d <- gls(log(range.elev.finer + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(f.r00a, f.r00b, f.r00c, f.r00d)
       # df      AICc
# f.r00a  6 1171.8338
# f.r00b  7  938.8000
# f.r00c  7  938.8000
# f.r00d  6  939.5454

summary(f.r00d)
# Generalized least squares fit by REML
  # Model: log(range.elev.finer + 0.1) ~ log(body.size) + log(rel.wing) +      log(range.size) + log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC     BIC    logLik
  # 939.2606 961.423 -463.6303

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          -7.233311 2.0797472 -3.477976  0.0006
# log(body.size)        0.117718 0.2531677  0.464980  0.6423
# log(rel.wing)         6.716837 2.2065442  3.044053  0.0025
# log(range.size)       0.671636 0.0486744 13.798551  0.0000
# log(bio1_mean + 100) -0.706061 0.1884199 -3.747273  0.0002

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.326                       
# log(rel.wing)         0.343  0.308                
# log(range.size)      -0.772 -0.105 -0.386         
# log(bio1_mean + 100) -0.668 -0.100 -0.278   0.318 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -5.66465961 -0.50700143 -0.02542197  0.62291407  1.87447867 

# Residual standard error: 1.119124 
# Degrees of freedom: 302 total; 297 residual

confint(f.r00d)
                           # 2.5 %     97.5 %
# (Intercept)          -11.3095403 -3.1570810
# log(body.size)        -0.3784817  0.6139174
# log(rel.wing)          2.3920903 11.0415845
# log(range.size)        0.5762360  0.7670361
# log(bio1_mean + 100)  -1.0753573 -0.3367647

###### min elevation using finer resolution - need to add 56.2 to make all response values positive
f.min00a <- gls(log(min.mean.elev.finer + 56.2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
f.min00b <- gls(log(min.mean.elev.finer + 56.2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
f.min00c <- gls(log(min.mean.elev.finer + 56.2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
f.min00d <- gls(log(min.mean.elev.finer + 56.2) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio1_mean + 100), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(f.min00a, f.min00b, f.min00c, f.min00d)
         # df     AICc
# f.min00a  6 1397.415
# f.min00b  7 1144.865
# f.min00c  7 1144.865
# f.min00d  6 1142.877

summary(f.min00d)
# Generalized least squares fit by REML
  # Model: log(min.mean.elev.finer + 56.2) ~ log(body.size) + log(rel.wing) +      log(range.size) + log(bio1_mean + 100) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 1142.593 1164.755 -565.2963

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                          # Value Std.Error   t-value p-value
# (Intercept)           20.423602 2.9286924  6.973625  0.0000
# log(body.size)         0.091143 0.3565098  0.255654  0.7984
# log(rel.wing)        -10.350086 3.1072474 -3.330950  0.0010
# log(range.size)       -0.528067 0.0685431 -7.704156  0.0000
# log(bio1_mean + 100)  -0.693547 0.2653323 -2.613882  0.0094

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.326                       
# log(rel.wing)         0.343  0.308                
# log(range.size)      -0.772 -0.105 -0.386         
# log(bio1_mean + 100) -0.668 -0.100 -0.278   0.318 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -4.3419595 -0.2149017  0.1612992  0.5162734  1.5053541 

# Residual standard error: 1.575947 
# Degrees of freedom: 302 total; 297 residual

confint(f.min00d)
                           # 2.5 %     97.5 %
# (Intercept)           14.6834709 26.1637339
# log(body.size)        -0.6076034  0.7898895
# log(rel.wing)        -16.4401794 -4.2599935
# log(range.size)       -0.6624086 -0.3937247
# log(bio1_mean + 100)  -1.2135890 -0.1735056


##### analyses using min temp of coldest month instead of mean annual temp 

##### maximum elevation
cold.max00a <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
cold.max00b <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
cold.max00c <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
cold.max00d <- gls(log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(cold.max00a, cold.max00b, cold.max00c, cold.max00d)
            # df     AICc
# cold.max00a  6 986.1218
# cold.max00b  7 745.7100
# cold.max00c  7 745.7100
# cold.max00d  6 755.0387

summary(cold.max00b)
# Generalized least squares fit by REML
  # Model: log(max.mean.elev) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio6_mean + 331) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 745.3291 771.1852 -365.6645

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
  # lambda 
# 0.442947 

# Coefficients:
                         # Value Std.Error   t-value p-value
# (Intercept)          -1.115110 1.5516657 -0.718653  0.4729
# log(body.size)        0.179690 0.2627179  0.683964  0.4945
# log(rel.wing)         5.063851 2.3652933  2.140898  0.0331
# log(range.size)       0.382744 0.0346001 11.061941  0.0000
# log(bio6_mean + 331) -0.431321 0.1028737 -4.192720  0.0000

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.564                       
# log(rel.wing)         0.239  0.002                
# log(range.size)      -0.599 -0.116 -0.062         
# log(bio6_mean + 331) -0.447 -0.088 -0.100   0.199 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -3.10865930 -0.53867504  0.01600551  0.54536581  1.96530506 

# Residual standard error: 1.01938 
# Degrees of freedom: 302 total; 297 residual

confint(cold.max00b)
                          # 2.5 %     97.5 %
# (Intercept)          -4.1563187  1.9260992
# log(body.size)       -0.3352280  0.6946073
# log(rel.wing)         0.4279611  9.6997404
# log(range.size)       0.3149295  0.4505594
# log(bio6_mean + 331) -0.6329496 -0.2296920

##### elevation breadth

cold.range00a <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
cold.range00b <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.9, form=~binom, phy = na.drags, fixed = FALSE))
cold.range00c <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
cold.range00d <- gls(log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(cold.range00a, cold.range00b, cold.range00c, cold.range00d)
              # df     AICc
# cold.range00a  6 1309.958
# cold.range00b  7 1043.677
# cold.range00c  7 1043.677
# cold.range00d  6 1042.527

summary(cold.range00d)
# Generalized least squares fit by REML
  # Model: log(range.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio6_mean + 331) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 1042.242 1064.405 -515.1212

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                          # Value Std.Error   t-value p-value
# (Intercept)          -11.709933 2.1760446 -5.381293  0.0000
# log(body.size)         0.069101 0.3003713  0.230054  0.8182
# log(rel.wing)          7.161350 2.5828591  2.772644  0.0059
# log(range.size)        0.773790 0.0560898 13.795558  0.0000
# log(bio6_mean + 331)  -0.343164 0.1594447 -2.152244  0.0322

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.400                       
# log(rel.wing)         0.300  0.304                
# log(range.size)      -0.768 -0.093 -0.358         
# log(bio6_mean + 331) -0.536 -0.086 -0.222   0.212 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.55958113 -0.43791748  0.07102217  0.62723696  1.81462138 

# Residual standard error: 1.329462 
# Degrees of freedom: 302 total; 297 residual

confint(cold.range00d)
                           # 2.5 %      97.5 %
# (Intercept)          -15.9749024 -7.44496451
# log(body.size)        -0.5196155  0.65781852
# log(rel.wing)          2.0990388 12.22366035
# log(range.size)        0.6638561  0.88372409
# log(bio6_mean + 331)  -0.6556698 -0.03065805


### minimum elevation
cold.min00a <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = TRUE))
cold.min00b <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 1, form=~binom, phy = na.drags, fixed = FALSE))
cold.min00c <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0.5, form=~binom, phy = na.drags, fixed = FALSE))
cold.min00d <- gls(log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) + log(bio6_mean + 331), data = drags.dat, na.action = na.exclude, correlation = corPagel(value = 0, form=~binom, phy = na.drags, fixed = TRUE))

AICc(cold.min00a, cold.min00b, cold.min00c, cold.min00d)
            # df      AICc
# cold.min00a  6 1190.4123
# cold.min00b  7  987.7298
# cold.min00c  7  987.7298
# cold.min00d  6  991.5375

summary(cold.min00b)
# Generalized least squares fit by REML
  # Model: log(min.mean.elev + 0.1) ~ log(body.size) + log(rel.wing) + log(range.size) +      log(bio6_mean + 331) 
  # Data: drags.dat 
       # AIC      BIC    logLik
  # 987.3489 1013.205 -486.6744

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.6645408 

# Coefficients:
                         # Value Std.Error    t-value p-value
# (Intercept)          19.092102  2.442915   7.815297  0.0000
# log(body.size)        0.328792  0.415691   0.790953  0.4296
# log(rel.wing)         2.797615  3.616273   0.773618  0.4398
# log(range.size)      -0.524290  0.051456 -10.189066  0.0000
# log(bio6_mean + 331) -0.536341  0.154303  -3.475900  0.0006

 # Correlation: 
                     # (Intr) lg(b.) lg(rl.) lg(rn.)
# log(body.size)       -0.569                       
# log(rel.wing)         0.226 -0.018                
# log(range.size)      -0.549 -0.123 -0.039         
# log(bio6_mean + 331) -0.407 -0.098 -0.085   0.186 

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -2.46567892 -0.51214683 -0.05925324  0.39990932  1.83827644 

# Residual standard error: 1.896811 
# Degrees of freedom: 302 total; 297 residual

confint(cold.min00b)

                          # 2.5 %     97.5 %
# (Intercept)          14.3040771 23.8801261
# log(body.size)       -0.4859476  1.1435322
# log(rel.wing)        -4.2901494  9.8853800
# log(range.size)      -0.6251418 -0.4234375
# log(bio6_mean + 331) -0.8387690 -0.2339132
