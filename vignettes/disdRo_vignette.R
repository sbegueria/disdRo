## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo=TRUE)
knitr::opts_chunk$set(comment=NA, fig.height=6, fig.width=7, fig.align='center')

# Remember that RStudio’s "Build & reload" does not build vignettes to
# save time. Use `devtools::build(vignettes=TRUE)` instead, or
# altrnatively: `tools::buildVignettes(dir = '.', tangle = TRUE)`.
# Similarly, `devtools::install_github()` will also not build vignettes by
# default, so it need to be forced with
# `devtools::install_github(build_vignettes = TRUE)`.

## ----source--------------------------------------------------------------
library(disdRo)

## ----thiesRawDataFiles, message=FALSE, warning=FALSE---------------------
f <- system.file('extdata', 'rawDataThies', package='disdRo')
files <- list.files(f, '.txt', full.names=TRUE, recursive=TRUE)
PSVD_T <- psvd_read(files, type='Thies')
dim(PSVD_T)

## ----globalCounts, message=FALSE, warning=FALSE--------------------------
day_T <- apply(PSVD_T, c(2,3), sum)
head(day_T)

barplot(colSums(day_T), main='Daily drop counts per size class')
barplot(rowSums(day_T), main='Daily drop counts per velocity class')

## ----logDSDscatter, message=FALSE, warning=FALSE-------------------------
plot(log10(colSums(day_T))~colnames(day_T), ylab='log10(N)', xlab='size (mm)', type='b')
grid()

## ----timeCounts, message=FALSE, warning=FALSE----------------------------
time <- apply(PSVD_T, 1, sum)
barplot(time, main='Number of drops per minute')

## ----dsdPlot, message=FALSE, warning=FALSE-------------------------------
psvd_plot(day_T, type='Thies')

## ----dsdPlot2, message=FALSE, warning=FALSE------------------------------
psvd_plot(day_T, model='Beard')

## ----dsdPlot3, message=FALSE, warning=FALSE------------------------------
psvd_plot(day_T, model='Beard', alt=2000)

## ----dsdPlot4, message=FALSE, warning=FALSE------------------------------
psvd_plot(day_T, model=NA, contour=TRUE)

## ----dsdPlot5, message=FALSE, warning=FALSE------------------------------
psvd_plot(day_T, theme='bw')

## ----parsivelExample, message=FALSE, warning=FALSE-----------------------
f <- system.file('extdata', 'rawDataParsivel', package='disdRo')
files <- list.files(f, '.txt.gz', full.names=TRUE, recursive=TRUE)
PSVD_P <- psvd_read(files, type='Parsivel')

day_P <- apply(PSVD_P, c(2,3), sum)
barplot(colSums(day_P), main='Daily drop counts per size class')
barplot(rowSums(day_P), main='Daily drop counts per velocity class')

plot(log10(colSums(day_P))~colnames(day_P), ylab='log10(N)', xlab='size (mm)', type='b')
grid()

time <- apply(PSVD_P, 1, sum)
barplot(time, main='Number of drops per minute')

psvd_plot(day_P, type='Parsivel', model='Beard')

## ----filter1, message=FALSE, warning=FALSE-------------------------------
flt <- psvd_filter(type='Thies', d=c(0.3,Inf))
image(flt)

## ----filter2, message=FALSE, warning=FALSE-------------------------------
flt <- psvd_filter(type='Thies', tau=0.5, alt=2000)
image(flt)

## ----psvdPlotFilter, message=FALSE, warning=FALSE------------------------
psvd_plot(day_T, model='Beard', filter=flt, alpha=0)

## ----psdPlot1, message=FALSE, warning=FALSE------------------------------
psd_plot(day_T)

## ----psdPlot2, message=FALSE, warning=FALSE------------------------------
psd_plot(day_T, filter=flt)

## ----pvdPlot1, message=FALSE, warning=FALSE------------------------------
pvd_plot(day_T)

## ----pvdPlot2, message=FALSE, warning=FALSE------------------------------
pvd_plot(day_T, filter=flt)

## ----perlscript, message=FALSE, warning=FALSE----------------------------
?dsd_integrate

f <- system.file('extdata', 'rawDataThies', package='disdRo')
int <- psvd_integrate(f)
summary(int)

## ----perlscript2, fig.height=12, fig.width=7, message=FALSE, warning=FALSE----
par(mfrow=c(3,1))
plot(int$r~int$time, type='l', xlab='', ylab='I (mm/h)',
     main='Precipitation rate')
plot(int$e~int$time, type='l', xlab='', ylab='ET (J m-2 mm-1)',
     main='Unit kinetic energy')
plot(int$z~int$time, type='l', xlab='', ylab='Z (dB mm6 m-3)',
     main='Radar reflectivity')

## ----perlscript3, fig.height=8, fig.width=7, message=FALSE, warning=FALSE----
par(mfrow=c(2,1))
plot(int$z~int$time, type='l', xlab='', ylab='I (mm/h)',
     main='Precipitation rate (computed)')
plot(int$z_meas~int$time, type='l', xlab='', ylab='I (mm/h)',
     main='Precipitation rate (measured)')

## ----perlscript4, message=FALSE, warning=FALSE---------------------------
par(mfrow=c(1,1))
plot(int$z~int$z_meas, ylab='Calculated', xlab='Measured')
grid()
abline(0,1)

## ----window, fig.height=8, fig.width=7, message=FALSE, warning=FALSE-----
library(zoo)
int <- zoo(int[,-c(1:3)], int$time)

event <- window(int, start='2013-06-07 28:00:00', end='2013-06-07 23:59:00')
plot(event[,c('r_meas','z_meas','d50')], type='l',
     xlab='', main='A precipitation event')

## ----echo=FALSE, out.width="400px"---------------------------------------
knitr::include_graphics('img/percentile_distribution.001.jpg')

## ----echo=FALSE, out.width="400px"---------------------------------------
knitr::include_graphics('img/percentile_distribution.002.jpg')

## ----echo=FALSE, out.width="400px"---------------------------------------
knitr::include_graphics('img/percentile_distribution.003.jpg')

