# compile with:
# roxygen2::roxygenise()
# or: Ctrl + Shift + D, if you’re using RStudio.

# Currently there is no way to build vignettes using devtools if you just use
# the RStudio button Build & Reload. You have to run:
# devtools::install(build_vignettes = TRUE)

#require(akima)

#' Read raw disdrometer data
#'
#' Retrieves (usually minutal) particle size and velocity distributio PSVD data,
#' that is the matrix of particle counts arranged by size and velocity classes,
#' from a list of raw disdrometer data files.
#'
#' @param files        A list of files to be processed.
#' @param type         Character vector designing the type of disdrometer,
#'                     currently one of 'Thies' or 'Parsivel' Defaults to
#'                     'Thies'.
#'
#' @return An array containing, for each minute, a particle size velocity
#' distribution (PSVD) matrix, i.e. a matrix of velocity (rows) vs.
#' size (columns) particle counts.
#'
#' @section References
#'
#' @examples
#'
#' f <- system.file('extdata', 'rawDataThies', package='disdRo')
#' files <- list.files(f, '.txt', full.names=TRUE, recursive=TRUE)
#' dsd <- dsd_read(files, type='Thies')
#' dim(dsd)
#' dimnames(dsd)
#'
#' @export
dsd_read <- function (files, type='Thies') {

  if (type!='Thies' & type!='Parsivel')
    stop('type must be one of c(Thies, Parsivel)')
  
  # particle size and velocity bins
  dia <- switch(type,
                Thies=c(0.125, 0.25, 0.375, 0.5, 0.750, 1, 1.250, 1.5, 1.75,
                        2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5),
                Parsivel=c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1,
                           1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4,
                           4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26))
  vel <- switch(type,
                Thies=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4,
                        4.2, 5, 5.8, 6.6, 7.4, 8.2, 9, 10, 11),
                Parsivel=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                           1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
                           4.8, 5.6, 6.4, 7.2, 8.0, 9.6, 11.2, 12.8, 14.4, 16.0,
                           19.2, 22.4))
  # particle size and velocity means
  dia_m <- switch(type,
                  Thies=(dia[1:22]+dia[2:23])/2,
                  Parsivel=(dia[1:32]+dia[2:33])/2)
  vel_m <- switch(type,
                  Thies=(vel[1:20]+vel[2:21])/2,
                  Parsivel=(vel[1:32]+vel[2:33])/2)
  
  # read the files
  dat <- NULL
  bas <- switch(type, Thies=80, Parsivel=95)
  imx <- switch(type, Thies=22, Parsivel=32)
  jmx <- switch(type, Thies=20, Parsivel=32)
  for (f in files) {
    dsd <- read.table(f, sep=';',header=FALSE)[,bas:{bas+jmx*imx-1}]
    date <- strsplit(rev(strsplit(f,'/')[[1]])[1],'.txt')[[1]][1]
  	ini <- strptime(date,'%Y%m%d%H')+60
  	n <- nrow(dsd)
  	#dates <- seq(ini,ini+3600-60,by='min')
  	#dates <- seq(ini,ini+3600-(61-n)*n,by='min')
  	dates <- seq(ini,ini+(n-1)*60,by='min')
  	dsd <- array(as.matrix(dsd),dim=c(n,jmx,imx))
  	if (type=='Parsivel') {
      dsd <- aperm(dsd, c(1,3,2))
      w <- which(dsd==256)
      dsd[w] <- 0
  	}
  	dimnames(dsd) <- list(time=as.character(dates),
  	                      velocity=vel_m, diameter=dia_m)
  	dat <- abind::abind(dat, dsd, along=1)
  }

  # return
  return(dat)
}





#' Plot disdrometer data
#'
#' Produce a PSVD plot: particle count velocity vs. size.
#'
#' @param x            A particle size velocity distribution (PSVD) matrix.
#' @param type         Character vector designing the type of disdrometer,
#'                     currently one of 'Thies' or 'Parsivel'. Defaults to
#'                     'Thies'.
#' @param model        Vector. Which theoretical models of V vs. DS curves to
#'                     plot on top of the PSVD. Defaults to c('Atlas',
#'                     'Uplinger','VanDijk')
#' @param contour      Logical: should 2d density estimate contour lines be
#'                     added to the plot? Defaults to FALSE.
#' @param theme       Character vector indicating a plotting theme to use.
#'                    Current options are 'color' (default) or 'bw' (black and
#'                    white).
#'
#' @return A DSD plot.
#'
#' @section References
#'
#' @examples
#'
#' f <- system.file('extdata', 'rawDataThies', package='disdRo')
#' files <- list.files(f, '.txt', full.names=TRUE, recursive=TRUE)
#' dsd <- dsd_read(files, type='Thies')
#' day <- apply(dsd, c(2,3), sum)
#' head(day)
#' # full plot
#' dsd_plot(day)
#' # choose one model
#' dsd_plot(day, model='VanDijk')
#' # no model
#' dsd_plot(day, model=NA)
#' # no model, add contour lines
#' dsd_plot(day, model=NA, contour=TRUE)
#'
#' @export
dsd_plot <- function(x, type='Thies',
                     model=c('Atlas','Uplinger','VanDijk'),
                     contour=FALSE, theme='color') {

  if (type!='Thies' & type!='Parsivel')
    stop('type must be one of c(Thies, Parsivel)')

  # particle size and velocity bins
  dia <- switch(type,
                Thies=c(0.125, 0.25, 0.375, 0.5, 0.750, 1, 1.250, 1.5, 1.75,
                        2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5),
                Parsivel=c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1,
                           1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 3.5, 4,
                           4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 26))
  vel <- switch(type,
                Thies=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.4, 1.8, 2.2, 2.6, 3, 3.4,
                        4.2, 5, 5.8, 6.6, 7.4, 8.2, 9, 10, 11),
                Parsivel=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                           1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
                           4.8, 5.6, 6.4, 7.2, 8.0, 9.6, 11.2, 12.8, 14.4, 16.0,
                           19.2, 22.4))
  # particle size and velocity bin widths
  dia_w <- switch(type,
                  Thies=-(dia[1:22]-dia[2:23]),
                  Parsivel=-(dia[1:32]-dia[2:33]))
  vel_w <- switch(type,
                  Thies=-(vel[1:20]-vel[2:21]),
                  Parsivel=-(vel[1:32]-vel[2:33]))
  # particle size and velocity means
  dia_m <- switch(type,
                  Thies=(dia[1:22]-dia[2:23]),
                  Parsivel=(dia[1:32]-dia[2:33]))
  vel_m <- switch(type,
                  Thies=(vel[1:20]-vel[2:21]),
                  Parsivel=(vel[1:32]-vel[2:33]))
  
  # transform to long format for ggplot
  x_m <- reshape2::melt(x)
  xx <- x
  dimnames(xx) <- list(vel_w, dia_w)
  x_m <- cbind(x_m, reshape2::melt(xx))[,-3]
  colnames(x_m) <- c('vel','dia','vel_w','dia_w','n')
  x_m$n <- log10(x_m$n)
  x_m <- x_m[x_m$n>1,]
  #summary(x_m)

  # heatmap - rectangular binding
  # original version, casts a warning in R CMD CHECK: g <- ggplot(x_m, aes(x=dia, y=vel, width=dia_w, height=vel_w, fill=n)) +
  g <- ggplot(x_m, aes_(x=~dia, y=~vel, width=~dia_w, height=~vel_w, fill=~n)) +
    geom_tile(alpha=1) +
    xlim(c(0,8)) + ylim(c(0,15)) +
    xlab('Diameter (mm)') + ylab('Velocity (m/s)') +
    theme_bw() +
    theme(panel.grid=element_blank())
  if (theme=='color') {
    g <- g + scale_fill_gradient2(name='NS', low='darkgreen',mid='yellow',high='darkred',
                         limits=c(0,5), midpoint=2,
                         na.value='darkgreen',labels=c(0,10,100,1000,10000,100000))
    if (length(model)>0) {
      if ('Atlas' %in% model) {
        g <- g + stat_function(aes(col='Atlas'), fun=function(x) 17.67*(x/10)^0.67)
      }
      if ('Uplinger' %in% model) {
        g <- g + stat_function(aes(col='Uplinger'), fun=function(x) 4.874*x*exp(-0.195*x))
      }
      if ('VanDijk' %in% model) {
        g <- g + stat_function(aes(col='VanDijk'), fun=function(x) 0.0561*x^3-0.912*x^2+5.03*x-0.254)
      }
      g <- g + guides(col=guide_legend(title='Model'))
      if (contour) {
        g <- g + geom_density_2d(alpha=0.5)
        #g <- g + stat_density_2d(aes(fill = ..level..), geom = "polygon")
      }
    }
  }
  if (theme=='bw') {
    g <- g + scale_fill_gradient2(name='NS', low='white',mid='grey50',high='black',
                                  limits=c(0,5), midpoint=2,
                                  na.value='white',labels=c(0,10,100,1000,10000,100000))
    if (length(model)>0) {
      if ('Atlas' %in% model) {
        g <- g + stat_function(aes(linetype='Atlas'), fun=function(x) 17.67*(x/10)^0.67)
      }
      if ('Uplinger' %in% model) {
        g <- g + stat_function(aes(linetype='Uplinger'), fun=function(x) 4.874*x*exp(-0.195*x))
      }
      if ('VanDijk' %in% model) {
        g <- g + stat_function(aes(linetype='VanDijk'), fun=function(x) 0.0561*x^3-0.912*x^2+5.03*x-0.254)
      }
      g <- g + guides(col=guide_legend(title='Model'))
      if (contour) {
        g <- g + geom_density_2d(color='black', alpha=0.5)
        #g <- g + stat_density_2d(aes(fill = ..level..), geom = "polygon")
      }
    }
  }
  # 
  # if (int) {
  #   # kernel density estimation (interpolated values)?
  #   x_int <- akima::interp(y=x_m$dia, x=x_m$vel, z=x_m$n,
  #                   xo=seq(min(vel),max(vel),length=100),
  #                   yo=seq(min(dia),max(dia),length=100))
  #   x_int <- matrix(x_int$z, 100, 100, dimnames=list(x_int$x, x_int$y))
  #   x_int <- reshape2::melt(x_int, na.rm=TRUE)
  #   x_int <- x_int[!is.na(x_int[,3]),]
  #   x_int <- as.data.frame(x_int)
  #   colnames(x_int) <- colnames(x_m[c(1,2,5)])
  #   #summary(x_int)
  #   
  #   # plot it
  #   # original version, casts a warning in R CMD CHECK:     g <- ggplot(x_int, aes(x=dia,y=vel,z=n)) + geom_tile(aes(fill=n), alpha=1) +
  #   g <- ggplot(x_int, aes_(x=~dia, y=~vel, z=~n)) +
  #     geom_tile(aes_(fill=~n), alpha=1) +
  #     stat_contour(col='black', alpha=0.2) +
  #     scale_fill_gradient2(name='NS', low='darkgreen',mid='yellow',high='darkred',
  #                          limits=c(0,5), midpoint=2,
  #                          na.value='darkgreen',labels=c(0,10,100,1000,10000,100000)) +
  #     stat_function(fun=function(x) 17.67*(x/10)^0.67, col='blue') + # Atlas
  #     stat_function(fun=function(x) 4.874*x*exp(-0.195*x), col='red') + # Uplinger
  #     stat_function(fun=function(x) 0.0561*x^3-0.912*x^2+5.03*x-0.254) + # Van Dijk
  #     xlim(c(0,8)) + ylim(c(0,15)) +
  #     xlab('Diameter (mm)') + ylab('Velocity (m/s)') +
  #     theme_bw() +
  #     theme(panel.grid=element_blank())
  # }
  
  return(g)
}



  
#' Compute PSVD integrated variables
#'
#' \code{dsd_integrate} reads raw disdrometer data and computes a series of
#' integrated variables.
#'
#' Currently, this is done via an external Perl script, so you need to have
#' Perl installed and working in your system. Beware: some users have reported
#' issues for running the Perl script in Windows.
#' It might be translated into a native R script in the future, if I find time
#' to do it.
#'
#' @param indir   A character vector with the url of the directory tree that
#'                contains the raw disdrometric data.
#' @param script  A character vector with the url of the Perl script. Defaults
#'                to the installation directory of the \code{disdRo} package.
#' @param outfile A character vector with the url of the output file. Defaults
#'                to a random file name in the session's temporary directory.
#' @param interp  A character vector indicating the interpolation method to
#'                use for inputing particle sizes and velocities within the bins
#'                of the PSVD matrix. One of `middle`, `uniform`, and `linear`,
#'                defaulting to `middle`.
#'
#' @return A data frame with the following items:
#' \describe{
#'    \item{type}{Disdrometer type, currently one of 'Thies' or 'Parsivel' (Factor)}
#'    \item{serial}{Sensor serial number (Numeric)}
#'    \item{time}{Date and time of the record (POSIXct)}
#'    \item{seconds}{Number of seconds since 1970-01-01 00:00:00 (Numerica)}
#'    \item{synop}{Synop code (Factor)}
#'    \item{r}{Precipitation intensity computed from the PSVD matrix, mm h−1 (Numeric)}
#'    \item{p}{Precipitation amount computed from the PSVD matrix, mm (Numeric)}
#'    \item{m}{Liquid water content computed from the PSVD matrix, g m-3 (Numeric)}
#'    \item{z}{Radar reflectivity computed from the PSVD matrix, dB mm6 m−3 (Numeric)}
#'    \item{e}{Kinetic energy computed from the PSVD matrix, J m−2 mm−1 (Numeric)}
#'    \item{mor}{Visibility computed from the PSVD matrix, m (Numeric)}
#'    \item{r_meas}{Precipitation intensity as reported by the disdrometer, mm h−1 (Numeric)}
#'    \item{z_meas}{Radar reflectivity as reported by the disdrometer, dB mm6 m−3 (Numeric)}
#'    \item{e_meas}{Kinetic energy as reported by the disdrometer, J m−2 h−1 (Numeric)}
#'    \item{mor_meas}{Visibility as reported by the disdrometer, m (Numeric)}
#'    \item{qual}{Data quality reported by the distrometer, 0-100 (Numeric)}
#'    \item{tmp}{Air temperature, ºC (Numeric)}
#'    \item{rh}{Relative humidity, 0-100 (Numeric)}
#'    \item{w}{Wind speed, m s-1 (Numeric)}
#'    \item{wd}{Wind direction, degrees (Numeric)}
#'    \item{np}{Number of particles detected computed from the PSVD matrix (Numeric)}
#'    \item{np_meas}{Number of particles detected as reported by the disdrometer (Numeric)}
#'    \item{lcurrent}{Laser control output, 1/100 mA (Numeric)}
#'    \item{ocontrol}{Optical control output, mV (Numeric)}
#'    \item{power}{Sensor power supply, V (Numeric)}
#'    \item{tmp_int}{Internal sensor temperature, ºC (Numeric)}
#'    \item{d10}{Particle diameter’s 10th percentile, mm (Numeric)}
#'    \item{d25}{Particle diameter’s 25th percentile, mm (Numeric)}
#'    \item{d50}{Particle diameter’s 50th percentile, mm (Numeric)}
#'    \item{d75}{Particle diameter’s 75th percentile, mm (Numeric)}
#'    \item{d90}{Particle diameter’s 90th percentile, mm (Numeric)}
#'    \item{dmean}{Mean particle diameter, mm (Numeric)}
#'    \item{v10}{Particle velocity’s 10th percentile, m s−1 (Numeric)}
#'    \item{v25}{Particle velocity’s 25th percentile, m s−1 (Numeric)}
#'    \item{v50}{Particle velocity’s 50th percentile, m s−1 (Numeric)}
#'    \item{v75}{Particle velocity’s 75th percentile, m s−1 (Numeric)}
#'    \item{v90}{Particle velocity’s 90th percentile, m s−1 (Numeric)}
#'    \item{vmean}{Mean particle velocity, m s−1 (Numeric)}
#'    \item{t_shift}{Telegram shift time, s (Numeric)}
#'    \item{nrow}{Telegram number of rows (Numeric)}
#'    \item{err}{Error status (Numeric)}
#'    \item{ncol}{Number of fields in the telegram (Numeric)}
#' }
#'
#' @section Bin interpolation:
#' Since the particle size and velocity distribution is not linear within the
#' bins of the PSVD matrix, different imputation methods exist. 'middle' will
#' assing the middle bin size and velocity to all the particles in the bin; 
#' 'uniform' assumes an uniform distribution of sizes and velocities within the
#' bin limits; and 'linear' assumes a linear distribution between the bin
#' limits.

#' @section Error codes:
#' \describe{
#'    \item{0}{No error}
#'    \item{1}{There is no telegram for that minute}
#'    \item{2}{Saturation: the limit of 9999 particles has been exceeded in at
#'    least one bin of the PSVD matrix (only for Thies)}
#'    \item{3}{Non conform characters found in SYNOP field}
#'    \item{4}{Non conform characters found in rain intensity field}
#'    \item{5}{9999.999 value found in rain intensity}
#'    \item{6}{The telegram only contains 'OK' or 'Version' (Parsivel)}
#'    \item{7}{Non conform characters found in the telegram}
#'    \item{21-23}{Parsivel error codes 1 to 3}
#'    \item{24-36}{Thies error codes 24 to 36}
#'    \item{37}{Multiple error codes in the telegram (Thies)}
#'}
#'
#' @examples
#'
#' f <- system.file('extdata/rawDataParsivel', package='disdRo')
#' x <- dsd_integrate(f)
#' head(x)
#' summary(x)
#'
#' @export
dsd_integrate <- function(indir, script=NULL, outfile=NULL, interp='linear') {
  if (is.null(outfile)) {
    outfile <- paste(tempdir(), do.call(paste0,
                     replicate(15, sample(LETTERS, 1, TRUE), FALSE)), sep='/')
  }
  if (is.null(script)) {
    #script <- './perl/process.pl'
    #script <- paste(find.package('disdRo'),'perl/process.pl',sep='/')
    script <- system.file('perl/process.pl', package='disdRo')
  }
  system(paste('perl', script, indir, outfile, interp))
  int <- read.table(outfile, sep=',', header=TRUE, na.strings='na')
  int$time <- as.POSIXct(int$time)
  return(int)
}
