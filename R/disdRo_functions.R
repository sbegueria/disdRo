# compile with:
# roxygen2::roxygenise()
# or: Ctrl + Shift + D, if you’re using RStudio.

# Currently there is no way to build vignettes using devtools if you just use
# the RStudio button Build & Reload. You have to run:
# devtools::install(build_vignettes = TRUE)

#require(akima)

#' Read raw disdrometer data
#'
#' Retrieves (usually minutal) particle size and velocity distribution (PSVD)
#' data, that is the matrix of particle counts arranged by size and velocity
#' binds, from a list of raw disdrometer data files.
#'
#' @param files        A list of files to be processed.
#' @param type         Character vector designing the type of disdrometer,
#'                     currently one of 'Thies' or 'Parsivel' Defaults to
#'                     'Thies'.
#'
#' @return An array containing, for each minute, a particle size velocity
#' distribution (PSVD) matrix, i.e. a matrix of size (rows) vs.
#' velocity (columns) particle counts.
#'
#' @section Note
#' So far, it is assumed that the data consists on the complete telegram is
#' recorded. The raw PSVD matrix is assumed to start on position 23 (Thies)
#' and 33 (Parsivel) of the telegram. This may not correspond to the factory
#' settings of these devices. Custom definition of the telegram needs to be
#' implemented.
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
  
  # particle size and velocity bin limits
#  dia <- switch(type, Thies=disdRo:::dia_t, Parsivel=disdRo:::dia_p)
#  vel <- switch(type, Thies=disdRo:::vel_t, Parsivel=disdRo:::vel_p)

  # particle size and velocity means
  dia_m <- switch(type, Thies=disdRo:::dia_m_t, Parsivel=disdRo:::dia_m_p)
  vel_m <- switch(type, Thies=disdRo:::vel_m_t, Parsivel=disdRo:::vel_m_p)
  
  # read the files
  dat <- NULL
  bas <- switch(type, Thies=80, Parsivel=95)
  imx <- switch(type, Thies=22, Parsivel=32) # raw DSD matrix starts at position 23 / 33
  jmx <- switch(type, Thies=20, Parsivel=32)
  for (f in files) {
    dsd <- read.table(f, sep=';',header=FALSE)[,bas:{bas+jmx*imx-1}]
    date <- strsplit(rev(strsplit(f,'/')[[1]])[1],'.txt')[[1]][1]
  	ini <- strptime(date,'%Y%m%d%H')+60
  	n <- nrow(dsd)
  	#dates <- seq(ini,ini+3600-60,by='min')
  	#dates <- seq(ini,ini+3600-(61-n)*n,by='min')
  	dates <- seq(ini,ini+(n-1)*60,by='min')
  	dsd <- array(as.matrix(dsd), dim=c(n,jmx,imx))
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
  dat <- aperm(dat, c(1,3,2))
  return(dat)
}





  
  
  
#' Terminal drop velocity theoretical models
#'
#' This function implements the Beard model of terminal drop velocity. This
#' is a physical approximation accounting for different factors such as viscous
#' drag, drop oblateness, etcetera. The function implemented here is valid for
#' drop sizes between 0.019 mm (19 µm) up to 7 mm. The function also contains
#' the approximations of Atlas, Uplinger and Van Dijk. The Beard model is
#' implemented by default assuming a standard atmosphere: P=101.325 kPa,
#' T=288.15 K, and rho= 1.225 kg/m3, but other values can be provided.
#'
#' @param size          Numeric. A vector of drop sizes (in mm).
#' @param model         Character. Name of the terminal velocity model to use.
#'                      One 'Beard', 'Atlas', 'Uplinger' or 'VanDijk'. Defaults
#'                      to 'Beard'.
#' @param eta           Numeric. A vector of air dynamic viscosity values.
#'                      Deafults to 1.818e-04 kg m-1 s-1. Optional, used only
#'                      in the Beard model.
#' @param P             Numeric. A vector of air pressure values. Defaults to
#'                      101.325 kPa. Optional, used  only in the Beard model.
#' @param T             Numeric. A vector of air temperature values. Defaults
#'                      to 288.15 K. Optional, used  only in the Beard model.
#' @param rho           Numeric. A vector of air densities. Defaults to
#'                      1.225 kg m-3. Optional, used  only in the Beard and
#'                      Atlas models.
#' @param alt           Numeric. Elevation, which can be used if provided to
#'                      estimate pressure, density and temperature of air in
#'                      according to the International Standard Atmosphere.
#'                      Defaults to 0 m above mean sea level. Optional, used 
#'                      in the Beard and Atlas models to estimate parameters
#'                      `eta`, `P`, `T` and `rho` if no values are provided.
#'
#' @return A vector of terminal fall velocities. It will cointain NA for any
#' drop sizes lower than 0.019 mm or larger than 7 mm (Beard model), or lower
#' than 0.1 and larger than 7 mm for the approximations.
#'
#' @section References
#' 
#' Atlas, D., Srivastava, R.C., Sekhon, R.S., 1973. Doppler radar
#' characteristics of precipitation at vertical incidence. Rev. Geophys.
#' Space Phys. 11, 1–35.
#' 
#' Beard, K.V., 1976. Terminal velocity and shape of cloud and precipitation
#' drops aloft. J. Atmos. Sci. 33, 851–864.
#' 
#' Uplinger, C.W. (1981) A new formula for raindrop terminal velocity. 20th
#' Conference of radar meteorology. American Meteorology Society, Boston (USA)
#' 389-391.
#' 
#' Van Dijk, A.I.J.M., Bruijnzeel, L.A., Rosewell, C.J., 2002. Rainfall
#' intensity-kinetic energy relationships: a critical literature appraisal.
#' J. Hydrol. 261, 1–23.
#' 
#' @examples
#' 
#' size <- seq(0, 6, 0.2)
#' psvd_model(size)
#' psvd_model(size, alt=2000)
#' psvd_model(size, model='Atlas')
#' psvd_model(size, model='Uplinger')
#' psvd_model(size, model='VanDijk')
#' 
#' @export
#' 
psvd_model <- function(size, model='Beard', eta=1.818e-05, P=101.325,
                       T=288.15, rho=1.225, alt=0) {
  if (length(eta)>1 & length(size) != length(eta))
    stop('Different parameter vector lengths.')
  if (length(P)>1 & length(size) != length(P))
    stop('Different parameter vector lengths.')
  if (length(T)>1 & length(size) != length(T))
    stop('Different parameter vector lengths.')
  if (!(model %in% c('Beard', 'Atlas', 'Uplinger','VanDijk')))
    stop('Model must be one of \'Beard\', \'Atlas\', \'Uplinger\', \'VanDijk\'.')

  # Calculate air parameters according to Standard Atmosphere
  # Following: https://en.wikipedia.org/wiki/Density_of_air for P, T, rho;
  # and the Sutherland's formula for eta (Crane Company. 1988. Flow of fluids
  # through valves, fittings, and pipe. Technical Paper No. 410 (TP 410)).
  if (alt != 0) {
    if (P == 101.325) {
      P <- 101.325 * (1 - (0.0065*alt)/288.15) ^(9.80665*0.0289644/8.31447/0.0065)
    }
    if (T == 288.15) {
      T <- 288.15 - (0.0065*alt)
    }
    if (rho == 1.225) {
      rho <- 1000 * (P*0.0289644) / (8.31447*T)
    }
    if (eta == 1.818e-05) {
      a <- 0.555*524.07 + 120
      b <- 0.555*1.8*T + 120
      eta <- 0.00001827 * (a/b) * (T/291.15)^(3/2)
    }
  }
  
  # Uplinger
  if (model=='Uplinger') {
    VInf <- 4.874 * size * exp(-0.195*size)
    w <- size <= 0.1 | size > 7
    VInf[w] <- NA
  }

  # Atlas
  if (model=='Atlas') {
    VInf <- (965 - 1030 * exp(-6*(size/10))) / 100
    # Correction for air density
    VInf <- VInf * (rho / 1.225)^0.4
    # equation of Atlas & Ulbrich (1977), used in Angulo (2016)
    #VInf <- 17.67*(size/10)^0.67
    w <- size <= 0.1 | size > 7
    VInf[w] <- NA
  }

  # Van Dijk
  if (model=='VanDijk') {
    VInf <- - 0.254 + 5.03*size - 0.912*size^2 + 0.0561*size^3
    w <- size <= 0.1 | size > 7
    VInf[w] <- NA
  }

  if (model=='Beard') {
    d0 <- size / 1000 # drop size, m
    deltaRho <- 1000 - rho # density difference (water - air)
    g <- 9.80665 # gravity acceleration, m s-2
    # for size <= 1.07 mm
    NDa <- 4 * rho * deltaRho * g * d0^3 / (3 * eta^2) # Davies number, adimens.
    X <- log(NDa)
    Y <- -0.318657e+01 + 0.992696*X - 0.153193e-02*X^2 - 0.987059e-03*X^3 -
      0.578878e-03*X^4 + 0.855176e-04*X^5 - 0.327815e-05*X^6
    l <- 6.62e-06 * (eta/1.818e-05) * (101.325/P) * (T/288.15)^(1/2) # mean free path, cm
    Csc <- 1 + 2.51 * l / d0 # slip correction factor - neglible for size>0.03 mm
    NRe <- Csc * exp(Y) # Reynolds number
    Vinf1 <- eta * NRe / (rho * d0) # terminal drop velocity, m s-1
    # for size > 1.07 mm
    sigma <- 0.073 # surface tension of water, N m-1 = kg s-2
    Bo <- 4 * deltaRho * g * d0^2 / (3 * sigma) # modified Bond number
    Np <- (sigma^3 * rho^2 / (eta^4 * deltaRho * g)) ^ (1/6) # physical property number (?)
    X <- log(Bo * Np)
    Y <- -0.500015e01 + 0.523778e01*X - 0.204914e01*X^2 + 0.475294*X^3 -
      0.542819e-01*X^4 + 0.238449e-02*X^5
    NRe <- Np * exp(Y) # Reynolds number
    Vinf2 <- eta * NRe / (rho * d0) # terminal drop velocity, m s-1
    #
    w <- size <= 1.07
    VInf <- c(Vinf1[which(w)],Vinf2[which(!w)])
    w <- size < 0.019
    VInf[w] <- NA
  }
  return(VInf)
}



#' Create a filter to remove outlier PSVD bins
#'
#' Produces a matrix that is used as a mask in further calculations that
#' involve using the PSVD matrix data. The purpose is to remove highly unlikely
#' drop size and velocity combinations, which are typical in disdrometer
#' records due to margin fallers, double drops, drop splashing, and other
#' issues.
#'
#' @param type          Character vector designing the type of disdrometer,
#'                      currently one of 'Thies' or 'Parsivel' Defaults to
#'                      'Thies'.
#' @param d             Numeric. A two-valued vector with the diameter limits
#'                      (inferior, superior). Defaults to (-Inf, Inf), so no
#'                      bins are removed.
#' @param v             Numeric. A two-valued vector with the velocity limits
#'                      (inferior, superior). Defaults to (-Inf, Inf), so no
#'                      bins are removed.
#' @param model         Character. Name of the terminal velocity model to use.
#'                      See `psvd_model()`.
#' @param tau           Numeric. A value between 0 and 1 that defines outlier
#'                      bins. Defaults to 0.5 (i.e., velocities that are 50%
#'                      off the theoretical model are removed). Defaults to
#'                      Inf, in which case all the bins are accepted and no
#'                      filtering is done according to a theoretical model.
#' @param eta           Numeric. Air dynamic viscosity, optional.
#'                      See `psvd_model()`.
#' @param P             Numeric. Air pressure, optional. See `psvd_model()`.
#' @param T             Numeric. Air temperature, optional. See `psvd_model()`.
#' @param rho           Numeric. Air density, optional. See `psvd_model()`.
#' @param alt           Numeric. Elevation, optional. See `psvd_model()`.
#'
#' @return A 22x20 (Thies) or 32x32 (Parsivel) logical matrix indicating
#' which PSVD bins to consider (TRUE) and which to remove (FALSE). Diameters
#' are stored as rows, and velocities as columns.
#'
#' @section References
#' 
#' 
#' @examples
#' 
#' # shown as images for easy visualization:
#' 
#' # filter only by drop size
#' image(psvd_filter(d=c(0.3,7), tau=Inf))
#' 
#' # filter only by theoretical velocity: low tolerance
#' image(psvd_filter(tau=0.5))
#' # higher tolerance
#' image(psvd_filter(tau=0.7))
#' # using other theoretical model
#' image(psvd_filter(model='Atlas', tau=0.5))
#' 
#' # filter by two criteria
#' image(psvd_filter(d=c(0.3,7), tau=0.5))
#' 
#' # filter for a Parsivel
#' image(psvd_filter(type='Parsivel', d=c(0.3,7), tau=0.5))
#' 
#' # used as input to `psvd_plot()`
#'  
#' 
#' @export
#' 
psvd_filter <- function(type='Thies', d=c(-Inf,Inf), v=c(-Inf,Inf),
                        model='Beard', tau=Inf, eta=1.818e-05,
                        P=101.325, T=288.15, rho=1.225, alt=0) {
  
  if (d[1]>d[2]) stop('Lower diameter limit is larger than higher.')
  if (v[1]>v[2]) stop('Lower velocity limit is larger than higher.')
  if (!is.infinite(tau) & (tau<0 | tau>1)) stop('tau must be between 0 and 1.')
  
  # particle size and velocity bin limits
  #  dia <- switch(type, Thies=disdRo:::dia_t, Parsivel=disdRo:::dia_p)
  #  vel <- switch(type, Thies=disdRo:::vel_t, Parsivel=disdRo:::vel_p)
  
  # particle size and velocity means
  dia_m <- switch(type, Thies=disdRo:::dia_m_t, Parsivel=disdRo:::dia_m_p)
  vel_m <- switch(type, Thies=disdRo:::vel_m_t, Parsivel=disdRo:::vel_m_p)
  
  
  fil <- matrix(TRUE, nrow=length(dia_m), ncol=length(vel_m))
  dimnames(fil) <- list(size=dia_m, velocity=vel_m)
  
  # Filter by diameter range
  w <- which(!(dia_m >= d[1] & dia_m <= d[2]))
  fil[w,] <- FALSE

  # Filter by velocity range
  w <- which(!(vel_m >= v[1] & vel_m <= v[2]))
  fil[,w] <- FALSE
  
  # Filter according to a theoretical terminal velocity model
  VInf <- psvd_model(dia_m, model)
  for (i in 1:length(dia_m)) {
    vratio <- vel_m/VInf[i]
    w <- which(vratio < (1-tau) | vratio > (1+tau))
    fil[i,w] <- FALSE
  }
  
  return(fil)
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
#'                     plot on top of the PSVD. One or more of c('Beard',
#'                     'Atlas', Uplinger','VanDijk'). Defaults to 'Beard'.
#' @param contour      Logical: should 2d density estimate contour lines be
#'                     added to the plot? Defaults to FALSE.
#' @param theme        Character vector indicating a plotting theme to use.
#'                     Current options are 'color' (default) or 'bw' (black and
#'                     white).
#' @param outlier      A value between 0 and 1. Removes outlier bins, i.e.
#'                     those that are between (1-value) and (1+value) far from
#'                     the Bear theoretical fall velocity model. Defaults to
#'                     Inf (no outliers are removed).
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
#' dsd_plot(day, model='Beard')
#' # no model
#' dsd_plot(day, model=NA)
#' # no model, add contour lines
#' dsd_plot(day, model=NA, contour=TRUE)
#'
#' filter <- psvd_filter(type='Thies', d=c(0.3,7), tau=0.5)
#' 
#' @export
dsd_plot <- function(x, type='Thies',
                     #model='Beard',
                     contour=FALSE, theme='color', filter=NULL) {

  if (type!='Thies' & type!='Parsivel')
    stop('type must be one of c(Thies, Parsivel)')

  # particle size and velocity bin limits
  dia <- switch(type, Thies=disdRo:::dia_t, Parsivel=disdRo:::dia_p)
  vel <- switch(type, Thies=disdRo:::vel_t, Parsivel=disdRo:::vel_p)
  
  # particle size and velocity bin means
  #dia_m <- switch(type, Thies=disdRo:::dia_m_t, Parsivel=disdRo:::dia_m_p)
  #vel_m <- switch(type, Thies=disdRo:::vel_m_t, Parsivel=disdRo:::vel_m_p)

  # particle size and velocity bin widths
  dia_w <- switch(type, Thies=disdRo:::dia_w_t, Parsivel=disdRo:::dia_w_p)
  vel_w <- switch(type, Thies=disdRo:::vel_w_t, Parsivel=disdRo:::vel_w_p)
  
  # transform data to long format for ggplot
  x_m <- reshape2::melt(x)
  xx <- x
  dimnames(xx) <- list(dia_w, vel_w)
  x_m <- cbind(x_m,
               reshape2::melt(xx))[,-3]
  colnames(x_m) <- c('dia','vel','dia_w','vel_w','n')
  x_m$n <- log10(x_m$n)
  
  # change opacity of outliers according to filter
  if (!is.null(filter)) {
    x_m$filter <- reshape2::melt(filter)[,3]
  } else {
    x_m$filter <- TRUE
  }
  x_m$alp <- x_m$filter
  x_m$alp <- pmin(1, (x_m$alp+0.75))
  
  x_m <- x_m[x_m$n>1,]
  #summary(x_m)
    
  # heatmap - rectangular binding
  # original version, casts a warning in R CMD CHECK: g <- ggplot(x_m, aes(x=dia, y=vel, width=dia_w, height=vel_w, fill=n)) +
  g <- ggplot(x_m, aes_(x=~dia, y=~vel, width=~dia_w, height=~vel_w, fill=~n, alpha=~alp)) +
    geom_tile() +
    xlim(c(0,8)) + ylim(c(0,15)) +
    xlab('Diameter (mm)') + ylab('Velocity (m/s)') +
    theme_bw() + 
    guides(alpha=FALSE) +
    theme(panel.grid=element_blank())
  if (theme=='color') {
    g <- g + scale_fill_gradient2(name='NS', low='darkgreen', mid='yellow', high='darkred',
                         limits=c(0,5), midpoint=2,
                         na.value='darkgreen',labels=c(0,10,100,1000,10000,100000))
    if (length(model)>0) {
      if ('Beard' %in% model) {
        g <- g + stat_function(aes(col='Beard'), fun=psvd_model,
                               args=list(model='Beard'))
      }
      if ('Atlas' %in% model) {
        g <- g + stat_function(aes(col='Atlas'), fun=psvd_model,
                               args=list(model='Atlas'))
      }
      if ('Uplinger' %in% model) {
        g <- g + stat_function(aes(col='Uplinger'), fun=psvd_model,
                               args=list(model='Uplinger'))
      }
      if ('VanDijk' %in% model) {
        g <- g + stat_function(aes(col='VanDijk'), fun=psvd_model,
                               args=list(model='VanDijk'))
      }
      g <- g + guides(col=guide_legend(title='Model'))
      if (contour) {
        g <- g + geom_density_2d(alpha=0.5, col='dark grey')
        #g <- g + stat_density_2d(aes(fill = ..level..), geom = "polygon")
      }
    }
  }
  if (theme=='bw') {
    g <- g + scale_fill_gradient2(name='NS', low='white',mid='grey50',high='black',
                                  limits=c(0,5), midpoint=2,
                                  na.value='white',labels=c(0,10,100,1000,10000,100000))
    if (length(model)>0) {
      if ('Beard' %in% model) {
        g <- g + stat_function(aes(linetype='Beard'), fun=psvd_model,
                               args=list(model='Beard'))
      }
      if ('Atlas' %in% model) {
        g <- g + stat_function(aes(linetype='Atlas'), fun=psvd_model,
                               args=list(model='Atlas'))
      }
      if ('Uplinger' %in% model) {
        g <- g + stat_function(aes(linetype='Uplinger'), fun=psvd_model,
                               args=list(model='Uplinger'))
      }
      if ('VanDijk' %in% model) {
        g <- g + stat_function(aes(linetype='VanDijk'), fun=psvd_model,
                               args=list(model='VanDijk'))
      }
      g <- g + guides(col=guide_legend(title='Model'))
      if (contour) {
        g <- g + geom_density_2d(alpha=0.5, col='dark grey')
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
#' 
#' @section Note
#' It is assumed that the data consists on the complete telegram is recorded,
#' which may not correspond to the default factory settings. Custom definition
#' of the telegram needs to be implemented.
#' 
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
