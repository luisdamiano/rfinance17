library(gridBase)
library(gridExtra)
library(lattice)
library(HH)

cols = c("#E41F26", "#2EA147", "#1D79B4", "#fdb863", "#E6E6E6", "#cfcfcf", "#63656A")

colMedians <- function(x) {
  f <- median # You can switch to 'mean' if you dare...
  
  if (length(dim(x)) == 4) {
    apply(x, c(1, 4), function(y) { f(y) })
  } else if (length(dim(x)) == 3) {
    apply(x, c(2, 3), function(y) { f(y) })
  } else if (length(dim(x)) == 2) {
    apply(x, 2, function(y) { f(y) } )
  } else {
    f(x)
  }
}

# Data manipulation for plots----------------------------------------------
# if(!exists('SR')) {
#   stop("Not so fast cowboy! This file should be executed following the order of the notebook.")
# }

queue <- list(
  list(
    getMat = function() {
      mat <- do.call(cbind, lapply(1:length(syms), function(i) {
        SR[[i]][[1]]
      }))
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      bquote(atop('Conditional volatility for a single regimen GARCH(1, 1)',
                  '(' * hat(sigma)[nt] * ')'))
    }),
  list(
    getMat = function() {
      mat <- t(colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, , 1, ]))
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      bquote(atop('Conditional volatility for State' ~ 1,
                  '(' * hat(sigma)[n1t] * ')'))
    }),
  list(
    getMat = function() {
      mat <- t(colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, , 2, ]))
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      bquote(atop('Conditional volatility for State' ~ 2,
                  '(' * hat(sigma)[n2t] * ')'))
    }),
  list(
    getMat = function() {
      s1 <- t(colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, , 1, ]))
      s2 <- t(colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, , 2, ]))
      p1 <- t(colMedians(extract(stan.fit, pars = 'F')[[1]][, , ]))
      p2 <- 1 - p1
      mat <- s1 * p1 + s2 * p2
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      bquote(atop('Conditional volatility weighted by filtered state probability',
                  '(' * hat(sigma)[n.t] * ')'))
    }),
  list(
    getMat = function() {
      mat <- t(colMedians(extract(stan.fit, pars = 'F')[[1]][, , ]))
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      k <- 2
      bquote(atop('Filtered probability for State' ~ 1 ~ 'at' ~ t,
                  '(' * hat(p)[n1t] * ')'))
    }),
  list(
    getMat = function() {
      mat <- t(colMedians(extract(stan.fit, pars = 'F')[[1]][, , ]))
      mat <- log(mat / ( 1 - mat))
      colnames(mat) <- syms
      mat
    },
    getMain = function() {
      k <- 2
      bquote(atop('Filtered probability for State' ~ 1 ~ 'at' ~ t,
                  '(' * hat(p)[n1t] * ')'))
    }),
  list(
    getMat = function() {
      n <- 2 # Ford
      sSR <- SR[[n]][[1]] # For single regime
      sRS1 <- colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 1, ])
      sRS2 <- colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 2, ])
      mat <- cbind(sSR, sRS1, sRS2)
      colnames(mat) <- c(bquote(sigma[.(n) ~ .t]), bquote(sigma[.(n) ~ '1'~t]), bquote(sigma[.(n) ~ '2'~t]))
      mat
    },
    getMain = function() {
      n <- 2
      bquote(atop('Conditional volatility for ' ~ .(syms[n]),
                  '(' * hat(sigma)[n.t] * ')'))
    }),
  list(
    getMat = function() {
      n <- 2 # Ford
      sRS1 <- colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 1, ])
      sRS2 <- colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 2, ])
      tsmat <- cbind(sRS1, sRS2)
      colnames(tsmat) <- c('Low', 'High')
#      colnames(tsmat) <- c(bquote(sigma[.(n) ~ .t]), bquote(sigma[.(n) ~ '1'~t]), bquote(sigma[.(n) ~ '2'~t]))
      
      cpmat <- colMedians(extract(stan.fit, pars = 'F')[[1]][, n, ])
      names(cpmat) <- c(bquote(pi[.(n) ~ .t]))
      
      list(tsmat, cpmat)
    },
    getMain = function() {
      n <- 2
      tsmain <- bquote(atop('Conditional volatility for ' ~ .(syms[n]),
                            '(' * hat(sigma)[.(n)~.t] * ')'))
      
      cpmain <- bquote(atop('Filtered probability for State' ~ 1 ~ 'at' ~ t,
                            '(' * hat(p)[.(n)~'1'~t] * ')'))
      
      list(tsmain, cpmain)
    }),
list(
  getMat = function() {
    n <- 2 # Ford
    mat1 <- cbind(x = t(colMedians(extract(stan.fit, pars = 'F')[[1]][, , ]))[, n], y = as.vector(y_t[, n]))
    attr(mat1, 'tag') <- list(bquote(pi[.(n) * '1' * t]), bquote(r_[.(n) * t]))
    
    mat2 <- t(colMedians(extract(stan.fit, pars = 'F')[[1]][, , ]))[, n]
    attr(mat2, 'tag') <- bquote(pi[.(n) * '1' * t])
    
    mat3 <- cbind(
      x = colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 1, ]^2), 
      y = colMedians(extract(stan.fit, pars = 'sigma_t')[[1]][, n, 2, ]^2)
      )
    attr(mat3, 'tag') <- list(bquote(sigma[.(n) * '1' * t]), bquote(sigma^2[.(n) * '2' * t]))

    list(mat1, mat2, mat3)
  },
  getMain = function() {
    n <- 2
    main1 <- bquote(atop('Conditional volatility for ' ~ .(syms[n]),
                          '(' * hat(sigma)[.(n)~.t] * ')'))
    
    main2 <- bquote(atop('Filtered probability for State' ~ 1 ~ 'at' ~ t,
                          '(' * hat(p)[.(n)~'1'~t] * ')'))
    
    main3 <- bquote(atop('Filtered probability for State' ~ 1 ~ 'at' ~ t,
                         '(' * hat(p)[.(n)~'1'~t] * ')'))

    list(main1, main2, main3)
  })
)

# Plots -------------------------------------------------------------------
tscsplot <- function(mat, main = '', mylim = c(floor(min(mat, 0)), ceiling(max(mat))), ...) {
  # Time series
  tsplot <- xyplot.ts(mat,
                      main = main, xlab = list('', cex = 0.8),
                      superpose = TRUE, outside = TRUE, ylim = mylim,
                      scales = list(cex = 0.8, tck = 0.5),
                      key = list(
                        text = list(colnames(mat)),
                        columns = ncol(mat), corner = c(0.5, 0.95), 
                        cex = 1, lines = list(lwd = 2, size = 2),
                        between = 0.2
                      ))
  
  # Cross sectional
  csplot <- splom(mat,
                  xlab = '',
                  superpanel = function(z, ...) {
                    panel.pairs(z, varname.cex = 1.5, axis.text.cex = 0.8, axis.line.tck = 0.5, 
                                prepanel.limits = function(x) { mylim }, ...)
                  },
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, pch = 21, col = cols[6], bg = cols[5], cex = 0.5, ...)
                    panel.abline(lm(y ~ x), lty = 1, lwd = 2, col = cols[4], ...)
                    panel.loess(x, y, span = 1/3, degree = 1, family = 'gaussian',
                                lty = 1, lwd = 2, col = cols[2])
                    panel.text(mylim[2] * 0.5, mylim[2] * 0.95,  cex = 1,
                               eval(parse(text = paste('bquote(hat(rho) == .(sprintf(\'%0.2f\', ', cor(x, y), ')))')))
                    )
                  },
                  key = list(text = list(c('Linear', 'Loess')), columns = 2, 
                             lines = list(col = cols[c(4, 2)]),
                             space = 'top', cex = 1, lwd = 2))
  
  grid.arrange(tsplot, csplot, 
               layout_matrix = matrix(
                 c(rep(1, 35), rep(2, 90)), byrow = TRUE, ncol = 5, nrow = 25)
  )
}

tscpplot <- function(mat, main = '', mylim = c(floor(min(mat[[1]], 0)), ceiling(max(mat[[1]]))), ...) {
  # Time series
  tsmat <- mat[[1]]
  tsplot <- xyplot.ts(tsmat,
                      main = main[[1]], xlab = list('', cex = 0.5),
                      superpose = TRUE, outside = TRUE, ylim = mylim,
                      scales = list(cex = 0.8, tck = 0.5),
                      col = cols[c(2, 1)],
                      key = list(
                        text = list(colnames(tsmat)),
                        columns = ncol(tsmat), corner = c(0.5, 0.95), 
                        cex = 1, lines = list(lwd = 2, size = 2),
                        between = 0.2, col = cols[c(2, 1)]
                      ))
  
  cpmat <- mat[[2]]
  cpplot <- xyplot.ts(cpmat,
                      panel = function(x, y, ...) {
                        panel.xyplot(x, y, col = cols[5], ...)
                        panel.abline(h = 0.5)
                      },
                      main = main[[2]], xlab = list('', cex = 0.8),
                      superpose = TRUE, outside = TRUE, ylim = c(0, 1),
                      scales = list(cex = 0.8, tck = 0.5)
  )
  grid.arrange(tsplot, cpplot, 
               layout_matrix = matrix(
                 1:2, byrow = TRUE, ncol = 2, nrow = 1)
  )
}

cstscpplot <- function(mat, main = '', mylim = c(floor(min(mat[[1]], 0)), ceiling(max(mat[[1]]))), ...) {
  # Cross sectional
  csmat <- mat[[1]]
  csplot1 <- xyplot(x ~ y, data = as.data.frame(csmat),
                  xlab = bquote(.(attr(csmat, 'tag')[[1]])),
                  ylab = bquote(.(attr(csmat, 'tag')[[2]])),
                  main = bquote(.(main[[1]])),
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, pch = 21, col = cols[6], bg = cols[5], cex = 0.5, ...)
                  })
  
  # Time series
  tsmat <- mat[[2]]
  tsplot <- tsacfplots(tsmat, strip = FALSE)
  tsplot <- update(acf.pacf.plot(tsmat)[1], layout=c(1,1), main = main[[2]], strip = FALSE) 
  #   
  #   # xyplot.ts(tsmat,
  #   #                   main = main[[2]], xlab = list('', cex = 0.5),
  #   #                   superpose = TRUE, outside = TRUE,
  #   #                   scales = list(cex = 0.8, tck = 0.5),
  #   #                   col = cols[c(2, 1)],
  #   #                   key = list(
  #   #                     text = list(colnames(tsmat)),
  #   #                     columns = ncol(tsmat), corner = c(0.5, 0.95), 
  #   #                     cex = 1, lines = list(lwd = 2, size = 2),
  #   #                     between = 0.2, col = cols[c(2, 1)]
  #   #                   ))
  
  # Cross sectional
  csmat <- mat[[3]]
  csplot2 <- xyplot(x ~ y, data = as.data.frame(csmat),
                   xlab = bquote(.(attr(csmat, 'tag')[[1]])),
                   ylab = bquote(.(attr(csmat, 'tag')[[2]])),
                   main = bquote(.(main[[3]])),
                   panel = function(x, y, ...) {
                     panel.xyplot(x, y, pch = 21, col = cols[6], bg = cols[5], cex = 0.5, ...)
                   })

  grid.arrange(csplot1, tsplot, csplot2, 
               layout_matrix = matrix(
                 1:3, byrow = TRUE, ncol = 3, nrow = 1)
  )
}