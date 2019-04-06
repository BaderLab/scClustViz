#' Internal fx to calc relative overlap of labels
#' 
#' @param i index of column of \code{ld}
#' @param ld labelData - rows are xy coords and xy sizes, columns are labels
#' @param oX original x coord
#' @param oY original y coord
#' 

overlap.proportion.point <- function(i,ld,oX,oY) {
  temp <- (ld[c("x","y"),i] - ld[c("x","y"),,drop=F]) / 
    (ld[c("sizeX","sizeY"),i] + ld[c("sizeX","sizeY"),,drop=F])
  temp[,i] <- (ld[c("x","y"),i] - c(oX[i],oY[i])) / ld[c("sizeX","sizeY"),i]
  return(temp)
}


#' Internal fx to calc dist from point to label end
#' 
#' @param i index of column of \code{ld}
#' @param ld labelData - rows are xy coords and xy sizes, columns are labels
#' @param oX original x coord
#' @param oY original y coord
#' 

origin.label.dist <- function(i,ld,oX,oY) {
  tempX <- (ld["x",i] + c(-1,1) * ld["sizeX",i]) - oX[i]
  return(c(x=tempX[which.min(abs(tempX))],ld["y",i] - oY[i]))
}


#' Internal fx to calc overlap of label and plot edge
#' 
#' @param i index of column of \code{ld}
#' @param ld labelData - rows are xy coords and xy sizes, columns are labels
#' 

overlap.edge <- function(i,ld) {
  rbind(ld["x",i] + c(-1,1) * ld["sizeX",i] - par("usr")[1:2],
        ld["y",i] + c(-1,1) * ld["sizeY",i] - par("usr")[3:4])
}


#' Prevent overlapping labels
#'
#' Determines clever placement of labels for plots made in base R
#'
#' @param x The set of x coordinates for the points to plot
#' @param y The set of y coordinates for the points to plot
#' @param label A character vector of labels for each x,y pair.
#' @param padding A numeric value representing the extra proportion of string
#'   size to add around each label. Not clear that this works as intended, so
#'   leave it at zero.
#' @param str.cex The cex value for the label. Defaults to
#'   \code{\link[graphics]{par}("cex")}. Passed to \code{\link[graphics]{strwidth}}.
#' @param str.font The font value for the label. Defaults to
#'   \code{\link[graphics]{par}("font")}. Passed to \code{\link[graphics]{strwidth}}.
#' 
#' @return A matrix with 2 columns corresponding to x and y coordinates of the
#'   label centre to be used downstream as an input e.g. to \code{text}
#' @export
#'
#' @references This code is inspired by
#'   \href{https://github.com/federicomarini/spreadLabels}{Federico Marini} and
#'   \href{https://github.com/ChristofferFlensburg/spreadLabels}{Christoffer
#'   Flensburg}.
#'   

spreadLabels2 <- function(x,y,label,padding=0,str.cex=par("cex"),str.font=par("font")) {
  #starting condition
  label.data <- rbind(x=x,y=y,
                      sizeX=strwidth(label,cex=str.cex,font=str.font) * .5 + 
                        strwidth("o",cex=str.cex,font=str.font) * padding,
                      sizeY=strheight(label,cex=str.cex,font=str.font) * .5 + 
                        strheight("o",cex=str.cex,font=str.font) * padding)
  # ^ distance from centre of label to x,y edge
  
  odp <- sapply(seq_len(ncol(label.data)),overlap.proportion.point,
                ld=label.data,oX=x,oY=y,simplify=F)
  is.overlap.point <- do.call(rbind,sapply(odp,function(X) apply(X,2,function(Y) all(abs(Y) <= 1)),simplify=F))
  ode <- sapply(seq_len(ncol(label.data)),overlap.edge,ld=label.data,simplify=F)
  is.overlap.edge <- sapply(ode,function(X) any(X[,1] <= 0) | any(X[,2] >= 0))
  
  #randomly move away from locally densest area until not overlapping 
  tempITER <- 0
  while (sum(is.overlap.point,is.overlap.edge) > 0 & tempITER < 200) {
    tempITER <- tempITER + 1
    
    # Push away from overlap with other labels
    for (i in which(apply(is.overlap.point,2,any))) {
      overlap.point.centroid <- rowMeans(odp[[i]][,is.overlap.point[,i],drop=F])
      sign.opc <-sign(overlap.point.centroid)
      sign.opc[sign.opc == 0] <- 1
      move.point.dist <- c(0,0)
      move.point.dist[which.min(abs(overlap.point.centroid))] <- 
        (sign.opc * (1 - abs(overlap.point.centroid)))[which.min(abs(overlap.point.centroid))]
      label.data[c("x","y"),i] <- label.data[c("x","y"),i] + 
        move.point.dist * label.data[c("sizeX","sizeY"),i] * rlnorm(2,meanlog=0,sdlog=0.3)
    }
    
    # Pull label end back to point
    od <- sapply(seq_len(ncol(label.data)),origin.label.dist,
                 ld=label.data,oX=x,oY=y)
    is.origin.dist <- apply(od,2,function(X)
      any(abs(X) > c(diff(par("usr")[1:2]) * 0.02,diff(par("usr")[3:4]) * 0.02)))
    for (i in which(is.origin.dist)) {
      label.data[c("x","y"),i] <- label.data[c("x","y"),i] -
        od[,i] * rlnorm(2,meanlog=-.3,sdlog=.1)
    }

    # Push away from overlap with edges
    ode <- sapply(seq_len(ncol(label.data)),overlap.edge,ld=label.data,simplify=F)
    is.overlap.edge <- sapply(ode,function(X) any(X[,1] <= 0) | any(X[,2] >= 0))
    for (i in which(is.overlap.edge)) {
      oe <- ode[[i]]
      oe[oe[,1] > 0,1] <- 0
      oe[oe[,2] < 0,2] <- 0
      move.edge <- c(0,0)
      move.edge <- move.edge - c(oe[1,which.max(abs(oe[1,]))],
                                 oe[2,which.max(abs(oe[2,]))])
      label.data[c("x","y"),i] <- label.data[c("x","y"),i] + 
        move.edge * rlnorm(2,meanlog=0,sdlog=0.5)
    }
    # points(t(label.data[c("x","y"),]),pch=20,col=viridis(100,d=-1,a=.5)[tempITER])
    # message(tempITER)
    
    odp <- sapply(seq_len(ncol(label.data)),overlap.proportion.point,
                  ld=label.data,oX=x,oY=y,simplify=F)
    is.overlap.point <- do.call(rbind,sapply(odp,function(X) apply(X,2,function(Y) all(abs(Y) <= 1)),simplify=F))
    ode <- sapply(seq_len(ncol(label.data)),overlap.edge,ld=label.data,simplify=F)
    is.overlap.edge <- sapply(ode,function(X) any(X[,1] <= 0) | any(X[,2] >= 0))
  }
  return(t(label.data[c("x","y"),]))
    
}
