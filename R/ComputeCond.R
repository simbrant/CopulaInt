compute.cond.y  <- function(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.type.y="c",var.types.x,sim=TRUE,n.sim=1000,delta=1e-3)
 {
   if(var.type.y == "d")
    {
      E.y.cond <- compute.cond.discr.y(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.types.x)
    }
   else
    {
      if(sim)
       {
      	 E.y.cond <-  compute.cond.cont.y.sim(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.types.x,n.sim=n.sim)
       }
      else
       {
      	 E.y.cond <- compute.cond.cont.y.int(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.types.x,delta=delta)
       }
     }

    E.y.cond
 }

compute.cond.discr.y  <- function(x, Matrix.xy, pair.copulas.xy, F.x, F.y, var.types.x)
 {
    d <- ncol(x)
    n <- nrow(x)
    Matrix.x <- Matrix.xy[2:(d+1),2:(d+1)]
    pair.copulas.x <- pair.copulas.xy
    for(i in 1:(d-1))
     {
       pair.copulas.x[[i]][[1]] <- NULL
     }
     pair.copulas.x[[d]] <- NULL
     cond.x <- compute.cond.x(x,Matrix.x,pair.copulas.x,F.x,var.types.x)
     ind.cond.x <- cond.x$ind.cond.x
     u.cond.x <- cond.x$u.cond.x
     ind.cond.x.discr <- cond.x$ind.cond.x.discr
     u.cond.x.discr <- cond.x$u.cond.x.discr
     top.nodes.x <- Matrix.x[1:2,1]
     ind.x <- which(top.nodes.x == Matrix.xy[2,1])
     ind.cond <- nrow(ind.cond.x)-(1:0)[ind.x]
     if(length(ind.cond.x.discr) > 0)
     {
        ind.cond.discr <-which((ind.cond.x.discr[,1]==ind.cond.x[ind.cond,1])&(apply(matrix(ind.cond.x.discr[,2:d],ncol=d-1),1,setequal,ind.cond.x[ind.cond,2:d])))
      }
     y.val <- sort(unique(F.y[[2]]))
     u.y <- F.y[[1]](y.val)
     u.x.tmp <- cbind(rep(u.y[1],n),u.cond.x[,ind.cond],rep(0,n))
     if(length(ind.cond.x.discr) > 0)
     {
        u.x.tmp <- cbind(u.x.tmp,u.cond.x.discr[,ind.cond.discr])
      }
     p.y.cond <-  matrix(rvinecopulib::hbicop(u.x.tmp,2,pair.copulas.xy[[d]][[1]]$family,pair.copulas.xy[[d]][[1]]$rotation,pair.copulas.xy[[d]][[1]]$par,var_types=c("d",var.types.x[Matrix.xy[2,1]])),ncol=1)
     for(i in 2:length(y.val))
      {
     	u.x.tmp <- cbind(rep(u.y[i],n),u.cond.x[,ind.cond],rep(u.y[i-1],n))
	if(length(ind.cond.x.discr) > 0)
          {
            u.x.tmp <- cbind(u.x.tmp,u.cond.x.discr[,ind.cond.discr])
	  }
        p.y.cond <-  cbind(p.y.cond,rvinecopulib::hbicop(u.x.tmp,2,pair.copulas.xy[[d]][[1]]$family,pair.copulas.xy[[d]][[1]]$rotation,pair.copulas.xy[[d]][[1]]$par,var_types=c("d",var.types.x[Matrix.xy[2,1]]))-p.y.cond[,i-1])
      }
    E.y.cond <- apply(t(t(p.y.cond)*y.val),1,sum)

    E.y.cond
 }

 compute.cond.cont.y.int <- function(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.types.x,delta=1e-3)
 {
    d <- ncol(x)
    n <- nrow(x)
    u.x.tmp <- compute.cond.y.prep(x,Matrix.xy,pair.copulas.xy,F.x,var.types.x,inv=FALSE)
    u <- seq(delta,1-delta,delta)
    y.u <- quantile(F.y[[2]],probs=u)
    f.cond.y.tmp <- matrix(compute.cond.density.y(rep(u,each=nrow(u.x.tmp)),matrix(rep(t(u.x.tmp),length(u)),ncol=ncol(u.x.tmp),byrow=TRUE),Matrix.xy,pair.copulas.xy,var.types.x),nrow=nrow(u.x.tmp))
    f.cond.y <- matrix(rep(y.u,each=nrow(u.x.tmp))*compute.cond.density.y(rep(u,each=nrow(u.x.tmp)),matrix(rep(t(u.x.tmp),length(u)),ncol=ncol(u.x.tmp),byrow=TRUE),Matrix.xy,pair.copulas.xy,var.types.x),nrow=nrow(u.x.tmp))
    E.y.cond <- apply(0.5*delta*(f.cond.y[,-length(u)]+f.cond.y[,-1]),1,sum)

    E.y.cond
  }

compute.cond.cont.y.sim <- function(x,Matrix.xy,pair.copulas.xy,F.x,F.y,var.types.x,n.sim=1000)
 {
   d <- ncol(x)
   n <- nrow(x)
   u.x.tmp <- compute.cond.y.prep(x,Matrix.xy,pair.copulas.xy,F.x,var.types.x)
   u.y <- runif(n*n.sim)
   u.y.inv <-  compute.cond.inv.y(u.y,
                                  matrix(rep(t(u.x.tmp),n.sim),
                                         ncol=ncol(u.x.tmp),
                                         byrow=TRUE),
                                  var.types.x,
                                  Matrix.xy,
                                  pair.copulas.xy)
   y.cond <- matrix(quantile(F.y[[2]],probs=u.y.inv),ncol=n.sim)
   E.y.cond <- apply(y.cond,1,mean)

   E.y.cond
 }

compute.cond.density.y <- function(u.y,u.x,Matrix.xy,pair.copulas.xy,var.types.x)
 {
 	d <- ncol(Matrix.xy)-1
 	u.x <- matrix(u.x,nrow=length(u.y))
 	ind.discr <- which(var.types.x == "d")
    is.discr <- as.numeric(var.types.x[Matrix.xy[d+1,1]] == "d")
    k <- 1
    density.cond.y <- h.func.deriv(u.y,u.x[,k],pair.copulas.xy[[1]][[1]]$family,pair.copulas.xy[[1]][[1]]$rotation,pair.copulas.xy[[1]][[1]]$par,var.types.x[Matrix.xy[d+1,1]],u.minus=u.x[,k+is.discr])
    for(j in 2:d)
     {
     	 u.y <- h.func(u.y,u.x[,k],2,pair.copulas.xy[[j-1]][[1]]$family,pair.copulas.xy[[j-1]][[1]]$rotation,pair.copulas.xy[[j-1]][[1]]$par,var.types.x[Matrix.xy[d-j+3,1]],u.minus=u.x[,k+is.discr])
     	 k <- k+1+is.discr
         is.discr <- as.numeric(var.types.x[Matrix.xy[d-j+2,1]] == "d")
     	 density.cond.y <- density.cond.y*h.func.deriv(u.y,u.x[,k],pair.copulas.xy[[j]][[1]]$family,pair.copulas.xy[[j]][[1]]$rotation,pair.copulas.xy[[j]][[1]]$par,var.types.x[Matrix.xy[d-j+2,1]],u.minus=u.x[,k+is.discr])
     }

    density.cond.y
 }

compute.cond.y.prep <- function(x,Matrix.xy,pair.copulas.xy,F.x,var.types.x,inv=TRUE)
 {
    d <- ncol(x)
    ind.discr <- which(var.types.x == "d")
    u.x <- F.all(x, var.types.x, F.x)$u.x
    Matrix.x <- Matrix.xy[2:(d+1),2:(d+1)]
    pair.copulas.x <- pair.copulas.xy
    for(i in 1:(d-1))
     {
        pair.copulas.x[[i]][[1]] <- NULL
     }
    pair.copulas.x[[d]] <- NULL
    cond.x <- compute.cond.x(x,Matrix.x,pair.copulas.x, F.x, var.types.x)
    ind.cond.x <- cond.x$ind.cond.x
    u.cond.x <- cond.x$u.cond.x
    ind.cond.x.discr <- cond.x$ind.cond.x.discr
    u.cond.x.discr <- cond.x$u.cond.x.discr
    top.nodes.x <- Matrix.x[1:2,1]
    ind.x <- nrow(ind.cond.x)-(1:0)[which(top.nodes.x == Matrix.xy[2,1])]
    u.x.tmp <- u.cond.x[,ind.x]
    if(length(ind.cond.x.discr) > 0)
     {
        ind.x.discr <-which((ind.cond.x.discr[,1]==ind.cond.x[ind.x,1])&(apply(matrix(ind.cond.x.discr[,2:d],ncol=d-1),1,setequal,ind.cond.x[ind.x,2:d])))
        u.x.tmp <- cbind(u.x.tmp,u.cond.x.discr[,ind.x.discr])
     }
    for(j in 1:(d-1))
     {
      	if(j == d-1)
      	 {
      	  	ind.x <- Matrix.xy[d+1,1]
      	  	if(inv)
      	  	 {
      	  	 	u.x.tmp <- cbind(u.x.tmp,u.x[,ind.x])
      	  	 }
      	  	else
      	  	 {
      	  	 	 u.x.tmp <- cbind(u.x[,ind.x],u.x.tmp)
      	  	 }
      	  	if(length(ind.discr) > 0)
      	  	 {
      	  	   	if(inv)
      	  	     {
      	  	 	    u.x.tmp <- cbind(u.x.tmp,u.x[,d+which(ind.discr == ind.x)])
      	  	     }
      	  	   else
      	  	    {
      	  	 	    u.x.tmp <- cbind(u.x.tmp[,1],u.x[,d+which(ind.discr == ind.x)],u.x.tmp[,2:ncol(u.x.tmp)])
      	  	    }

      	  	 }
      	  }
      	 else
      	  {
      	  	  ind.x <- which((ind.cond.x[,1]==Matrix.xy[j+2,1])&(apply(matrix(ind.cond.x[,2:(d-j)],ncol=d-j-1),1,setequal,Matrix.xy[(j+3):(d+1),1])))
      	  	  if(inv)
      	  	   {
      	  	 	  u.x.tmp <- cbind(u.x.tmp,u.cond.x[,ind.x])
      	  	   }
      	  	  else
      	  	   {
      	  	 	  u.x.tmp <- cbind(u.cond.x[,ind.x],u.x.tmp)
      	  	   }
              if(length(ind.cond.x.discr) > 0)
               {
                   ind.x.discr <- which((ind.cond.x.discr[,1]==ind.cond.x[ind.x,1])&(apply(matrix(ind.cond.x.discr[,2:(d-j)],ncol=d-j-1),1,setequal,ind.cond.x[ind.x,2:(d-j)])))
                   if(inv)
      	  	        {
      	  	 	       u.x.tmp <- cbind(u.x.tmp,u.cond.x.discr[,ind.x.discr])
      	  	         }
      	  	       else
      	  	        {
      	  	 	       u.x.tmp <- cbind(u.x.tmp[,1],u.cond.x.discr[,ind.x.discr],u.x.tmp[,2:ncol(u.x.tmp)])
      	  	        }
                }
            }
        }

     u.x.tmp
 }

compute.cond.inv.y <- function(u.y,u.x,var.types.x, Matrix.xy, pair.copulas.xy)
 {
    d <- length(var.types.x)
    k <- 0
    u.y.inv <- u.y
    for(j in 1:d)
     {
         u.yx <- cbind(u.y.inv,u.x[,k+1])
         k <- k+1
         if(var.types.x[Matrix.xy[j+1,1]] == "d")
          {
              u.yx <- cbind(u.yx,u.x[,k+1])
              k <- k+1
          }
         u.y.inv <- rvinecopulib::hbicop(u.yx, 2,
                                         pair.copulas.xy[[d-j+1]][[1]]$family,
                                         pair.copulas.xy[[d-j+1]][[1]]$rotation,
                                         pair.copulas.xy[[d-j+1]][[1]]$par,
                                         inverse=TRUE,
                                         var_types=c("c", var.types.x[Matrix.xy[j+1,1]]))
    }

  u.y.inv
 }

compute.cond.x <- function(x,Matrix.x,pair.copulas.x,F.x,var.types.x)
 {
   d <- ncol(x)
   ind.discr <- which(var.types.x == "d")
   u.x <- F.all(x, var.types.x, F.x)$u.x

   ind.cond.x <- matrix(0,d*(d-1),d)
   ind.cond.x.discr <- c()
   u.cond.x <- c()
   u.cond.x.discr <- c()
   ind.cond.discr.1 <- ind.cond.discr.2 <- NULL
   for(i in 1:(d-1))
    {
      ind.cond.x[2*(i-1)+1:2,1] <- c(Matrix.x[i,i],Matrix.x[d,i])
      ind.cond.x[2*(i-1)+1:2,2] <- c(Matrix.x[d,i],Matrix.x[i,i])
      u.x.tmp <- u.x[,c(Matrix.x[i,i],Matrix.x[d,i])]
      u.cond.x.1 <- h.func(u.x.tmp[,1],u.x.tmp[,2],2,pair.copulas.x[[1]][[i]]$family,pair.copulas.x[[1]][[i]]$rotation,pair.copulas.x[[1]][[i]]$par,var.type=var.types.x[Matrix.x[d,i]],u.minus=u.x[,d+which(ind.discr == Matrix.x[d,i])])
      u.cond.x.2 <- h.func(u.x.tmp[,1],u.x.tmp[,2],1,pair.copulas.x[[1]][[i]]$family,pair.copulas.x[[1]][[i]]$rotation,pair.copulas.x[[1]][[i]]$par,var.type=var.types.x[Matrix.x[i,i]],u.minus=u.x[,d+which(ind.discr == Matrix.x[i,i])])
      u.cond.x <- cbind(u.cond.x,u.cond.x.1,u.cond.x.2)
      if(var.types.x[Matrix.x[i,i]] == "d")
       {
         ind.cond.x.discr <- rbind(ind.cond.x.discr,ind.cond.x[2*(i-1)+1,])
	     u.x.tmp <- u.x[,c(d+which(ind.discr == Matrix.x[i,i]),Matrix.x[d,i])]
         u.cond.x.discr.1 <- h.func(u.x.tmp[,1],u.x.tmp[,2],2,pair.copulas.x[[1]][[i]]$family,pair.copulas.x[[1]][[i]]$rotation,pair.copulas.x[[1]][[i]]$par,var.type=var.types.x[Matrix.x[d,i]],u.minus=u.x[,d+which(ind.discr == Matrix.x[d,i])])
	    u.cond.x.discr <- cbind(u.cond.x.discr,u.cond.x.discr.1)
       }
      if(var.types.x[Matrix.x[d,i]] == "d")
       {
         ind.cond.x.discr <- rbind(ind.cond.x.discr,ind.cond.x[2*(i-1)+2,])
	     u.x.tmp <- u.x[,c(Matrix.x[i,i],d+which(ind.discr == Matrix.x[d,i]))]
         u.cond.x.discr.2 <- h.func(u.x.tmp[,1],u.x.tmp[,2],1,pair.copulas.x[[1]][[i]]$family,pair.copulas.x[[1]][[i]]$rotation,pair.copulas.x[[1]][[i]]$par,var.type=var.types.x[Matrix.x[i,i]],u.minus=u.x[,d+which(ind.discr == Matrix.x[i,i])])
	 u.cond.x.discr <- cbind(u.cond.x.discr,u.cond.x.discr.2)
       }
   }
   if(d > 2)
    {
      for(j in 2:(d-1))
       {
         for(i in 1:(d-j))
          {
            ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1:2,1] <- c(Matrix.x[i,i],Matrix.x[d-j+1,i])
            ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1:2,2:(j+1)] <- matrix(c(Matrix.x[(d-j+1):d,i],Matrix.x[c(i,(d-j+2):d),i]),ncol=j,byrow=TRUE)
            ind.cond.1 <- which((ind.cond.x[,1]==ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1,1])&(apply(matrix(ind.cond.x[,2:j],ncol=j-2+1),1,setequal,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1,3:(j+1)])))
            ind.cond.2 <- which((ind.cond.x[,1]==ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+2,1])&(apply(matrix(ind.cond.x[,2:j],ncol=j-2+1),1,setequal,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+2,3:(j+1)])))
            if(length(ind.cond.x.discr) > 0)
             {
	       ind.cond.discr.1 <- which((ind.cond.x.discr[,1]==ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1,1])&(apply(matrix(ind.cond.x.discr[,2:j],ncol=j-2+1),1,setequal,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1,3:(j+1)])))
               ind.cond.discr.2 <- which((ind.cond.x.discr[,1]==ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+2,1])&(apply(matrix(ind.cond.x.discr[,2:j],ncol=j-2+1),1,setequal,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+2,3:(j+1)])))
             }
	    u.x.tmp <- u.cond.x[,c(ind.cond.1,ind.cond.2)]
            u.cond.x.1 <- h.func(u.x.tmp[,1],u.x.tmp[,2],2,pair.copulas.x[[j]][[i]]$family,pair.copulas.x[[j]][[i]]$rotation,pair.copulas.x[[j]][[i]]$par,var.type=var.types.x[Matrix.x[d-j+1,i]],u.minus=u.cond.x.discr[,ind.cond.discr.2])
            u.cond.x.2 <- h.func(u.x.tmp[,1],u.x.tmp[,2],1,pair.copulas.x[[j]][[i]]$family,pair.copulas.x[[j]][[i]]$rotation,pair.copulas.x[[j]][[i]]$par,var.type=var.types.x[Matrix.x[i,i]],u.minus=u.cond.x.discr[,ind.cond.discr.1])
            u.cond.x <- cbind(u.cond.x,u.cond.x.1,u.cond.x.2)
            if(var.types.x[Matrix.x[i,i]] == "d")
       	     {
               ind.cond.x.discr <- rbind(ind.cond.x.discr,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+1,])
	       u.x.tmp <- cbind(u.cond.x.discr[,ind.cond.discr.1],u.cond.x[,ind.cond.2])
               u.cond.x.discr.1 <- h.func(u.x.tmp[,1],u.x.tmp[,2],2,pair.copulas.x[[j]][[i]]$family,pair.copulas.x[[j]][[i]]$rotation,pair.copulas.x[[j]][[i]]$par,var.type=var.types.x[Matrix.x[d-j+1,i]],u.minus=u.cond.x.discr[,ind.cond.discr.1])
	       u.cond.x.discr <- cbind(u.cond.x.discr,u.cond.x.discr.1)
       	     }
      	    if(var.types.x[Matrix.x[d-j+1,i]] == "d")
       	     {
               ind.cond.x.discr <- rbind(ind.cond.x.discr,ind.cond.x[(j-1)*(2*d-j)+2*(i-1)+2,])
	       u.x.tmp <- cbind(u.cond.x[,ind.cond.1],u.cond.x.discr[,ind.cond.discr.2])
               u.cond.x.discr.2 <- h.func(u.x.tmp[,1],u.x.tmp[,2],1,pair.copulas.x[[j]][[i]]$family,pair.copulas.x[[j]][[i]]$rotation,pair.copulas.x[[j]][[i]]$par,var.type=var.types.x[Matrix.x[i,i]],u.minus=u.cond.x.discr[,ind.cond.discr.2])
	       u.cond.x.discr <- cbind(u.cond.x.discr,u.cond.x.discr.2)
       	     }
          }
       }
    }

list(ind.cond.x=ind.cond.x,u.cond.x=u.cond.x,ind.cond.x.discr=ind.cond.x.discr,u.cond.x.discr=u.cond.x.discr)
 }

h.func <- function(u1,u2,cond.var,family,rotation,par,var.type="c",u.minus=NULL)
 {
   u1 <- pmin(pmax(u1,0),1)
   u2 <- pmin(pmax(u2,0),1)
   u <- cbind(u1,u2)
   if(var.type == "d")
    {
      if(cond.var == 1)
       {
         u.tmp <- cbind(u.minus,u[,2])
       }
      else
       {
         u.tmp <- cbind(u[,1],u.minus)
       }
      h.func <- (rvinecopulib::pbicop(u,family,rotation,par)-rvinecopulib::pbicop(u.tmp,family,rotation,par))/(u[,cond.var]-u.minus)
    }
   else
    {
      h.func <- rvinecopulib::hbicop(u,cond.var,family,rotation,par)
    }

   pmin(pmax(h.func,0),1)
 }

 h.func.deriv <- function(u1,u2,family,rotation,par,var.type="c",u.minus=NULL)
 {
   u1 <- pmin(pmax(u1,0),1)
   u2 <- pmin(pmax(u2,0),1)
   u <- cbind(u1,u2)
   if(var.type == "d")
    {
       u.tmp <- cbind(u[,1],u.minus)
       h.func.deriv <- (h.func(u[,1],u[,2],1,family,rotation,par)-h.func(u.tmp[,1],u.tmp[,2],1,family,rotation,par))/(u[,2]-u.minus)
    }
   else
    {
      h.func.deriv <- rvinecopulib::dbicop(u,family,rotation,par)
    }

   h.func.deriv
 }

