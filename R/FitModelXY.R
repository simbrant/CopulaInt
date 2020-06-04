
fit.model.xy <- function(y, x, var.type.y, var.types.x,family.set=c("gaussian","clayton","gumbel"), x.ind=NULL)
{
  d <- ncol(x)
  n <- length(y)
  y <- matrix(y, ncol=1)
  if(length(x.ind) == 0)
  {
    x.ind <- 1:n
  }
  ind.discr.x <- which(var.types.x == "d")
  F.u.x <- F.all(x, var.types.x=var.types.x)
  F.x <- F.u.x$F.x
  u.x <- F.u.x$u.x
  mod.x <- fit.model.x(x,F.x,var.types.x,family.set=family.set)
  Matrix.x <- mod.x$Matrix.x
  pair.copulas.x <- mod.x$pair.copulas.x
  cond.x <- mod.x$cond.x
  cond.x$u.cond.x <- cond.x$u.cond.x[x.ind,]
  if(any(var.types.x == "d"))
  {
    cond.x$u.cond.x.discr <- cond.x$u.cond.x.discr[x.ind,]
  }
  u.x <- u.x[x.ind,]
  mod.xy <- fit.model.y(y,x,u.x,var.type.y,var.types.x,Matrix.x,pair.copulas.x,cond.x,family.set=c("gaussian","clayton","gumbel"))
  Matrix.xy <- mod.xy$Matrix.xy
  pair.copulas.xy <- mod.xy$pair.copulas.xy
  F.y <- mod.xy$F.y

  list(Matrix.xy=Matrix.xy,pair.copulas.xy=pair.copulas.xy,F.x=F.x,F.y=F.y)
}

fit.model.y <- function(y, x, u.x, var.type.y, var.types.x,Matrix.x, pair.copulas.x,cond.x,family.set=c("gaussian","clayton","gumbel"))
{
  d <- ncol(x)
  ind.discr.x <- which(var.types.x == "d")
  ind.cond.x <- cond.x$ind.cond.x
  u.cond.x <- cond.x$u.cond.x
  ind.cond.x.discr <- cond.x$ind.cond.x.discr
  u.cond.x.discr <- cond.x$u.cond.x.discr
  Matrix.xy <- cbind(rep(0,d+1),rbind(rep(0,d),Matrix.x))
  pair.copulas.xy <- vector("list",length=d)
  top.node <- Matrix.x[1:2,1]
  tau.1 <- cor(y,x[,top.node[1]],method="kendall")
  tau.2 <- cor(y,x[,top.node[2]],method="kendall")
  if(abs(tau.1) > abs(tau.2))
  {
    Matrix.xy[,1] <- c(d+1,Matrix.xy[3:(d+1),2],Matrix.xy[2,2])
  }
  else
  {
    if(Matrix.xy[3,3] == top.node[2])
    {
      Matrix.xy[,1] <- c(d+1,top.node[1],Matrix.xy[4:(d+1),3],top.node[2])
    }
    else
    {
      Matrix.xy[,1] <- c(d+1,c(setdiff(cbind(diag(Matrix.xy)[-c(1,d+1)],diag(Matrix.xy[3:(d+1),2:d])),top.node[2]),top.node[2]))
    }
  }
  ind.discr.y <- 2*as.numeric(var.type.y == "d")
  ind.cond.discr <- NULL
  F.u.y <- F.all(y, var.type.y)
  F.y <- list(F.u.y$F.x[[1]],y)
  u.y <- F.u.y$u.x
  u.xy <- cbind(u.y[,1],u.x[,Matrix.xy[d+1,1]],u.y[,ind.discr.y],u.x[,d++which(ind.discr.x == Matrix.xy[d+1,1])])
  cop.1 <- rvinecopulib::bicop(u.xy,c(var.type.y,var.types.x[Matrix.xy[d+1,1]]),family.set)
  pair.copulas.xy[[1]] <- append(pair.copulas.x[[1]],list(cop.1),0)
  u.y <- matrix(h.func(u.xy[,1],u.xy[,2],2,pair.copulas.xy[[1]][[1]]$family,pair.copulas.xy[[1]][[1]]$rotation,pair.copulas.xy[[1]][[1]]$par,var.type=var.types.x[Matrix.xy[d+1,1]],u.minus=u.x[,d+which(ind.discr.x == Matrix.xy[d+1,1])]),ncol=1)
  if(var.type.y == "d")
  {
    u.y.minus <- h.func(u.xy[,3],u.xy[,2],2,pair.copulas.xy[[1]][[1]]$family,pair.copulas.xy[[1]][[1]]$rotation,pair.copulas.xy[[1]][[1]]$par,var.type=var.types.x[Matrix.xy[d+1,1]],u.minus=u.x[,d+which(ind.discr.x == Matrix.xy[d+1,1])])
    u.y <- cbind(u.y,u.y.minus)
  }
  ind.cond <- which((ind.cond.x[,1] == Matrix.xy[d,1])&(apply(matrix(ind.cond.x[,2],ncol=1),1,setequal,Matrix.xy[d+1,1])))
  if(length(ind.cond.x.discr) > 0)
  {
    ind.cond.discr <- which((ind.cond.x.discr[,1]==ind.cond.x[ind.cond,1])&(apply(matrix(ind.cond.x.discr[,2],ncol=1),1,setequal,ind.cond.x[ind.cond,2])))
  }
  u.xy <- cbind(u.y[,1],u.cond.x[,ind.cond],u.y[,ind.discr.y],u.cond.x.discr[,ind.cond.discr])
  for(j in 2:d)
  {
    cop.j <- rvinecopulib::bicop(u.xy,c(var.type.y,var.types.x[Matrix.xy[d-j+2,1]]),family.set)
    if(j < d)
    {
      pair.copulas.xy[[j]] <- append(pair.copulas.x[[j]],list(cop.j),0)
      u.y <- matrix(h.func(u.xy[,1],u.xy[,2],2,pair.copulas.xy[[j]][[1]]$family,pair.copulas.xy[[j]][[1]]$rotation,pair.copulas.xy[[j]][[1]]$par,var.type=var.types.x[Matrix.xy[d-j+2,1]],u.minus=u.cond.x.discr[,ind.cond.discr]),ncol=1)
      if(var.type.y == "d")
      {
        u.y.minus <- h.func(u.xy[,3],u.xy[,2],2,pair.copulas.xy[[j]][[1]]$family,pair.copulas.xy[[j]][[1]]$rotation,pair.copulas.xy[[j]][[1]]$par,var.type=var.types.x[Matrix.xy[d-j+2,1]],u.minus=u.cond.x.discr[,ind.cond.discr])
        u.y <- cbind(u.y,u.y.minus)
      }
      ind.cond <- which((ind.cond.x[,1] == Matrix.xy[d-j+1,1])&(apply(matrix(ind.cond.x[,2:(j+1)],ncol=j),1,setequal,Matrix.xy[(d-j+2):(d+1),1])))
      if(length(ind.cond.x.discr) > 0)
      {
        ind.cond.discr <- which((ind.cond.x.discr[,1]==ind.cond.x[ind.cond,1])&(apply(matrix(ind.cond.x.discr[,2:(j+1)],ncol=j),1,setequal,ind.cond.x[ind.cond,2:(j+1)])))
      }
      u.xy <- cbind(u.y[,1],u.cond.x[,ind.cond],u.y[,ind.discr.y],u.cond.x.discr[,ind.cond.discr])
    }
    else
    {
      pair.copulas.xy[[j]] <- list(cop.j)
    }
  }

  list(Matrix.xy=Matrix.xy,pair.copulas.xy=pair.copulas.xy,F.y=F.y)
}

fit.model.x <- function(x,F.x,var.types.x,family.set=c("gaussian","clayton","gumbel"))
{
  d <- ncol(x)
  ind.discr.x <- which(var.types.x == "d")
  u.x <- F.all(x,var.types.x=var.types.x,F.x=F.x)$u.x
  mod.x <- rvinecopulib::vinecop(data=u.x,var_types=var.types.x,family_set=family.set)
  Matrix.x <- rvinecopulib::as_rvine_matrix(mod.x$structure)[d:1,]
  pair.copulas.x <- mod.x$pair_copulas
  cond.x <- compute.cond.x(x,Matrix.x,pair.copulas.x,F.x,var.types.x)

  list(Matrix.x=Matrix.x,pair.copulas.x=pair.copulas.x,cond.x=cond.x)
}

