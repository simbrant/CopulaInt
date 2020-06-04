F.all <- function(x,var.types.x,F.x=NULL)
 {
 	d <- ncol(x)
 	n <- nrow(x)
 	if(length(F.x) == 0)
 	 {
 	 	 F.x <- list()
 	 	 for(i in 1:d)
 	 	  {
 	 	  	  F.x[[i]] <- ecdf(x[,i])
 	 	  }
 	 }
   u.x <- x
   u.x.discr <- c()
   for(i in 1:d)
    {
    	u.x[,i] <- F.x[[i]](x[,i])
    	if(var.types.x[i] == "d")
    	 {
    	 	u.x.discr <- cbind(u.x.discr,F.x[[i]](x[,i]-1))
    	 }
    	else
    	 {
    	 	 u.x[,i] <- u.x[,i]*n/(n+1)
    	 }
    }
    u.x <- cbind(u.x,u.x.discr)

    list(u.x=as.matrix(u.x,nrow=nrow(x)),F.x=F.x)
 }
