# Dimension function which works for not only matrix and data.frame
dimension <- function( ... ){
	args <- list(...)
	lapply(args, function(x){
			   if(is.null(dim(x)))
				  return(length(x))
				  dim(x)})
}
