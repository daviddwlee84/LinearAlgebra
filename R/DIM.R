# Dimension
#
# Dimension function which works for not only matrix and data.frame
# @param ... R objects
# @return dimension of any R object
dimension <- function( ... ){
	args <- list(...)
	lapply(args, function(x){
			   if(is.null(dim(x)))
				  return(length(x))
				  dim(x)})
}
