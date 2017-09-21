# Back Substitution (only works for square matrix)
BackSubstitution <- function(coefMatrix, attachVector, print = FALSE){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NA)
	}
	if(dim(coefMatrix)[1] != dim(coefMatrix)[2]){
		cat("Error: Please Input Square Matrix\n")
		return(NA)
	}
	for(row in c(1:(dim(coefMatrix)[1]-1))){
		if(coefMatrix[row,row] == 0){
			# make sure the pivot element is in pivotal row
			for(irow in c((row+1):dim(coefMatrix)[1])){
				if(coefMatrix[irow,row] != 0){
					# interchange rows => pivotal row is the first row
					if(print){
						cat("Interchange ", irow, "rows with", row, "rows\n")
					}
					coefMatrix[c(row, irow),] <- coefMatrix[c(irow, row),]
					attachVector[c(row, irow)] <- attachVector[c(irow, row)]
					break
				}
				if(irow == dim(coefMatrix)[1]){
					cat("Error: The ", row, " column is all 0\n")
					return(NA)
				}
			}
		}
		pivot = coefMatrix[row,row]
		pivotal_row = coefMatrix[row,]
		for(irow in c((row+1):dim(coefMatrix)[1])){
			times = coefMatrix[irow, row]/pivot 
			coefMatrix[irow,] <- coefMatrix[irow,] - times* pivotal_row
			attachVector[irow] <- attachVector[irow] - times * attachVector[row]
		}
		#coefMatrix[row,] <- pivotal_row/pivot
		#attachVector[row] <- attachVector[row]/pivot
		if(print){
			cat("Step ", row, " result:\n")
			cat("pivot: ", pivot, "\n")
			print(coefMatrix)
		}
	}
	return(list(coefMatrix, matrix(attachVector)))
}

