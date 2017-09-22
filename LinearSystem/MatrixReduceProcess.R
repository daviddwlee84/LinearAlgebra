library(MASS) # used for fractions() function

# Reduce Augmented Matrix to Strictly Triangular Form
# (only works for square matrix)
# Using Row Operations I, III
ReduceAugmentedMatrix <- function(coefMatrix, attachVector, frac = TRUE, print = FALSE){
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
					# Row Operation I
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
			# Row Operation III
			times = coefMatrix[irow, row]/pivot 
			coefMatrix[irow,] <- coefMatrix[irow,] - times * pivotal_row
			attachVector[irow] <- attachVector[irow] - times * attachVector[row]
		}
		if(print){
			cat("Step ", row, " result:\n")
			cat("pivot: ", pivot, "\n")
			print(coefMatrix)
		}
	}
	answer <- SolveSTF(coefMatrix, attachVector)
	if(frac){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
		answer <- fractions(answer)
	}
	return(list(StrictlyTriangularForm=coefMatrix, AttachVector=matrix(attachVector), Answer=answer))
}

# Solve Strictly Triangular Form by Back Substitution
# (only works for square matrix)
SolveSTF <- function(STFMatrix, attachVector){
	pivot <- STFMatrix[dim(STFMatrix)[1],dim(STFMatrix)[2]]
	answer <- attachVector[dim(STFMatrix)[1]]/pivot
	for(row in c((dim(STFMatrix)[1]-1):1)){
		pivot <- STFMatrix[row,row]
		answer <- c(answer, (attachVector[row] - STFMatrix[row,c(dim(STFMatrix)[1]:(row+1)) ] %*% answer) / pivot)
	}
	answer <- answer[c(length(answer):1)]
	return(answer)
}





# Gaussian Elimination: Reduce Augmented Matrix to Row Echelon Form
# Definition: A matrix is said to be in row echelon form
# (i)	If the first nonzero entry in each nonzero row is 1.
# (ii)	If row k does not consist entirely of zeros, the number of leading zero entries in row k + 1 is greater than the number of leading zero entries in row k.
# (iii)	If there are rows whose entries are all zero, they are below the row shaving nonzero entries.
#
# Using Row Operations I, II, III
GaussianElimination <- function(coefMatrix, attachVector, frac = TRUE, print = FALSE){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NA)
	}
	col <- 1
	for(row in c(1:(dim(coefMatrix)[1]))){
		findPivot <- FALSE
		while(!findPivot){
			if(coefMatrix[row,col] == 0){
				# make sure the pivot element is in pivotal row
				for(irow in c((row+1):dim(coefMatrix)[1])){
					if(coefMatrix[irow,col] != 0){
						# Row Operation I
						# interchange rows => pivotal row is the first row
						if(print){
							cat("Interchange ", irow, "rows with", row, "rows\n")
						}
						coefMatrix[c(row, irow), ] <- coefMatrix[c(irow, row), ]
						attachVector[c(row, irow)] <- attachVector[c(irow, row)]
						findPivot <- TRUE
						break
					}
					if(irow == dim(coefMatrix)[1]){
						# all the possible choices for a pivot element in a given column are 0
						if(print){
							cat("Going to next column (", col, ")\n")
						}
						col <- col + 1
					}
				}
			} else {
				findPivot <- TRUE
			}
		}
		pivot = coefMatrix[row,col]
		pivotal_row = coefMatrix[row,]
		if(row < dim(coefMatrix)[1]){
			for(irow in c((row+1):dim(coefMatrix)[1])){
				# Row Operation III
				times = coefMatrix[irow, col]/pivot 
				coefMatrix[irow,] <- coefMatrix[irow,] - times * pivotal_row
				attachVector[irow] <- attachVector[irow] - times * attachVector[row]
			}
		}

		# Row Operation II
		# it's necessary in order to scale the rows so the leading coefficients are all 1
		coefMatrix[row,] <- pivotal_row/pivot
		attachVector[row] <- attachVector[row]/pivot
		
		if(print){
			cat("Step ", row, " result:\n")
			cat("pivot: ", pivot, "\n")
			print(coefMatrix)
		}
		if(col < dim(coefMatrix)[2]){
			col <- col + 1
		} else {
			break
		}
	}
	if(frac){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
		#answer <- fractions(answer)
	}
	return(list(RowEchelonForm=coefMatrix, AttachVector=matrix(attachVector)))
}

