library(MASS) # used for fractions() function
#' Reduce Augmented Matrix
#' 
#' Reduce Augmented Matrix to Strictly Triangular Form
#' (only works for square matrix)
#' Using Row Operations I, III
#' @param coefMatrix ceofficient matrix
#' @param attachVector An additional column which attach to the coefficient matrix => augmented matrix
#' @param FRAC Flag of showing result by fraction (Default is TRUE)
#' @param PRINT Flag of printing process detail (Default is FALSE)
#' @return Return a List with StrictlyTriangularForm; AttachVector; Answer
#' @export
ReduceAugmentedMatrix <- function(coefMatrix, attachVector, FRAC = TRUE, PRINT = FALSE){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NULL)
	}
	if(dim(coefMatrix)[1] != dim(coefMatrix)[2]){
		cat("Error: Please Input Square Matrix\n")
		return(NULL)
	}
	for(row in c(1:(dim(coefMatrix)[1]-1))){
		if(coefMatrix[row,row] == 0){
			# make sure the pivot element is in pivotal row
			for(irow in c((row+1):dim(coefMatrix)[1])){
				if(coefMatrix[irow,row] != 0){
					# Row Operation I
					# interchange rows => pivotal row is the first row
					if(PRINT){
						cat("Interchange ", irow, "rows with", row, "rows\n")
					}
					coefMatrix[c(row, irow),] <- coefMatrix[c(irow, row),]
					attachVector[c(row, irow)] <- attachVector[c(irow, row)]
					break
				}
				if(irow == dim(coefMatrix)[1]){
					cat("Error: The ", row, " column is all 0\n")
					return(NULL)
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
		if(PRINT){
			cat("Step ", row, " result:\n")
			cat("pivot: ", pivot, "\n")
			print(coefMatrix)
		}
	}
	answer <- SolveSTF(coefMatrix, attachVector)
	if(FRAC){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
		answer <- fractions(answer)
	}
	return(list(StrictlyTriangularForm=coefMatrix, AttachVector=matrix(attachVector), Answer=answer))
}

# Solve Strictly Triangular Form by Back Substitution
# (only works for square matrix)
# TODO: consistent: multiple answer & inconsistent
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





#' Gaussian Elimination
#' 
#' Gaussian Elimination: Reduce Augmented Matrix to Row Echelon Form
#' Definition: A matrix is said to be in row echelon form
#' (i)	If the first nonzero entry in each nonzero row is 1.
#' (ii)	If row k does not consist entirely of zeros, the number of leading zero entries in row k + 1 is greater than the number of leading zero entries in row k.
#' (iii)	If there are rows whose entries are all zero, they are below the row shaving nonzero entries.
#' Using Row Operations I, II, III
#' @param coefMatrix ceofficient matrix
#' @param attachVector An additional column which attach to the coefficient matrix => augmented matrix
#' @param FRAC Flag of showing result by fraction (Default is TRUE)
#' @param PRINT Flag of printing process detail (Default is FALSE)
#' @return Return a List with StrictlyTriangularForm; AttachVector
#' @export
GaussianElimination <- function(coefMatrix, attachVector, FRAC = TRUE, PRINT = FALSE){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NULL)
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
						if(PRINT){
							cat("Interchange ", irow, "rows with", row, "rows\n")
						}
						coefMatrix[c(row, irow), ] <- coefMatrix[c(irow, row), ]
						attachVector[c(row, irow)] <- attachVector[c(irow, row)]
						findPivot <- TRUE
						break
					}
					if(irow == dim(coefMatrix)[1]){
						# all the possible choices for a pivot element in a given column are 0
						if(PRINT){
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
		
		if(PRINT){
			cat("Step ", row, " result:\n")
			print(coefMatrix)
		}
		if(col < dim(coefMatrix)[2]){
			col <- col + 1
		} else {
			break
		}
	}
	answer <- SolveREF(coefMatrix, attachVector)
	if(FRAC){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
		answer <- fractions(answer)
	}
	return(list(RowEchelonForm=coefMatrix, AttachVector=matrix(attachVector), Answer=answer))
}

# Solve Reduce Echelon Form by Back Substitution
# TODO: consistent: multiple answer & inconsistent
SolveREF <- function(REFMatrix, attachVector){
	pivot <- FindNextPivot(REFMatrix, c(1, 1)) # pivot is the location not value
	if(is.null(pivot)){
	  return(NULL)
	} else {
	  leadingVariables <- pivot[2] # leadingVariables store column
	}
	while(TRUE){
	  # finding all leading variable
		pivot <- FindNextPivot(REFMatrix, pivot)
		if(is.null(pivot)){
			break
		}
		leadingVariables <- c(leadingVariables, pivot[2])
	}
	return(answer)
}

# Find Next Pivot (Leading Variable)
FindNextPivot <- function(REFMatrix, currentPivot){
  print(currentPivot)
	if(currentPivot[1] == dim(REFMatrix)[1] || currentPivot[2] == dim(REFMatrix)[2]){
	  return(NULL)
	}
  row <- currentPivot[1] + 1
	for(col in c((currentPivot[2]+1):dim(REFMatrix)[2])){
	  if(REFMatrix[row, col] != 0){
			return(c(row, col))
		}
	  if(col == dim(REFMatrix)[2]) {
  			return(NULL)
		}
	}
}

