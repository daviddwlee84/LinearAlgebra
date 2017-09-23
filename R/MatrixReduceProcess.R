library(MASS) # used for fractions() function
#' Reduce Augmented Matrix
#' 
#' Reduce Augmented Matrix to Strictly Triangular Form
#' Using Row Operations I, III
#' (only works for square matrix)
#' (will have the same result as solve(coefMatrix, attachVector))
#' @param coefMatrix ceofficient matrix
#' @param attachVector An additional column which attach to the coefficient matrix => augmented matrix
#' @param FRAC Flag of showing result by fraction (Default is TRUE)
#' @param PRINT Flag of printing process detail (Default is FALSE)
#' @return Return a List with StrictlyTriangularForm; AttachVector; Answer
#' @export
ReduceAugmentedMatrix <- function(coefMatrix, attachVector = NA, FRAC = TRUE, PRINT = FALSE, SOLVE = TRUE){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NULL)
	}
	if(dim(coefMatrix)[1] != dim(coefMatrix)[2]){
		cat("Error: Please Input Square Matrix\n")
		return(NULL)
	}
  if(is.na(attachVector[1])){
    SOLVE <- FALSE
    attachVector <- matrix(0, dim(coefMatrix)[1], 1)
    onlyCoefMatrix <- TRUE
  } else {
    onlyCoefMatrix <- FALSE
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
			if(!onlyCoefMatrix){
			  print(attachVector)			  
			}
		}
	}
  
  if(SOLVE){
    answer <- SolveSTF(coefMatrix, attachVector)
  } else {
    answer <- NULL
  }
  
	if(FRAC){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
		if(SOLVE){
		  answer <- fractions(answer)
		}
	}
  if(!onlyCoefMatrix){
	  return(list(StrictlyTriangularForm=coefMatrix, AttachVector=matrix(attachVector), Answer=answer))
  } else {
    return(StrictlyTriangularForm=coefMatrix)
  }
}

# Solve Strictly Triangular Form by Back Substitution
# (only works for square matrix)
# TODO: infinity solution
SolveSTF <- function(STFMatrix, attachVector){
  if(!isConsistent(STFMatrix, attachVector)){
    cat("System is inconsistent\n")
    return(NULL)
  }
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
#' 
#' @param coefMatrix ceofficient matrix
#' @param attachVector An additional column which attach to the coefficient matrix => augmented matrix
#' @param FRAC Flag of showing result by fraction (Default is TRUE)
#' @param PRINT Flag of printing process detail (Default is FALSE)
#' @param accuracy Number lower than accuracy equal 0 (set NA to unable)
#' @return Return a List with StrictlyTriangularForm; AttachVector; Answer (answer will show in string if there is any free variable)
#' @export
GaussianElimination <- function(coefMatrix, attachVector = NA, FRAC = TRUE, PRINT = FALSE, SOLVE = TRUE, accuracy=1e-10){
	if(is.null(dim(coefMatrix))){
		cat("Error: Please Input Matrix\n")
		return(NULL)
	}
  if(is.na(attachVector[1])){
    attachVector <- matrix(0, dim(coefMatrix)[1], 1)
    SOLVE <- FALSE
    onlyCoefMatrix <- TRUE
  } else {
    onlyCoefMatrix <- FALSE
  }
  
	col <- 1
	for(row in c(1:(dim(coefMatrix)[1]-1))){
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
			if(!onlyCoefMatrix){
			  print(attachVector)
			}
		}
		if(col < dim(coefMatrix)[2]){
			col <- col + 1
		} else {
			break
		}
	}
	
	if(SOLVE){
	  answer <- SolveREF(coefMatrix, attachVector, FRAC, accuracy)
	} else {
	  answer <- NULL
	}
	
	if(FRAC){
		coefMatrix <- fractions(coefMatrix)
		attachVector <- fractions(attachVector)
	}
	if(!onlyCoefMatrix){
	  return(list(RowEchelonForm=coefMatrix, AttachVector=matrix(attachVector), Answer=answer))
	} else {
	  return(RowEchelonForm=coefMatrix)
	}
}

# Solve Reduce Echelon Form by Back Substitution
SolveREF <- function(REFMatrix, attachVector, FRAC, accuracy=NA){
  if(!isConsistent(REFMatrix, attachVector)){
    cat("System is inconsistent\n")
    return(NULL)
  }
  temp <- FindVariables(REFMatrix, accuracy)
  variables <- temp[[1]]
  rows <- temp[[2]]
  varcount <- temp[[3]]
  if(is.null(variables)){
    return(NULL)
  }
  STFMatrix <- matrix(NA, rows, varcount[1])
  RHSMatrix <- matrix(NA, rows, varcount[2]+1) # Right-hand side
  RHSMatrix[,1] <- attachVector[c(1:rows)]
  for(row in c(1:rows)){
    stfcol <- 1
    rhscol <- 2
    for(col in c(1:dim(REFMatrix)[2])){
      switch(variables[col],
             "L"={
               STFMatrix[row, stfcol] <- REFMatrix[row, col]
               stfcol <- stfcol + 1
             },
             "F"={
               RHSMatrix[row, rhscol] <- -REFMatrix[row, col]
               rhscol <- rhscol + 1
             })
    }
  }
  
  answerMatrix <- MSolveSTF(STFMatrix, RHSMatrix)
  
  if(FRAC){
    answerMatrix <- fractions(answerMatrix)
  }
  
  # Find free variables column
  freeVarCol <- NULL
  for(col in c(1:dim(REFMatrix)[2])){
    if(variables[col] == "F"){
      freeVarCol <- c(freeVarCol, col)
    }
  }
  
  if(!is.na(accuracy)){
    for(row in c(1:dim(answerMatrix)[1])){
      for(col in c(1:dim(answerMatrix)[2])){
        if(abs(answerMatrix[row, col]) < accuracy){
          answerMatrix[row, col] <- 0
        }
      }
    }
  }
  
  # Listing the answer
  if(is.null(freeVarCol)){
    # Without free variable
    answer <- answerMatrix
  } else {
    # Show answer in strings and represent free variable in xn
    answer <- NULL
    for(row in c(1:dim(answerMatrix)[1])){
      tempAnswer <- NULL
      firstVal <- TRUE
      for(col in c(1:dim(answerMatrix)[2])){
        if(answerMatrix[row,col] != 0){
          if(col == 1){
            tempAnswer <- paste0(answerMatrix[row,col])
            firstVal <- FALSE
          } else {
            if(firstVal){
              if(answerMatrix[row,col] == 1){
                tempAnswer <- paste0("x", freeVarCol[col-1])
                firstVal <- FALSE
              } else {
                tempAnswer <- paste0(answerMatrix[row,col], " * x", freeVarCol[col-1])
              }
            } else{
              if(abs(answerMatrix[row,col]) == 1){
                if(answerMatrix[row,col] > 0){
                  tempAnswer <- paste0(tempAnswer, " + ", "x", freeVarCol[col-1])
                } else {
                  tempAnswer <- paste0(tempAnswer, " - ", "x", freeVarCol[col-1])
                }
              } else {
                if(answerMatrix[row,col] > 0){
                  tempAnswer <- paste0(tempAnswer, " + ", answerMatrix[row,col], " * x", freeVarCol[col-1])
                } else {
                  tempAnswer <- paste0(tempAnswer, " - ", -answerMatrix[row,col], " * x", freeVarCol[col-1])
                }
              }
            }
          }
        }
      }
      answer <- c(answer, tempAnswer)
    }
    answer <- matrix(answer)
  }
  return(answer)
}

# Find Leading Variables and Free Variables
# Less than accuracy = 0 (default = NA)
# return Leading Variables as "L"
# return Free Variables as "F"
FindVariables <- function(REFMatrix, accuracy=NA){

  variables <- NULL
  varcol <- 0
  Lcount <- 0
  Fcount <- 0
  for(row in c(1:dim(REFMatrix)[1])){
    for(col in c(varcol+1:dim(REFMatrix)[2])){
      varcol <- col
      if(varcol > dim(REFMatrix)[2]){
        break
      }
      if(is.na(accuracy)){
        if(REFMatrix[row, col] != 0){
          variables[varcol] <- "L"
          Lcount <- Lcount + 1
          break # next row
        } else {
          variables[varcol] <- "F"
          Fcount <- Fcount + 1
          next # next col
        }
      } else {
        if(abs(REFMatrix[row, col]) > accuracy){
          variables[varcol] <- "L"
          Lcount <- Lcount + 1
          break # next row
        } else {
          variables[varcol] <- "F"
          Fcount <- Fcount + 1
          next # next col
        }
      }
      
    }
  }
  return(list(variables, row, c(Lcount, Fcount)))
}

# Modified Solve Strictly Triangular Form by Back Substitution
# (only works for square matrix)
# (will have the same result as solve(STFMatrix, RHSMatrix))
MSolveSTF <- function(STFMatrix, RHSMatrix){
  answerMatrix <- matrix(NA, dim(RHSMatrix)[1], dim(RHSMatrix)[2])
  for(row in c((dim(STFMatrix)[1]):1)){
    answerMatrix[row, ] <- RHSMatrix[row, ]
    for(col in c(row:dim(STFMatrix)[2])){
      if(col < dim(STFMatrix)[2]){
        answerMatrix[row, ] <- answerMatrix[row, ] - STFMatrix[row, col+1]*answerMatrix[col+1,]
      }
    }
  }
  return(answerMatrix)
}


# Find Next Pivot (Leading Variable)
# return next pivot location
# [[Deprecated]]
FindNextPivot <- function(REFMatrix, currentPivot){
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

#' TestConsistent
#' 
#' Test if a linear system is consistent or inconsistent
#' @param Matrix any matrix
#' @param attachMatrix matrix or vector (only support vector for current version)
#' @return consistent: TRUE; inconsisitent: FALSE
#' @export
#' 
# TODO: attachMatrix support matrix
isConsistent <- function(Matrix, attachMatrix){
  if(is.vector(attachMatrix)){
    attachMatrix <- matrix(attachMatrix)
  }
  if(dim(Matrix)[1] != dim(attachMatrix)[1]){
    return(FALSE)
  }
  
  temp <- GaussianElimination(Matrix, attachMatrix, FRAC = FALSE, SOLVE = FALSE)
  REFMatrix <- temp$RowEchelonForm
  attachVector <- temp$AttachVector
  
  for(row in c(dim(REFMatrix)[1]:1)){
    if(allZero(REFMatrix[row, ])){
      if(!allZero(attachVector[row, ])){
        return(FALSE)
      }
    } else {
      return(TRUE)
    }
  }
}

# Test if all element is zero
allZero <- function(Vector){
  for(i in c(1:length(Vector))){
    if(Vector[i] != 0){
      return(FALSE)
    }
  }
  return(TRUE)
}