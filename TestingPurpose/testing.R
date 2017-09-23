## testing script for MatrixReduceProcess.R

source("R/MatrixReduceProcess.R")

# loading testing data
load(file = "TestingPurpose/testingMatrix.rda")
matrix1 <- testingMatrix[[1]]
matrix2 <- testingMatrix[[2]]

attvec1 <- testingMatrix[[3]]
attvec2 <- testingMatrix[[4]]

## testing functions ##
# test for ReduceAugmentedMatrix
print(matrix1)
print(attvec1)
print(ReduceAugmentedMatrix(matrix1, attvec1))

print(matrix2)
print(attvec2)
print(ReduceAugmentedMatrix(matrix2, attvec2))

# test for GaussianElimination

print(matrix1)
print(attvec1)
print(GaussianElimination(matrix1, attvec1))

print(matrix2)
print(attvec2)
print(GaussianElimination(matrix2, attvec2))
