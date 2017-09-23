## testing script for MatrixReduceProcess.R

source("R/MatrixReduceProcess.R")

# loading testing data
load(file = "TestingPurpose/testingMatrix.rda")
matrix1 <- testingMatrix[[1]]
matrix2 <- testingMatrix[[2]]

attvec1 <- testingMatrix[[3]]
attvec2 <- testingMatrix[[4]]

sunnym1 <- testingMatrix[[5]]
sunnym2 <- testingMatrix[[6]]
sunnyv1 <- testingMatrix[[7]]
sunnyv2 <- testingMatrix[[8]]

sunnym3 <- testingMatrix[[9]]
sunnyv3 <- testingMatrix[[10]]


## testing functions ##
# test for ReduceAugmentedMatrix
#print(matrix1)
#print(attvec1)
#print(ReduceAugmentedMatrix(matrix1, attvec1))

#print(matrix2)
#print(attvec2)
#print(ReduceAugmentedMatrix(matrix2, attvec2))

#print(sunnym3)
#print(sunnyv3)
#print(ReduceAugmentedMatrix(sunnym3, sunnyv3, PRINT=TRUE))

# test for GaussianElimination

#print(matrix1)
#print(attvec1)
#print(GaussianElimination(matrix1, attvec1))

#print(matrix2)
#print(attvec2)
#print(GaussianElimination(matrix2, attvec2))

#print(sunnym1)
#print(sunnyv1)
#print(GaussianElimination(sunnym1, sunnyv1, PRINT=TRUE))

#print(sunnym2)
#print(sunnyv2)
#print(GaussianElimination(sunnym2, sunnyv2, PRINT=TRUE))

print(sunnym3)
print(sunnyv3)
print(GaussianElimination(sunnym3, sunnyv3, PRINT=TRUE))