library(PADIGM)
context("tests user functions")

allSPRigthOutput = data.frame("C00042" = as.vector(c("hsa:5105"=3,"hsa:47"=Inf)),
                   "C00036"= as.vector(c("hsa:5105"=Inf, "hsa:47"=Inf)));
testWrongInputDF = data.frame(c("a"));
testWrigthInputDF = data.frame("gene" = as.vector(c("hsa:5105","hsa:47")),
                               "metabolites"= as.vector(c("C00042","C00036")));

outputWrong <- getAllShortestDistances("hsa01100",testWrongInputDF);
outputRight <- getAllShortestDistances("hsa01100", testWrigthInputDF);
test_that("test-getAllShortestDistances", {
    expect_identical(outputWrong, allSPRigthOutput)
    expect_identical(outputRight, allSPRigthOutput)


})




