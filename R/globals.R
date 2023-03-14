# Global variables - required to solve CRAN check NOTE
# Author: Hugo Pedder
# Date created: 2018-09-10

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("studyID", "time", "treatment",
                                                       "fup", "fups", "study", "arm", "fupcount", "y",
                                                       "2.5%", "50%", "97.5%", "ref.median",
                                                       "Var1", "value", "Parameter", ".",
                                                       "int.fun", ".data",
                                                       "bins", "jagsdata", "n", "narm",
                                                       "tbin", "tmax", "tmin"))
