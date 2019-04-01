
Run all Table7.R, Table7Lambda.R, Table7MVNKnownLambda.R, Table7MVNUnknownLambda.R to get the results for Table 7.

Table7Function.R contains function to run simulation for sparse biclustering and MVN biclustering with lambda being specified.  Other methods such as the ssvd and plaid are also run.

BiclusterLambda.R contains function to run the simulation for sparse biclustering with lambda automatically chosen using BIC.

MVNLambdaKnown.R contains function to run the simulation for MVN biclustering with known Sigma and Delta, and lambda automatically chosen using BIC.

MVNLambdaKnown.R contains function to run the simulation for MVN biclustering with unknown Sigma and Delta, and lambda automatically chosen using BIC.


Table7.R gives the results for all methods except for sparse biclustering with lambda automatically chosen and MVN biclustering with lambda automatically chosen.

Table7BiclusterLambda.R gives results for sparse biclustering with lambda automatically chosen

Table7MVNKnownLambda.R gives results for MVN biclustering with known Sigma and Delta, and with lambda automatically chosen

Table7MVNUnknownLambda.R gives results for MVN biclustering with unknown Sigma and Delta, and lambda automatically chosen.  (This R script only provide result for one iteration.  The user needs to run it for 50 iterations and take the average of the provided results.)
