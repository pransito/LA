## LA options and prep
# currently centered-only
use.z    <- 1 
# use absolute value of loss as predictor
loss.abs <- 1
# use only some trials; i.e. use only trials where subjects where overall uncertain
only.some <- 0 
# ranges to get only some trials
uncertainty.up <- 0.6
uncertainty.lo <- 0.4
# calculate ist?
do_lmlist <- 0

# if there are more than x% percent missings in response then drop-out
missing.cutoff <- 0.07 

# prep the final matrix; will be a multiple column matrix
subject <- c()
p       <- c()
r       <- c()