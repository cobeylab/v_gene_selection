est_error_rate <- 0.0018
typical_length <- 300

# Prob both sequences identical
p00 <- dbinom(0,typical_length,est_error_rate) * dbinom(0,typical_length,est_error_rate)

# Prob one sequence has a mutation, the other none
p10 <- 2* dbinom(0,typical_length,est_error_rate) * dbinom(1, typical_length, est_error_rate)

# Prob each seq has a mutation (assuming it's never an identical mutation, i.e. prob that seqs differ by two mutations)
p11 <-  dbinom(1,typical_length,est_error_rate)^2

# Prob a sequence has two mutations, the other none
p20 <- 2 * dbinom(0,typical_length,est_error_rate) * dbinom(2, typical_length, est_error_rate)

p00 + p10 + p11 + p20




