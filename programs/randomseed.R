############################## Set seeds #######################################
library(digest)
library(dplyr)
library(purrr)

# =====================================================
# (1) Convert input string into MD5 hash → integer → seed
# =====================================================
md5_to_seed <- function(input_string) {
  # Get MD5 as raw bytes (length 16)
  hash_raw <- digest(input_string, algo = "md5", serialize = FALSE, raw = TRUE)
  
  # Take first 4 bytes and combine into an integer
  raw_int <- sum(as.integer(hash_raw[1:4]) * 256^(0:3))
  
  # Map into R's valid seed range (0 … 2^31 - 1)
  seed_num <- raw_int %% (2^31 - 1)
  
  return(seed_num)
}

input_string <- "Wenxiao"
number <- md5_to_seed(input_string)
cat("MD5-based safe seed for", input_string, "is:", number, "\n")

# Set the seed value
set.seed(number)
cat("Initial seed number: ", number, "\n", sep = "")

# Generate a list of 10 random numbers
a <- 0
b <- 2^31 - 1
random_numbers <- sample(a:b, 10, replace = TRUE)

# Print the list
cat("Seed", paste(random_numbers, collapse = " "), "\n")
###############################################################################