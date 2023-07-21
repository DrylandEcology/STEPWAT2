#!/bin/bash


# STEPWAT2 should meet the following expectations as formulated with
# PR #528 (https://github.com/DrylandEcology/STEPWAT2/pull/528) that was
# merged into the main branch on Sep 18, 2022 with commit
# (https://github.com/DrylandEcology/STEPWAT2/commit/7dce3c7a21a7ff3eee563765319301fa952bbaa6)
# "Exactly reproduce random number sequences"

echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"# This script checks that STEPWAT2 meets the following expected/desired behavior:"$'\n'\
$'\n'\
"#  * if seed != 0 (output is reproduced among runs)"$'\n'\
"#      ** weather is identical among runs and cells"$'\n'\
"#      ** weather is different among years, iterations and seeds"$'\n'\
$'\n'\
"# * if seed == 0 (output cannot be reproduced among runs)"$'\n'\
"#      ** weather is different among cells, years, iterations and runs"$'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"



#--- Setup ---------
# Filename tag of gridded biomass per cell k: '${tag_gridded_biomass_cellk}k.csv'
#  * branch 'Seed_Dispersal': "bmass"
#  * branch 'master': "g_bmassavg"
tag_gridded_biomass_cellk="bmass"

# Filename tag of gridded biomass averaged across cells: '${fname_gridded_biomass_meancell}'
#  * branch 'Seed_Dispersal': "all_cell_average_biomass.csv"
#  * branch 'master': "g_bmass_cell_avg.csv"
fname_gridded_biomass_meancell="all_cell_average_biomass.csv"


# Paths to gridded and nongridded example runs
dirname_gridded="testing.sagebrush.master"
dirname_nongridded="testing.sagebrush.master/Stepwat_Inputs"


#--- Function definitions ------
#--- Functionality to update model.in
fname_modelin="${dirname_nongridded}""/Input/model.in"
fname_modelin_backup="${fname_modelin}""_backup"

# create backup copy of modelin
cp "${fname_modelin}" "${fname_modelin_backup}"

# Function to set niter nyrs  seed
# $1 niter
# $2 nyrs
# $3 seed
set_modelinput () {
  printf "%d %d %d\n" "$1" "$2" "$3" > "${fname_modelin}"
}

reset_modelinput() {
  cp "${fname_modelin_backup}" "${fname_modelin}"
}

# Re-instate model.in from backup when script exists (failure or successfully)
trap reset_modelinput EXIT


#--- Variables that count overall failures and successes
count_failures=0 # initialize variable to count failures
count_successes=0 # initialize variable to count successes



#--- Function to check for errors produced by make targets (but avoid false hits)
# $1 Text string
# Returns lines that contain "failed", "abort", "trap", " not ", or "error"
# ("error" besides "Werror", "Wno-error")
check_target_error() {
  res_error=$(echo "$1" | awk 'tolower($0) ~ /[^-w]error|failed|abort|trap|( not )/')

  if [ "${res_error}" ]; then
    echo "Target: warnings/errors:"
    echo "${res_error}"
  else
    echo "Target: success."
  fi
}


#--- Function to check for successful/failed expectations
# $1 text string produced by Rscript: TRUE signals success; FALSE signals failure
# On failure, increment count_failures; on success, increment count_successes
# Print $1
check_Rscript_output() {
  res_failure=$(echo "$1" | awk 'tolower($0) ~ /false/')

  if [ "${res_failure}" ]; then
    count_failures=$((count_failures + 1))
  else
    count_successes=$((count_successes + 1))
  fi

  echo "$1"
}

# $1 text string produced by diff: empty string signals success; content signals failure
# On failure, increment count_failures; on success, increment count_successes
# Print $1
check_diff_output() {
  if [ "${res}" ]; then
    count_failures=$((count_failures + 1))
    echo "${res}"
  else
    count_successes=$((count_successes + 1))
    echo "TRUE"
  fi
}




#------------ Checks -----------------------------------------------------------
echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- Compile STEPWAT2 ------"
res=$(make clean stepwat 2>&1)
check_target_error "${res}"


echo $'\n'$'\n'$'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- Test nongridded mode ----------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"

echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- * if seed != 0 (output is reproduced among runs) ------"$'\n'\
"# Set model.in: 2 100 7 (niter, nyrs, seed)"

set_modelinput 2 100 7
dirtest="Output_i2y100s7_nongr"

rm -r "${dirname_nongridded}"/"${dirtest}"
res=$(make bint_testing_nongridded 2>&1)
check_target_error "${res}"

cp -r "${dirname_nongridded}"/Output "${dirname_nongridded}"/"${dirtest}"

res=$(make bint_testing_nongridded 2>&1)
check_target_error "${res}"


echo $'\n'\
"#--- ** weather (and all STEPWAT2 output) is exactly identical among runs ------"$'\n'

res=$(diff -aqr "${dirname_nongridded}"/Output "${dirname_nongridded}"/"${dirtest}" 2>&1)
check_diff_output "${res}"


echo $'\n'\
"#--- ** weather is different among years ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'x <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_nongridded}"/Output \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among iterations ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'x <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("StdDev", "StdDev.1")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_nongridded}"/Output \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among seeds ------"$'\n'\
"# Set model.in: 2 100 6 (niter, nyrs, seed) # or any seed other than 0 and 7"

set_modelinput 2 100 6

res=$(make bint_testing_nongridded 2>&1)
check_target_error "${res}"

# Expect differences
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'dir_out2 <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x1 <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e 'x2 <- read.csv(file.path(dir_out2, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e '!isTRUE(all.equal(x1, x2))' \
  "${dirname_nongridded}"/Output \
  "${dirname_nongridded}"/"${dirtest}" \
  2>&1)
check_Rscript_output "${res}"





echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- * if seed == 0 (output cannot be reproduced among runs) ------"$'\n'\
"# Set model.in: 2 100 0 (niter, nyrs, seed)"

set_modelinput 2 100 0
dirtest="Output_i2y100s0_nongr"

rm -r "${dirname_nongridded}"/"${dirtest}"
res=$(make bint_testing_nongridded 2>&1)
check_target_error "${res}"

cp -r "${dirname_nongridded}"/Output "${dirname_nongridded}"/"${dirtest}"

res=$(make bint_testing_nongridded 2>&1)
check_target_error "${res}"


echo $'\n'\
"#--- ** weather is different among years ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'x <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_nongridded}"/Output \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among iterations ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'x <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("StdDev", "StdDev.1")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_nongridded}"/Output \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among runs ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'dir_out2 <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x1 <- read.csv(file.path(dir_out, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e 'x2 <- read.csv(file.path(dir_out2, "bmassavg.csv"))[, c("PPT", "Temp")]' \
  -e '!isTRUE(all.equal(x1, x2))' \
  "${dirname_nongridded}"/Output \
  "${dirname_nongridded}"/"${dirtest}" \
  2>&1)
check_Rscript_output "${res}"




echo $'\n'$'\n'$'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- Test gridded mode -------------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"

echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- * if seed != 0 (output is reproduced among runs) ------"$'\n'\
"# Set model.in: 2 100 7 (niter, nyrs, seed)"

set_modelinput 2 100 7
dirtest="Output_i2y100s7_gr"

rm -r "${dirname_gridded}"/"${dirtest}"
res=$(make bint_testing_gridded 2>&1)
check_target_error "${res}"

cp -r "${dirname_gridded}"/Output "${dirname_gridded}"/"${dirtest}"

res=$(make bint_testing_gridded 2>&1)
check_target_error "${res}"


echo $'\n'\
"#--- ** weather (and all STEPWAT2 output) is exactly identical among runs ------"

res=$(diff -aqr "${dirname_gridded}"/Output "${dirname_gridded}"/"${dirtest}" 2>&1)
check_diff_output "${res}"


echo $'\n'\
"#--- ** weather (and all STEPWAT2 output) is identical among cells ------"

res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'tag_gridded_biomass_cellk <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- lapply(seq_len(4), function(k) read.csv(file.path(dir_out, paste0(tag_gridded_biomass_cellk, k - 1, ".csv")))[, c("PPT", "Temp")])' \
  -e 'sapply(seq_len(4), function(k) all.equal(x[[k]], x[[1]]))' \
  "${dirname_gridded}"/Output \
  "${tag_gridded_biomass_cellk}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among years ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'fname_gridded_biomass_meancell <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- read.csv(file.path(dir_out, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_gridded}"/Output \
  "${fname_gridded_biomass_meancell}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among iterations ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'tag_gridded_biomass_cellk <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- lapply(seq_len(4), function(k) read.csv(file.path(dir_out, paste0(tag_gridded_biomass_cellk, k - 1, ".csv")))[, c("StdDev", "StdDev.1")])' \
  -e 'sapply(x, function(xk) apply(xk, 2, sd) > 0)' \
  "${dirname_gridded}"/Output \
  "${tag_gridded_biomass_cellk}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among seeds ------"$'\n'\
"# Set model.in: 2 100 6 (niter, nyrs, seed) # or any seed other than 0 and 7"

set_modelinput 2 100 6

res=$(make bint_testing_gridded 2>&1)
check_target_error "${res}"

# Expect differences
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'dir_out2 <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'fname_gridded_biomass_meancell <- as.character(commandArgs(TRUE)[[3L]])' \
  -e 'x1 <- read.csv(file.path(dir_out, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e 'x2 <- read.csv(file.path(dir_out2, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e '!isTRUE(all.equal(x1, x2))' \
  "${dirname_gridded}"/Output \
  "${dirname_gridded}"/"${dirtest}" \
  "${fname_gridded_biomass_meancell}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- * if seed == 0 (output cannot be reproduced among runs) ------"$'\n'\
"# Set model.in: 2 100 0 (niter, nyrs, seed)"

set_modelinput 2 100 0
dirtest="Output_i2y100s0_gr"

rm -r "${dirname_gridded}"/"${dirtest}"
res=$(make bint_testing_gridded 2>&1)
check_target_error "${res}"

cp -r "${dirname_gridded}"/Output "${dirname_gridded}"/"${dirtest}"

res=$(make bint_testing_gridded 2>&1)
check_target_error "${res}"


echo $'\n'\
"#--- ** weather is different among cells ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'tag_gridded_biomass_cellk <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- lapply(seq_len(4), function(k) read.csv(file.path(dir_out, paste0(tag_gridded_biomass_cellk, k - 1, ".csv")))[, c("PPT", "Temp")])' \
  -e 'sapply(seq_len(4)[-1], function(k) !isTRUE(all.equal(x[[k]], x[[1]])))' \
  "${dirname_gridded}"/Output \
  "${tag_gridded_biomass_cellk}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among years ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'fname_gridded_biomass_meancell <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- read.csv(file.path(dir_out, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e 'apply(x, 2, sd) > 0' \
  "${dirname_gridded}"/Output \
  "${fname_gridded_biomass_meancell}" \
  2>&1)
check_Rscript_output "${res}"


echo $'\n'\
"#--- ** weather is different among iterations ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'tag_gridded_biomass_cellk <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'x <- lapply(seq_len(4), function(k) read.csv(file.path(dir_out, paste0(tag_gridded_biomass_cellk, k - 1, ".csv")))[, c("StdDev", "StdDev.1")])' \
  -e 'sapply(x, function(xk) apply(xk, 2, sd) > 0)' \
  "${dirname_gridded}"/Output \
  "${tag_gridded_biomass_cellk}" \
  2>&1)
check_Rscript_output "${res}"



echo $'\n'\
"#--- ** weather is different among runs ------"
res=$(Rscript \
  -e 'dir_out <- as.character(commandArgs(TRUE)[[1L]])' \
  -e 'dir_out2 <- as.character(commandArgs(TRUE)[[2L]])' \
  -e 'fname_gridded_biomass_meancell <- as.character(commandArgs(TRUE)[[3L]])' \
  -e 'x1 <- read.csv(file.path(dir_out, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e 'x2 <- read.csv(file.path(dir_out2, fname_gridded_biomass_meancell))[, c("PPT", "Temp")]' \
  -e '!isTRUE(all.equal(x1, x2))' \
  "${dirname_gridded}"/Output \
  "${dirname_gridded}"/"${dirtest}" \
  "${fname_gridded_biomass_meancell}" \
  2>&1)
check_Rscript_output "${res}"



echo $'\n'$'\n'$'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#--- Overall summary ---------------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"$'\n'

if [ "${count_failures}" -gt 0 ]; then
  echo "Failure: ""${count_failures}"" expectations failed and ""${count_successes}"" expectations confirmed."
else
  echo "Success: ""${count_successes}"" expectations confirmed."
fi

echo $'\n'\
"#-----------------------------------------------------------------------"$'\n'\
"#-----------------------------------------------------------------------"$'\n'


#--- Cleanup ---
reset_modelinput
