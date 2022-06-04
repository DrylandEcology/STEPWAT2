| Unix | Windows | Release | License | Coverage | Downloads |
| :---- | :---- | :---- | :---- | :---- | :---- |
[ ![Travis build status][1]][2] | [![Appveyor build status][3]][4] | [ ![github release][5]][6] | [![license][7]][8] | [![codecov status][9]][10] | [![github downloads][11]][12] |

[1]: https://travis-ci.org/DrylandEcology/STEPWAT2.svg?branch=master
[2]: https://travis-ci.org/DrylandEcology/STEPWAT2
[3]: https://ci.appveyor.com/api/projects/status/qx2o3j4jpp0ej0on/branch/master?svg=true
[4]: https://ci.appveyor.com/project/DrylandEcologyGit/stepwat2/branch/master
[5]: https://img.shields.io/github/release/DrylandEcology/STEPWAT2.svg?label=current+release
[6]: https://github.com/DrylandEcology/STEPWAT2/releases
[7]: https://img.shields.io/github/license/DrylandEcology/STEPWAT2.svg
[8]: https://www.gnu.org/licenses/gpl.html
[9]: https://codecov.io/gh/DrylandEcology/STEPWAT2/branch/master/graph/badge.svg
[10]: https://codecov.io/gh/DrylandEcology/STEPWAT2
[11]: https://img.shields.io/github/downloads/DrylandEcology/STEPWAT2/total.svg
[12]: https://github.com/DrylandEcology/STEPWAT2

<br>

# STEPWAT2

<br>

## Basic instructions for working with the code
1) Clone (obtain code from online github repository):
`git clone --recursive https://github.com/DrylandEcology/STEPWAT2.git`

2) Change working directory to STEPWAT2 folder: `cd STEPWAT2`

3) Compile the code (make an executable application): `make`

Other tasks:
* Clean local repository, e.g., after having run some tests: `git clean -d -f`

* Reset files to their online github repository version, i.e., overwrite changes made to
  local files: `git reset --hard`

* Pull down/sync latest commits/version from github repository: `git pull`


#### Code development and pull-requests
* Please follow instructions in our
[guidelines](https://github.com/Burke-Lauenroth-Lab/workflow_guidelines)


<br>

## Basic instruction for running a test project

* Print a brief explanation of options:
```
./stepwat --help
```

```
>   Usage : steppe [-d startdir] [-f files.in] [-q] [-e] [-o] [-g]
>      -d : supply working directory (default=.)
>      -f : supply list of input files (default=files.in)
>      -q : quiet mode, don't print message to check logfile.
>      -p : prints progress bar
>      -e : echo initialization results to logfile
>      -o : write SOILWAT output to output files. Contains average over all iterations and standard deviation.
>      -g : use gridded mode
>      -i : write SOILWAT output to output files for each iteration
>-STdebug : generate sqlite database with STEPWAT information
```

* Run the gridded version (-g) of STEPWAT2:

```
cd testing.sagebrush.master/
./stepwat -f files.in -g
```


* Run the non-gridded version of STEPWAT2 from the Stepwat_Inputs/ folder using SOILWAT2 to drive the water cycle:

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in
```


* Run non-gridded version of STEPWAT using SOILWAT and output all variables passed between Stepwat and SOILWAT:

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in -ssxwdebug.in
```


* Run STEPWAT2:

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in
```

* Run STEPPE (get SOILWAT output for all iterations and average)

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in -o -i
```

* Alternatively, use `makefile` targets (compile, copy, and execute in one command)
  - Non-gridded version, with SOILWAT2 water cycle, and with iteration and aggregated SOILWAT2-output:
    ```
    make bint_testing_nongridded
    ```
  - Gridded version with SOILWAT2 water cycle
    ```
    make bint_testing_gridded
    ```
  - If you need to clean up first, then run, e.g.,
    ```
    make cleanall bint_testing_nongridded
    ```



<br>

## When switching development branches:

* Checkout the new branch:
```
git checkout -b [branch name]
```
* To ensure all submodules are updated correctly:
```
git submodule update --init --recursive
```


<br>

## Compare output of a new branch to output from a previous (reference) release

Depending on the purpose of the development branch
the new output should be exactly the same as reference output or
differ in specific ways in specific variables.

The following steps provide a starting point for such comparisons:

```
# Simulate on master branch and copy output to "Output_ref"
git checkout master
make clean bint_testing_nongridded
cp -r testing.sagebrush.master/Stepwat_Inputs/Output testing.sagebrush.master/Stepwat_Inputs/Output_ref

# Switch to development branch <branch_xxx> and run the same simulation
git checkout <branch_xxx>
make clean bint_testing_nongridded

# Compare the two sets of outputs
#   * Lists all output files and determine if they are exactly they same
diff -qsr testing.sagebrush.master/Stepwat_Inputs/Output testing.sagebrush.master/Stepwat_Inputs/Output_ref
#   * Produce two figures (scatter plots and time series of biomass)
#     that are saved as PDFs at testing.sagebrush.master/
Rscript tools/compare_bmassavg_against_ref.R
```



<br>

## Doxygen documentation:
STEPWAT2 uses Doxygen to automatically generate documentation for the code base. For information on how to generate the documentation see our [GitHub wiki](https://github.com/DrylandEcology/STEPWAT2/wiki/For-Developers).

<br>

## Note: repository renamed from StepWat to STEPWAT2 on Feb 23, 2017

All existing information should [automatically be redirected](https://help.github.com/articles/renaming-a-repository/) to the new name.

Contributors are encouraged, however, to update local clones to [point to the new URL](https://help.github.com/articles/changing-a-remote-s-url/), i.e.,
```
git remote set-url origin https://github.com/Burke-Lauenroth-Lab/STEPWAT2.git
```
