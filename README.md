# STEPWAT2

<br>

## Basic instructions for working with the code
1) Clone (obtain code from online github repository):
`git clone --single-branch --recursive https://github.com/Burke-Lauenroth-Lab/STEPWAT2.git`

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
> Usage: stepwat [-d startdir] [-f files.in] [-q] [-s] [-e] [-o] [-g]
>   -d : supply working directory (default=.)
>   -f : supply list of input files (default=files.in)
>   -q : quiet mode, don't print message to check logfile.
>   -s : use SOILWAT model for resource partitioning.
>   -e : echo initialization results to logfile
>   -o : print all the soilwat output in addition to the stepwat output
>   -g : use gridded mode
```

* Run the gridded version (-g) of STEPWAT2 using SOILWAT2 to drive the water cycle (-s):

```
cd testing.sagebrush.master/
./stepwat -f files.in -g -s
```


* Run the non-gridded version of STEPWAT2 from the Stepwat_Inputs/ folder using SOILWAT2
  to drive the water cycle (-s):

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in -s
```


* Run STEPPE (no SOILWAT2)

```
cd testing.sagebrush.master/Stepwat_Inputs/
./stepwat -f files.in
```


<br>

## Note: repository renamed from StepWat to STEPWAT2 on Feb 23, 2017

All existing information should [automatically be redirected](https://help.github.com/articles/renaming-a-repository/) to the new name.

Contributors are encouraged, however, to update local clones to [point to the new URL](https://help.github.com/articles/changing-a-remote-s-url/), i.e.,
```
git remote set-url origin https://github.com/Burke-Lauenroth-Lab/STEPWAT2.git
```

