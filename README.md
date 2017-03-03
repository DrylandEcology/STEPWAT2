# StepWat
folder

#To run the gridded version of STEPWAT from the testing.sagebrush.MT_drs folder:
This call will save all SOILWAT2 output
 ./stepwat    -f  files.in  -g -s  -o ../sw_src/testing/files_step_soilwat_grid.in

This call will not save any SOILWAT2 output
./stepwat -f files.in -g -s
 
#To run the non-gridded version of STEPWAT from the Stepwat_Inputs folder:
This call will save all SOILWAT2 output
 ./stepwat    -f  files.in -s -o ../../sw_src/testing/files_step_soilwat.in

This call will not save any SOILWAT2 output
./stepwat -f files.in -s


## Note: repository renamed from StepWat to STEPWAT2 on Feb 23, 2017

All existing information should [automatically be redirected](https://help.github.com/articles/renaming-a-repository/) to the new name.

Contributors are encouraged, however, to update local clones to [point to the new URL](https://help.github.com/articles/changing-a-remote-s-url/), i.e., 
```
git remote set-url origin https://github.com/Burke-Lauenroth-Lab/STEPWAT2.git
```

If using GitHub Desktop you may want to remove existing clones, update the local clones via command line, re-start GitHub Desktop and add
a new clone of the updated repository.
