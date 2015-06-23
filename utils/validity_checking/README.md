Semi-Automated checking
=======================

This should give a slightly simpler approach at running all the possible
communication frameworks after a new feature was added to see if anything
broke.  You can either do everything automatically (this will run all the
frameworks) or modify the generated files a little.

### Step 1
Compile HASEonGPU (ideally on one of the compute nodes)
 - This should automatically add the new binary file to the MATLAB example.

### Step 2
On the HeadNode (hypnos5), run the script
`utils/validity_checking/integration_testing.sh` (should be located in the
`HASEonGPU` folder)
 - This script needs as a parameter the folder where the comparison should be
   located. This location must be empty (empty folder is OK, nothing is also
   OK)
 - The script will copy the experiment 3 times (graybat,mpi,threaded) to the
   comparison folder
 - The laserPumpCladding example will be modified to use the correct
   communication framework

### Step 3
Run the experiments
 - Switch to the comparison folder, it has sub-folders for all possible
   communication frameworks
 - edit the file `integration_testing.m` so that it runs (and checks) only the
   experiments you want
 - if you are satisfied, run `qsub submit.sh`

### Step 4
Wait until all 3 experiments are done
 - The MATLAB script `integration_testing.m` will automatically run all defined
   experiments in sequence (this approach needs only 1 MATLAB licencse)
 - Afterwards, the script `gain_analysis.m` is called inside each of the
   experiment folders


### Step 5
Look at the output
 - `gain_analysis.m` will create the files `result.txt` and `gain_curve_*.txt`
   that contain information.


### Troubleshooting
If the program crashes (usually because the experiment itself crashed), there
is no really nice way to continue. At least in the context of this framework.

What you CAN do: if an experiment already completed, you can start the
`gain_analysis.m` manually. To do that, go into the experiment folder
(graybat,mpi, or threaded) and execute 
`matlab -r "addpath .. , gain_analysis ,exit"` 
This will also place files in the parent folder.
