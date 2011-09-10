To run the code, first load from afs the same version of ROOT that I used:
source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.30.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh

Then produce a .root file with the model of the backgrounds and signal:
root.exe
.L create_tight.C+
new_RA4()

Then you run the CLs hypothesis testing code (with one-sided Profile Likelihood test statistic) on the resulting root file:
root.exe
.L StandardHypoTestInvDemo.C+
StandardHypoTestInvDemo("tight.root", "wspace", "SbModel","BModel", "data", 0, 3, true, 12, 6, 17, 5000)


