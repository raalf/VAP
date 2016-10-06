# VAP
The MATLAB version of the higher-order potential flow code FreeWake.


KNOWN ISSUES WITH VAP1.0:
-If relaxing the wake, the very first row of elements emitted (timestep 0) are not connected to the second row. 
-Splits in the wings are not handled correctly.

TO-DO:

- Remove structures
- Vectorize functions

We can implement some of the low-level functions such as vortex sheet induction from HVAP which are written well. 

Areas to focus on:

- DVE creation
    - generateDVEs_v2
    - generateSingfct
    - panelBoundary
    - generateWakeDVEs

- D-Matrix
    - fcnAssembleWingD
    - fcnDVEKinCond
    - fcnNewWakeVorticity
    - fcnSolveDWakeMatrix

