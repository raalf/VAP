# VAPTOR
The MATLAB version of the higher-order potential flow code FreeWake applied for ROTORS

Branch:
- Devin (VAPTOR)
    - Modified version of VAP 
    - Focuses on multirotor vehicle application

Current State:
- Reads related VAPTOR input
- Generates DVE using the VAP Generating DVE
- Calculate appropriate inflow velocity vector
- Timesteping
    - Move rotor
    - Generate new wake elements
	- Create and solve WD-Matrix for new elements
	- Solve wing D-Matrix with wake-induced velocities
    - Generate rotor resultant
	- Solve entire WD-Matrix
- Incorporate multiple blades

Next To-Do:
- Calculate forces
- Relax wake