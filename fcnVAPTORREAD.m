function [flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, ...
    valDELTAE, seqALPHAR, valJ, valDENSITY, valKINV, valAREA, valDIA,...
    vecROTAX, valPANELS, matGEOM, vecAIRFOIL, vecN, vecM, vecSYM, ...
    valINTERF] = fcnVAPTORREAD(strFILE)

% INPUT:
%   strFILE - file name of input text file in the local directory (or if not, with the appropriate path in the name)
% OUTPUT:
%   flagRELAX - 0 if fixed wake, 1 if relaxed
%   flagSTEADY - 0 if unsteady, 1 if steady

%   valMAXTIME - maximum number of timesteps
%   valMINTIME - minimum number of timestep
%   valAZNUM  - number of azimuth locations
%   valDELTAE - convergence criteria of change in span efficiency between timesteps

%   seqALPHAR - sequence of rotor plane angle of attacks
%   valJ - Advance ratio
%   valDENSITY - fluid density, kg/m^3
%   valKINV - kinematic viscosity (1.46e-05 as standard)

%   valAREA - projected wing area (m^2)
%   valDIA - propeller diameter (m)
%   vecROTAX - rotation axis [x y z] (m)

%   valPANELS - number of wing panels
%   matGEOM - 2 x 5 x valPANELS matrix, with (x,y,z) coords of edge points, and chord and twist at each edge
%   vecSYM - valPANELS x 1 vector of 0, 1, or 2 which denotes the panels with symmetry condition (1 or 2 being local edge number)
%   vecAIRFOIL - valPANELS x 1 vector of airfoil numbers for the panels

%   vecN - valPANELS x 1 vector of spanwise elements per DVE
%   vecM - valPANELS x 1 vector of chordwise elements per DVE

%   valINTERF - interference drag value (%)

fp = fopen(strFILE);

%% Reading header flags
% Reading relaxed wake flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
flagRELAX = fscanf(fp,'%d');

% Reading steady or unsteady flag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
flagSTEADY = fscanf(fp,'%d');

%% Reading time step information
% Reading maximum number of time steps
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valMAXTIME = fscanf(fp,'%d');

% Reading minimum number of time steps
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valMINTIME = fscanf(fp,'%d');

% Reading number of azmith locations
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valAZNUM = fscanf(fp,'%lf');

% Reading deltae
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valDELTAE = fscanf(fp,'%lf');

%% Reading flow conditions

% Reading sequence of rotor alphas to analyze
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
seqALPHAR = fscanf(fp,'%lf');

% Reading advanced ratio to be considered
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valJ = fscanf(fp,'%lf');

% Reading density
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valDENSITY = fscanf(fp,'%lf');

% Reading density
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valKINV = fscanf(fp,'%lf');
%% Reading Rotor Reference Values
% Reading reference aero
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valAREA = fscanf(fp,'%lf');

% Reading rotor diameter
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valDIA = fscanf(fp,'%lf');

ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
vecROTAX = fscanf(fp,'%lf');
vecROTAX = vecROTAX';

%% Reading panel/rotor/lifting line information
% Reading No. of panels
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valPANELS = fscanf(fp,'%lf');

%% Reading panel information and geometry

vecN = zeros(valPANELS,1);
vecM = zeros(valPANELS,1);
vecAIRFOIL = zeros(valPANELS,1);
vecSYM = zeros(valPANELS,1);

for i = 1:valPANELS
    % Reading spanwise 'n'
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecN(i) = fscanf(fp,'%lf',1);
    
    % Reading chordwise 'm'
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecM(i) = fscanf(fp,'%lf',1);
    
    % Reading airfoil number
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecAIRFOIL(i) = fscanf(fp,'%lf',1);
    
    % Reading symmetry information
    ch = fscanf(fp,'%c',1);
    while(ch~='=');
        ch = fscanf(fp,'%c',1);
    end
    vecSYM(i) = fscanf(fp,'%lf',1);
    
    % Skipping geometry column header
    fgets(fp);
    fgets(fp);
    
    % Reading geometry
    % Explanation below:
    %{
        info_geometry(x,y,z)
            x is for the left or right point
                1 left
                2 right
            y is for the values
                1 x
                2 y
                3 z
                4 chord
                5 epsilon
                6 boundary condition
            z is panel number
    %}
    
    matGEOM(1,:,i) = fscanf(fp,'%lf',5);
    matGEOM(2,:,i) = fscanf(fp,'%lf');
    
end

% Reading intereference drag
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end

valINTERF = fscanf(fp,'%lf',1);

fclose(fp);

clear ans ch i j fp idx1





















