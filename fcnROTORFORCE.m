function [valCT, valCQ, vecCTCONV] = fcnROTORFORCE(thrustind, thrustfree, axialind, axialfree, inddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV)
%   This function uses the previously calculated induced and freestream
%   forces and calculates non-dimensionalized values.
%
%   Inputs:
% 
%   Outputs:
%   valCT - Thrust Coefficient
%   valCQ - Torque Coefficient
%   vecCTCONV - Convergence thrust, vector the length of valAZNUM

% Calculate total force values per density
thrust = sum(thrustind)+sum(thrustfree);
torque = sum(axialind)+sum(axialfree)+sum(inddrag);

% Calculate non-dimensionalized coefficient using US customary definitions
% CT = T/(rho*Omega^2*D^4)
valCT = thrust/(((valRPM/60)^2)*((valDIA)^4));

% CQ = Q/(rho*A*Omega^2*R^3)
valCQ = torque/(((valRPM/60)^2)*((valDIA)^5));

% Convergence Thrust (average thrust across 1 full rotation)
temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
vecCTCONV(temp)= valCT;
valCTCONV = mean(vecCTCONV);

end

