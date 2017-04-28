function [valCT, valCQ, valCP, vecCTCONV, vecCQCONV, vecCPCONV, vecDISTHRUST, vecDISTORQUE, vecDISNORM] = fcnROTORFORCE(nind, nfree, nfreecs, thrustind, thrustfree, thrustCFfree, axialind, axialfree, inddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV,  vecCQCONV, vecCPCONV, vecCPRADI)
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
thrust = sum(thrustind)+sum(thrustfree)+sum(thrustCFfree);
torque = sum(axialind.*vecCPRADI+axialfree.*vecCPRADI+inddrag.*vecCPRADI);
power = torque*2.*pi.*(valRPM./60);

% Calculate non-dimensionalized coefficient using US customary definitions
% CT = T/(rho*Omega^2*D^4)
valCT = thrust/(((valRPM/60)^2)*((valDIA)^4));

% CQ = Q/(rho*A*Omega^2*D^3)
valCQ = torque/(((valRPM/60)^2)*((valDIA)^5));

% CP = P/(rho(RPM/60)^3*D^5);
valCP = power/((valRPM/60)^3*(valDIA^5));

% Convergence Thrust (average thrust across 1 full rotation)
temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
vecCTCONV(temp)= valCT;
vecCQCONV(temp) = valCQ;
vecCPCONV(temp) = valCP;

% Flow distributions
vecDISTHRUST = thrustind + thrustfree + thrustCFfree;
vecDISTORQUE = axialind + axialfree + inddrag;
vecDISNORM = nind + nfree + nfreecs;
end

