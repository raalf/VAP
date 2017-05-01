function [valCT, valCQ, valCP, vecCTCONV, vecCQCONV, vecCPCONV, vecDISTHRUST, vecDISNORM] = fcnROTORFORCE(nind, nfree, nfreecs, thrustind, thrustfree, thrustCFfree,  thrustinddrag, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV,vecCQCONV, vecCPCONV, vecQARM, vecDVETE, vecTHETA)
%   This function uses the previously calculated induced and freestream
%   forces and calculates non-dimensionalized values.
%
%   Inputs:
% 
%   Outputs:
%   valCT - Thrust Coefficient
%   valCQ - Torque Coefficient
%   vecCTCONV - Convergence thrust, vector the length of valAZNUM
tempDi = zeros(size(axialfree,1),1);
tempTi = zeros(size(axialfree,1),1);
tempSi = zeros(size(axialfree,1),1);
tempDi(vecDVETE==3) = inddrag;
tempTi(vecDVETE==3) = thrustinddrag;
tempSi(vecDVETE==3) = sideinddrag;

% Flow distributions
vecDISNORM = nind + nfree + nfreecs;
vecDISTHRUST = thrustind + thrustfree + thrustCFfree + tempTi;
vecDISTAXIAL = axialind + axialfree + tempDi + axialCFfree; 
vecDISTSIDE = sideind + sidefree + sideCFfree + tempTi;

% Calculate total force and moment values per density
thrust = sum(thrustind)+sum(thrustfree)+sum(thrustCFfree) + sum(thrustinddrag);
Py = sum(vecDISTSIDE.*sin(vecTHETA) + vecDISTAXIAL.*sin(pi - vecTHETA));
Px = sum(vecDISTSIDE.*cos(vecTHETA) + vecDISTAXIAL.*cos(pi - vecTHETA));

torque = sum(axialind.*vecQARM)+sum(axialfree.*vecQARM)+sum(inddrag.*vecQARM(vecDVETE==3)) + sum(axialCFfree.*vecQARM);
power = torque*2.*pi.*(valRPM./60);
Mx = vecDISTHRUST.*(vecQARM.*sin(vecTHETA));
My = vecDISTHRUST.*(vecQARM.*cos(vecTHETA));

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

end

