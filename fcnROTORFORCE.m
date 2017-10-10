function [ thrustind, thrustfree, thrustCFfree, tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE] = fcnROTORFORCE(difthrustP, diffsideP, diffaxialP, nind, nfree, nfreecs, Pthrust, thrustind, thrustfree, thrustCFfree,  thrustinddrag, Ptorque, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCQCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecQARM, vecDVETE, vecTHETA, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE)
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
vecDISTHRUST = thrustind + thrustfree + thrustCFfree + tempTi + Pthrust + difthrustP;
vecDISAXIAL = axialind + axialfree + tempDi + axialCFfree; 
vecDISSIDE = sideind + sidefree + sideCFfree + tempSi;

% Calculate total force and moment values per density
thrust = sum(thrustind)+sum(thrustfree)+sum(thrustCFfree) + sum(thrustinddrag) + sum(Pthrust) + sum(difthrustP);
Fy = sum(vecDISSIDE.*sin(vecTHETA) + vecDISAXIAL.*sin(pi - vecTHETA));
Fx = sum(vecDISSIDE.*cos(vecTHETA) + vecDISAXIAL.*cos(pi - vecTHETA)); % ADD viscous

torque = sum(axialind.*vecQARM)+sum(axialfree.*vecQARM)+sum(inddrag.*vecQARM(vecDVETE==3)) + sum(axialCFfree.*vecQARM) + sum(Ptorque) + sum(diffaxialP.*vecQARM); %sum(diffaxialP);
power = torque*2.*pi.*(valRPM./60);
Mx = sum(vecDISTHRUST.*(vecQARM.*sin(vecTHETA)));
My = sum(vecDISTHRUST.*(vecQARM.*cos(vecTHETA)));

% Calculate non-dimensionalized coefficient using US customary definitions
% CT = T/(rho*Omega^2*D^4)
valCT = thrust/(((valRPM/60)^2)*((valDIA)^4));
valFy = Fy/(((valRPM/60)^2)*((valDIA)^4));
valFx = Fx/(((valRPM/60)^2)*((valDIA)^4));

valCQ = torque/(((valRPM/60)^2)*((valDIA)^5));
valCP = power/((valRPM/60)^3*(valDIA^5));
valCMy = My/(((valRPM/60)^2)*((valDIA)^5));
valCMx = Mx/(((valRPM/60)^2)*((valDIA)^5));

% Convergence Forces and Moments (average thrust across 1 full rotation)
temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
vecCTCONV(temp)= valCT;
vecCFyCONV(temp)= valFy;
vecCFxCONV(temp) = valFx;

vecCQCONV(temp) = valCQ;
vecCPCONV(temp) = valCP;
vecCMyCONV(temp) = valCMy;
vecCMxCONV(temp) = valCMx;

matDISNORM(:,temp) = vecDISNORM;
matDISTHRUST(:,temp) = vecDISTHRUST;
matDISAXIAL(:,temp) = vecDISAXIAL;
matDISSIDE(:,temp) = vecDISSIDE;

% if valTIMESTEP == 225
%     save('Forces256')
% end
% if valTIMESTEP == 229
%     save('Forces228')
% end
% if valTIMESTEP == 231
%     save('Forces232')
% end
% if valTIMESTEP == 237
%     save('Forces236')
% end

end

