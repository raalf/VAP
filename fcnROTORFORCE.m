function [ thrustind, thrustfree, thrustCFfree, tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE] = fcnROTORFORCE(difthrustP, diffsideP, diffaxialP, nind, nfree, nfreecs, Pthrust, thrustind, thrustfree, thrustCFfree,  thrustinddrag, Ptorque, axialCFfree, axialind, axialfree, inddrag, sideind, sidefree, sideCFfree, sideinddrag, valRPM, valDIA, valAZNUM, valTIMESTEP, vecCTCONV, vecCFxCONV, vecCFyCONV, vecCQCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecQARM, vecDVETE, vecTHETA, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, valNUMRO, vecDVEROTOR)
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

for i = 1:valNUMRO
    idx = vecDVEROTOR == i;
    idx2 = (vecDVETE==3 & vecDVEROTOR == i);
    % Calculate total force and moment values per density
    thrust(i) = sum(thrustind(idx))+sum(thrustfree(idx))+sum(thrustCFfree(idx)) + sum(thrustinddrag(idx)) + sum(Pthrust(idx)) + sum(difthrustP(idx));
    Fy(i) = sum(vecDISSIDE(idx).*sin(vecTHETA(idx)) + vecDISAXIAL(idx).*sin(pi - vecTHETA(idx)));
    Fx(i) = sum(vecDISSIDE(idx).*cos(vecTHETA(idx)) + vecDISAXIAL(idx).*cos(pi - vecTHETA(idx))); % ADD viscous

    torque(i) = sum(axialind(idx).*vecQARM(idx))+sum(axialfree(idx).*vecQARM(idx))+sum(inddrag(idx).*vecQARM(idx2)) + sum(axialCFfree(idx).*vecQARM(idx)) + sum(Ptorque(idx)) + sum(diffaxialP(idx).*vecQARM(idx)); %sum(diffaxialP);
    
    Mx(i) = sum(vecDISTHRUST(idx).*(vecQARM(idx).*sin(vecTHETA(idx))));
    My(i) = sum(vecDISTHRUST(idx).*(vecQARM(idx).*cos(vecTHETA(idx))));

end

power = torque.*2.*pi.*(valRPM./60);
% Calculate non-dimensionalized coefficient using US customary definitions
% CT = T/(rho*Omega^2*D^4)
valCT = thrust./(((valRPM/60)^2)*((valDIA)^4));
valFy = Fy./(((valRPM/60)^2)*((valDIA)^4));
valFx = Fx./(((valRPM/60)^2)*((valDIA)^4));

valCQ = torque./(((valRPM/60)^2)*((valDIA)^5));
valCP = power./((valRPM/60)^3*(valDIA^5));
valCMy = My./(((valRPM/60)^2)*((valDIA)^5));
valCMx = Mx./(((valRPM/60)^2)*((valDIA)^5));

% Convergence Forces and Moments (average thrust across 1 full rotation)
temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
vecCTCONV(temp,:)= valCT;
vecCFyCONV(temp,:)= valFy;
vecCFxCONV(temp,:) = valFx;

vecCQCONV(temp,:) = valCQ;
vecCPCONV(temp,:) = valCP;
vecCMyCONV(temp,:) = valCMy;
vecCMxCONV(temp,:) = valCMx;

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

