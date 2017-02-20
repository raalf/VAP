function [vecLIFTDIST, vecMOMDIST, valVINF] = fcnFORCEDIST(vecCLDIST,matQTRCRD,vecSPNWSEAREA,valDENSITY,matCLDIST,valWEIGHT,valCL,vecDVEHVSPN,vecLEDVES)

% This function computes the dimensional force and moment distribution
% across the wing, resolved to the aerodynamic center line
%
% INPUT:
%
% OUTPUT:
% vecLIFTDIST - 1 x sum(vecN) matrix of the total lift at each spanwise
% station
% vecMOMDIST - 1 x sum(vecN) matrix of the total pitching moment at each
% spanwise station (+ve nose up)

% Calculate Vinf required for steady level flight. This will be used for
% load calculations
% q_inf = valWEIGHT/(valCL*sum(vecSPNWSEAREA));
% valVINF = sqrt(2*q_inf/valDENSITY);
valVINF = 30;

vecLIFTDIST = (2*0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*vecCLDIST)./(2*vecDVEHVSPN(vecLEDVES));

% matMOMDIST = 0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*(sum(vecSPNWSECRD)/size(vecSPNWSECRD,1)).*matCLDIST.*matQTRCRD;
matMOMDIST = (0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*matCLDIST.*matQTRCRD)./(2*vecDVEHVSPN(vecLEDVES));

vecMOMDIST = sum(matMOMDIST,2);

end