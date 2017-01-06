function [vecLIFTDIST, vecMOMDIST] = fcnFORCEDIST(vecCLDIST,matQTRCRD,vecSPNWSEAREA,valDENSITY,matCLDIST,vecSPNWSECRD,valVINF,valPANELS)

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

vecLIFTDIST = 0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*vecCLDIST;
temp1 = repmat(vecLIFTDIST,1,2*valPANELS-1);

% matMOMDIST = 0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*(sum(vecSPNWSECRD)/size(vecSPNWSECRD,1)).*matCLDIST.*matQTRCRD;
matMOMDIST = 0.5*valDENSITY*valVINF*valVINF.*vecSPNWSEAREA.*matCLDIST.*matQTRCRD;

vecMOMDIST = sum(matMOMDIST,2);

end