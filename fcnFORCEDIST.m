function [vecLIFTDIST, vecMOMDIST] = fcnFORCEDIST(vecCLDIST,matQTRCRD,vecUINF,vecSPNWSEAREA,valDENSITY,matCLDIST,vecSPNWSECRD)

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

% vecLIFTDIST = 0.5*valDENSITY*norm(vecUINF)*norm(vecUINF).*vecSPNWSEAREA.*vecCLDIST;
vecLIFTDIST = 0.5*valDENSITY*60*60.*vecSPNWSEAREA.*vecCLDIST;

% matMOMDIST = 0.5*valDENSITY*norm(vecUINF)*norm(vecUINF).*vecSPNWSEAREA.*(sum(vecSPNWSECRD)/size(vecSPNWSECRD,1)).*matCLDIST.*matQTRCRD;
matMOMDIST = 0.5*valDENSITY*60*60.*vecSPNWSEAREA.*(sum(vecSPNWSECRD)/size(vecSPNWSECRD,1)).*matCLDIST.*matQTRCRD;

vecMOMDIST = sum(matMOMDIST,2);

end