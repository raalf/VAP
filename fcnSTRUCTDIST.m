function [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecSPANDIST] = fcnSTRUCTDIST(vecDVEHVSPN,vecDVELE,vecEIxCOEFF,vecGJtCOEFF,vecEACOEFF,vecCGCOEFF, vecJTCOEFF, vecLMCOEFF)

% Find DVEs on LE of wing
[ledves, ~, ~] = find(vecDVELE > 0);

% Determine spanwise location (y coordinate) of DVEs
tempSPANDIST = 2*vecDVEHVSPN(ledves);
matSPANDIST = repmat(tempSPANDIST,1,size(tempSPANDIST));

temp = triu(matSPANDIST);

vecSPANDIST = sum(temp(:,1:size(matSPANDIST,2)))';
vecSPANDIST = [0; vecSPANDIST];

% Spanwise bending stiffness distribution. Cols 2 and 3 are the first and
% second derivatives
matEIx(:,1) = (vecEIxCOEFF(1).*vecSPANDIST.^2 + vecEIxCOEFF(2).*vecSPANDIST + vecEIxCOEFF(3))';
matEIx(:,2) = (2*vecEIxCOEFF(1).*vecSPANDIST + vecEIxCOEFF(2))';
matEIx(:,3) = (repmat(2*vecEIxCOEFF(1),1,size(vecSPANDIST)))';

% Spanwise torsional stiffness distribution. Col 2 is the first derivative
matGJt(:,1) = (vecGJtCOEFF(1).*vecSPANDIST.^2 + vecGJtCOEFF(2).*vecSPANDIST + vecGJtCOEFF(3))';
matGJt(:,2) = (2*vecGJtCOEFF(1).*vecSPANDIST + vecGJtCOEFF(2))';

vecEA = vecEACOEFF(1).*vecSPANDIST.^2 + vecEACOEFF(2).*vecSPANDIST + vecEACOEFF(3);
vecCG = vecCGCOEFF(1).*vecSPANDIST.^2 + vecCGCOEFF(2).*vecSPANDIST + vecCGCOEFF(3);
vecJT = vecJTCOEFF(1).*vecSPANDIST.^2 + vecJTCOEFF(2).*vecSPANDIST + vecJTCOEFF(3);
vecLM = vecLMCOEFF(1).*vecSPANDIST.^2 + vecLMCOEFF(2).*vecSPANDIST + vecLMCOEFF(3);

vecLSM = vecEA - vecCG;

end


