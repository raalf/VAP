function [valPANELS, vecN, vecM, vecAIRFOIL, matGEOMNULTI ] = fcnGEOMMULTIROTOR(valNUMRO, valPANELS, vecRODIR, vecROTAXLOC, vecN, vecM, vecAIRFOIL, matGEOM, matROTAX)
% This function applies multiple mutliple rotors to the input geometry


%% Create new geometry for rotors
% Create a local matGEOM wrt the local rotation axis for each rotor
matGEOMLOC = matGEOM(:,1:3,:)-vecROTAXLOC;
matGEOMLOC = repmat(matGEOMLOC, 1, 1, valNUMRO);
tempRODIR = -2*repmat(reshape(repmat((-0.5*(vecRODIR-1)),1,valPANELS)',1,1,valPANELS*valNUMRO),2,1,1)+1;

matGEOMLOC = [matGEOMLOC(:,1,:), tempRODIR.*matGEOMLOC(:,2,:), matGEOMLOC(:,3,:)];


% Create a matrix with the new rotaion axis
%tempROTAX = repmat(permute(reshape(repmat(reshape(permute(matROTAX,[2,1,3]),1,3,valNUMRO),valPANELS,1,1),valPANELS*valNUMRO,3),[3,2,1]),2,1,1);
tempROTAX = repmat(permute(((reshape(repmat(permute(reshape(matROTAX',1,3,valNUMRO),[2,1,3]),1,valPANELS,1),3,valPANELS*valNUMRO,1))'),[3,2,1]),2,1,1);
% Add new new rotation axis with local rotor geometry
matGEOMMULTI = tempROTAX + matGEOMLOC;

% Apply twist and chord distribution
matGEOMNULTI = horzcat(matGEOMMULTI,repmat(matGEOM(:,4:5,:),1,1,valNUMRO));

%% Repmat all required variables
vecN = repmat(vecN,valNUMRO,1,1);
vecM = repmat(vecM,valNUMRO,1,1);
vecAIRFOIL = repmat(vecAIRFOIL,valNUMRO,1,1);
valPANELS = valPANELS*valNUMRO;
end