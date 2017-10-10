function [valNELE, matNEWNPVLST, vecAIRFOIL, vecDVELE, vecDVETE, ...
    vecDVEYAW, vecDVEPANEL, vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, ...
    vecDVEROLL, vecDVEAREA, vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, vecDVELESWP, ...
    vecDVEHVCRD, vecDVEHVSPN, vecSYM, vecQARM, matADJE, matNEWCENTER, matNEWVLST, matDVE, matNEWDVENORM, matVLST] = fcnDVEMULTIROTORNEW(valNUMRO, valNELE, valNUMB, vecDVETIP, vecDVETESWP, vecDVEPITCH, vecDVEWING, vecDVEMCSWP, vecM, vecN, vecDVEPANEL, vecDVEROLL, vecDVELESWP, vecDVEYAW, vecDVEHVCRD, vecDVEHVSPN, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, matROTAX, vecAIRFOIL, matNPVLST, matDVE, matADJE, matVLST, matCENTER, matDVENORM)
% This function modifies the created DVEs and all required input values
%	for multiple rotor blades.
%
% Inputs:
%   Parameters from input file and fcnGENERATEDVES that require
%   modification
%
% Outputs:
%   Each input has been modified for the given number of blades:
%   From input file
%       valAREA - multiply by valNUMB
%       vecSYM - Add (will be deleted)
%       vecAIRFOIL - add
%       vecN - add
%       venM - add
%
%   From Generate DVEs
%		matCENTER - rotate values center values
%		vecDVEHVSPN - add to vector for each blade
%		vecDVEHVCRD - add to vector for each blade
%		vecDVELESWP - add to vector for each blade
%		vecDVEMCSWP - add to vector for each blade
%		vecDVETESWP - add to vector for each blade
%		vecDVEROLL - add to vector for each blade
%		vecDVEPITCH - add to vector for each blade
%		vecDVEYAW - add to vector for each blade
%		vecDVEAREA - add to vector for each blade
%		matDVENORM - rotate values normal vectors
%		matVLST - rotate vertices
%		valNELE - multiply by valNUMB
%		matDVE - matDVE + max(matDVE) for each blade
%		matADJE - matADJE + max(matADJE) for each blade
%		vecDVESYM - delete in later versions
%		vecDVETIP - add to vector for each blade
%		vecDVEWING - vecDVEWING + max(vecDVEWING) for each blade
%		vecDVELE - add to vector for each blade
%		vecDVETE - add to vector for each blade
%		vecDVEPANEL - vecDVEPANEL + max(vecDVEPANEL) for each blade


%% Create rotation matrix
% Angle of each blade from original (+ve CCW & in rad)
tempTHETA = (0:2*pi/valNUMB:2*pi)';
tempTHETA(valNUMB+1) = [];
tempTHETA(1) = [];

tempTOP = [cos(tempTHETA), -sin(tempTHETA), zeros([valNUMB-1,1])];
tempMID = [sin(tempTHETA), cos(tempTHETA), zeros([valNUMB-1,1])];
tempLOW = [zeros([valNUMB-1, 1]), zeros([valNUMB-1, 1]), ones([valNUMB-1, 1])];
tempROTATE = reshape(tempTOP',[1, 3, valNUMB-1]);
tempROTATE(2,:,:) = reshape(tempMID',[1, 3, valNUMB-1]);
tempROTATE(3,:,:) = reshape(tempLOW',[1, 3, valNUMB-1]);

% 2D rotation matrix
tempROTATE2D = (reshape(permute(tempROTATE,[2,1,3]),[3 (valNUMB-1)*3]))';


%% Parameters the must Rotate
% Make Each point relative to rotation axis
% tempCENTER = matCENTER - vecROTAX;
matNEWVLST = matVLST;
matNEWNPVLST = matNPVLST;
vecNEWAIRFOIL = vecAIRFOIL;
vecNEWSYM = vecSYM;
vecNEWN = vecN;
vecNEWM = vecM;
vecNEWDVEHVSPN = vecDVEHVSPN;
vecNEWDVEHVCRD = vecDVEHVCRD;
vecNEWDVELESWP = vecDVELESWP;
vecNEWDVEMCSWP = vecDVEMCSWP;
vecNEWDVETESWP = vecDVETESWP;
vecNEWDVEAREA = vecDVEAREA;
vecNEWDVETIP = vecDVETIP;
vecNEWDVELE = vecDVELE;
vecNEWDVETE = vecDVETE; 
matNEWDVE = matDVE;

for i = 1:valNUMRO
    idxDVEBLADE = vecDVEWING==i;
    matDVEBLADE = matDVE(idxDVEBLADE,:);
    [idxVLSTBLADE,~,c] = unique(matDVEBLADE); % This line is critical thanks to T.D.K.
    vecROTAX = matROTAX(i,:);
    
    tempVLST = matVLST(idxVLSTBLADE,:) - vecROTAX;
    tempNPVLST = matNPVLST(idxVLSTBLADE,:) - vecROTAX;
    
% Rotate values
% tempNEWCENTER =(tempROTATE2D*tempCENTER')';
tempNEWVLST = (tempROTATE2D*tempVLST')';
tempNEWNPVLST = (tempROTATE2D*tempNPVLST')';
% tempDVENORM = (tempROTATE2D*matDVENORM')';

% Calculated below using fcnDVECORNER2PARAM now
% Reshape into orgininal matCENTER format
% temp = reshape(tempNEWCENTER,[numel(tempCENTER)/3,3,valNUMB]);
% matNEWCENTER = reshape(permute(temp,[2,1,3]),[3,valNUMB*(numel(tempCENTER)/3)])' + vecROTAX;

temp = reshape(tempNEWVLST,[numel(tempVLST)/3,3,valNUMB-1]);
tempVLSTADD = reshape(permute(temp,[2,1,3]),[3,(valNUMB-1)*(numel(tempVLST)/3)])' + vecROTAX;
matNEWVLST = [matNEWVLST; tempVLSTADD];
temp = reshape(tempNEWNPVLST,[numel(tempNPVLST)/3,3,(valNUMB-1)]);
tempVLSTADD = reshape(permute(temp,[2,1,3]),[3,(valNUMB-1)*(numel(tempNPVLST)/3)])' + vecROTAX;
matNEWNPVLST = [matNEWNPVLST; tempVLSTADD];

%% Parameters to increase vector for number of blades
vecNEWAIRFOIL = [vecNEWAIRFOIL; repmat(vecAIRFOIL(idxDVEBLADE),[valNUMB-1,1])];
vecNEWN = [vecNEWN; repmat(vecN(idxDVEBLADE),[valNUMB-1,1])];
vecNEWM = [vecNEWM; repmat(vecM(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVEHVSPN =[vecNEWDVEHVSPN; repmat(vecDVEHVSPN(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVEHVCRD =[vecNEWDVEHVCRD; repmat(vecDVEHVCRD(idxDVEBLADE),[valNUMB-1,1])]; 
vecNEWDVELESWP = [vecNEWDVELESWP; repmat(vecDVELESWP(idxDVEBLADE),[valNUMB-1,1])]; 
vecNEWDVEMCSWP = [vecNEWDVEMCSWP; repmat(vecDVEMCSWP(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVETESWP = [vecNEWDVETESWP; repmat(vecDVETESWP(idxDVEBLADE),[valNUMB-1,1])];

vecNEWDVEAREA = [vecNEWDVEAREA; repmat(vecDVEAREA(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVETIP = [vecNEWDVETIP; repmat(vecDVETIP(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVELE = [vecNEWDVELE; repmat(vecDVELE(idxDVEBLADE),[valNUMB-1,1])];
vecNEWDVETE = [vecNEWDVETE; repmat(vecDVETE(idxDVEBLADE),[valNUMB-1,1])];

%% Parameters to be multiplied by constant

% New matDVE list
tempADDI = ((reshape(repmat(1:(valNUMB-1),size(matDVE(idxDVEBLADE,:),1),4),[(valNUMB-1)*size(matDVE(idxDVEBLADE,:),1) 4]))*(max(max(matNEWDVE))));
matNEWDVE = [matNEWDVE; repmat(matDVE(idxDVEBLADE,:),[valNUMB-1 1])+tempADDI]; 
end

% New matADJE
if isempty(matADJE) == 1 % If there is only one panel per wing, matADJE will be empty
    tempADDI = (reshape(repmat(1:valNUMB,numel(matADJE)/4,2),[valNUMB*numel(matADJE)/4 2]))*0;
else
    tempADDI = (reshape(repmat(1:valNUMB,numel(matADJE)/4,2),[valNUMB*numel(matADJE)/4 2]))*(max(max(matADJE(:,1),max(matADJE(:,3)))))- (max(max(matADJE(:,1),max(matADJE(:,3)))));
end
tempADJEADD = repmat([matADJE(:,1) matADJE(:,3)],[valNUMB 1]) + tempADDI;
tempADJDUP = repmat([matADJE(:,2) matADJE(:,4)],[valNUMB,1]);
    
matADJE = [tempADJEADD(:,1),tempADJDUP(:,1),tempADJEADD(:,2),tempADJDUP(:,2)];

% New vecDVEWING
tempADDI = (reshape(repmat(1:valNUMB,numel(vecDVEWING),1),[valNUMB*numel(vecDVEWING) 1]))*(max(max(vecDVEWING))) - max(max(vecDVEWING));
vecDVEWING = repmat(vecDVEWING,[valNUMB 1])+tempADDI;

% New vecDVEPANEL
tempADDI = (reshape(repmat(1:valNUMB,numel(vecDVEPANEL),1),[valNUMB*numel(vecDVEPANEL) 1]))*(max(max(vecDVEPANEL))) - max(max(vecDVEPANEL));
vecDVEPANEL = repmat(vecDVEPANEL,[valNUMB 1])+tempADDI;

valNELE = valNELE*valNUMB;

%% Updating roll, pitch, yaw, etc
[~, ~, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ~, ~, ~, ~, matNEWDVENORM, ...
    ~, ~, matNEWCENTER] = fcnVLST2DVEPARAM( matDVE, matNEWVLST);

% Chordwise radial distances
vecQARM = abs(matNEWCENTER(:,2)-vecROTAX(2));

%% Scatter plot of centers and verticies for validation
% figure(1)
% clf(1)
% hold on
% scatter3(matNEWCENTER(:,1),matNEWCENTER(:,2),matNEWCENTER(:,3),'+')
% scatter3(matNEWVLST(:,1),matNEWVLST(:,2),matNEWVLST(:,3),'*')
% quiver3(matNEWCENTER(:,1),matNEWCENTER(:,2),matNEWCENTER(:,3), matNEWDVENORM(:,1),matNEWDVENORM(:,2),matNEWDVENORM(:,3))
% axis equal
% xlabel('X-Dir','FontSize',15);
% ylabel('Y-Dir','FontSize',15);
% zlabel('Z-Dir','FontSize',15);
% box on
% grid on
% axis equal
% axis tight
% hold off

end

