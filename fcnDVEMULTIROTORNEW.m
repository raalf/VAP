function [valDELTIME, vecAZNUM, vecDVEVLSTROTOR, vecDVEROTOR, valNEWNELE, matNEWNPVLST, vecNEWAIRFOIL, vecNEWDVELE, vecNEWDVETE, ...
    vecNEWDVEYAW, vecNEWDVEPANEL, vecNEWDVETIP, vecNEWDVEWING, vecDVESYM, vecNEWM, vecNEWN, ...
    vecNEWDVEROLL, vecNEWDVEAREA, vecNEWDVEPITCH, vecNEWDVEMCSWP, vecNEWDVETESWP, vecNEWDVELESWP, ...
    vecNEWDVEHVCRD, vecNEWDVEHVSPN, vecSYM, vecQARM, matNEWADJE, matNEWCENTER, matNEWVLST, matNEWDVE, matNEWDVENORM] = fcnDVEMULTIROTORNEW(valNUMRO, valNELE, valNUMB, vecDVETIP, vecDVETESWP, vecDVEWING, vecDVEMCSWP, vecM, vecN, vecDVEPANEL, vecDVELESWP, vecDVEHVCRD, vecDVEHVSPN, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, matROTAX, vecAIRFOIL, matNPVLST, matDVE, matADJE, matVLST, vecRPM, valAZNUM)
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
%matNEWVLST = matVLST;
% vecNEWAIRFOIL = vecAIRFOIL;
% vecNEWSYM = vecSYM;
% vecNEWN = vecN;
% vecNEWM = vecM;
vecNEWDVEHVSPN = vecDVEHVSPN;
vecNEWDVEHVCRD = vecDVEHVCRD;
vecNEWDVELESWP = vecDVELESWP;
vecNEWDVEMCSWP = vecDVEMCSWP;
vecNEWDVETESWP = vecDVETESWP;
vecNEWDVEAREA = vecDVEAREA;
vecNEWDVETIP = vecDVETIP;
vecNEWDVELE = vecDVELE;
vecNEWDVETE = vecDVETE; 
matNEWVLST0 = [];
matNEWNPVLST0 = [];
vecDVEROTOR = [];
m = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
for i = 1:valNUMRO
	idxDVEBLADE = vecDVEWING==i;
    matDVEBLADE = matDVE(idxDVEBLADE,:);
    [idxVLSTBLADE,A,B] = unique(matDVEBLADE); % This line is critical thanks to T.D.K.
    
    matNEWVLST0 = [matNEWVLST0;matVLST(idxVLSTBLADE,:)];
    matNEWNPVLST0 = [matNEWNPVLST0;matNPVLST(idxVLSTBLADE,:)];
end

for i = 1:size(matDVE,1)
    for j = 1:size(matDVE,2)
        tempVLST = matVLST(matDVE(i,j),:);
        idx = (sum(tempVLST(1,:) == matNEWVLST0(:,:),2) == 3);
        matNEWDVE0(i,j) = find(idx == 1);
    end
end
matNEWVLST = matNEWVLST0;
matNEWNPVLST = matNEWNPVLST0;
matNEWDVE = matNEWDVE0;

% clf
for j = 1:(valNUMB-1)
for i = 1:valNUMRO
    idxDVEBLADE = vecDVEWING==i;
    matDVEBLADE = matNEWDVE0(idxDVEBLADE,:);
    [idxVLSTBLADE,~,~] = unique(matDVEBLADE); % This line is critical thanks to T.D.K.
    vecROTAX = matROTAX(i,:);
    
    tempVLST = matNEWVLST0(idxVLSTBLADE,:) - vecROTAX;
    tempNPVLST = matNEWNPVLST0(idxVLSTBLADE,:) - vecROTAX;
    
% Rotate values
% tempNEWCENTER =(tempROTATE2D*tempCENTER')';
tempNEWVLST = (tempROTATE(:,:,j)*tempVLST')';
tempNEWNPVLST = (tempROTATE(:,:,j)*tempNPVLST')';
% tempDVENORM = (tempROTATE2D*matDVENORM')';

% Calculated below using fcnDVECORNER2PARAM now
% Reshape into orgininal matCENTER format
% temp = reshape(tempNEWCENTER,[numel(tempCENTER)/3,3,valNUMB]);
% matNEWCENTER = reshape(permute(temp,[2,1,3]),[3,valNUMB*(numel(tempCENTER)/3)])' + vecROTAX;

% temp = reshape(tempNEWVLST,[numel(tempVLST)/3,3,valNUMB-1]);
% tempVLSTADD = reshape(permute(temp,[2,1,3]),[3,(valNUMB-1)*(numel(tempVLST)/3)])' + vecROTAX;
% matNEWVLST = [matNEWVLST; tempVLSTADD];
matNEWVLST = [matNEWVLST; tempNEWVLST+vecROTAX];
% d = tempNEWVLST+vecROTAX;
% d2 = matNEWVLST0(idxVLSTBLADE,:);
% hold on
% plot3(d(:,1),d(:,2),d(:,3),'d','Color',m(i,:))
% plot3(d2(:,1),d2(:,2),d2(:,3),'d','Color',m(i,:))
% axis equal
% temp = reshape(tempNEWNPVLST,[numel(tempNPVLST)/3,3,(valNUMB-1)]);
% tempVLSTADD = reshape(permute(temp,[2,1,3]),[3,(valNUMB-1)*(numel(tempNPVLST)/3)])' + vecROTAX;
% matNEWNPVLST = [matNEWNPVLST; tempVLSTADD];
matNEWNPVLST = [matNEWNPVLST; tempNEWNPVLST+vecROTAX];

%% Parameters to increase vector for number of blades
% vecNEWAIRFOIL = [vecNEWAIRFOIL; repmat(vecAIRFOIL(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWN = [vecNEWN; repmat(vecN(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWM = [vecNEWM; repmat(vecM(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVEHVSPN =[vecNEWDVEHVSPN; repmat(vecDVEHVSPN(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVEHVCRD =[vecNEWDVEHVCRD; repmat(vecDVEHVCRD(idxDVEBLADE),[valNUMB-1,1])]; 
% vecNEWDVELESWP = [vecNEWDVELESWP; repmat(vecDVELESWP(idxDVEBLADE),[valNUMB-1,1])]; 
% vecNEWDVEMCSWP = [vecNEWDVEMCSWP; repmat(vecDVEMCSWP(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVETESWP = [vecNEWDVETESWP; repmat(vecDVETESWP(idxDVEBLADE),[valNUMB-1,1])];

% vecNEWAIRFOIL = [vecNEWAIRFOIL;vecAIRFOIL(idxDVEBLADE)];
% vecNEWN = [vecNEWN; vecN(idxDVEBLADE)];
% vecNEWM = [vecNEWM; vecM(idxDVEBLADE)];
vecNEWDVEHVSPN =[vecNEWDVEHVSPN; vecDVEHVSPN(idxDVEBLADE)];
vecNEWDVEHVCRD =[vecNEWDVEHVCRD; vecDVEHVCRD(idxDVEBLADE)]; 
vecNEWDVELESWP = [vecNEWDVELESWP; vecDVELESWP(idxDVEBLADE)]; 
vecNEWDVEMCSWP = [vecNEWDVEMCSWP; vecDVEMCSWP(idxDVEBLADE)];
vecNEWDVETESWP = [vecNEWDVETESWP; vecDVETESWP(idxDVEBLADE)];


% vecNEWDVEAREA = [vecNEWDVEAREA; repmat(vecDVEAREA(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVETIP = [vecNEWDVETIP; repmat(vecDVETIP(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVELE = [vecNEWDVELE; repmat(vecDVELE(idxDVEBLADE),[valNUMB-1,1])];
% vecNEWDVETE = [vecNEWDVETE; repmat(vecDVETE(idxDVEBLADE),[valNUMB-1,1])];

vecNEWDVEAREA = [vecNEWDVEAREA;vecDVEAREA(idxDVEBLADE)];
vecNEWDVETIP = [vecNEWDVETIP; vecDVETIP(idxDVEBLADE)];
vecNEWDVELE = [vecNEWDVELE; vecDVELE(idxDVEBLADE)];
vecNEWDVETE = [vecNEWDVETE; vecDVETE(idxDVEBLADE)];


%% Parameters to be multiplied by constant

% New matDVE list
% tempADDI = ((reshape(repmat(1:(valNUMB-1),size(matDVE(idxDVEBLADE,:),1),4),[(valNUMB-1)*size(matDVE(idxDVEBLADE,:),1) 4])));
% temp = repmat(matDVE(idxDVEBLADE,:),[valNUMB-1 1]);
% for j = 1:(valNUMB-1)
%     tempD = (1+temp-min(temp(tempADDI==j))).*(tempADDI == j)+(tempADDI == j)*(max(max(matNEWDVE)));
%     matNEWDVE = [matNEWDVE;tempD(any(tempD,2),:)]
% end
% hold on
% if i == 1
% scatter3(matNEWVLST(matNEWDVE(1:8,:),1),matNEWVLST(matNEWDVE(1:8,:),2),matNEWVLST(matNEWDVE(1:8,:),3),'k')
% elseif i == 2
% scatter3(matNEWVLST(matNEWDVE(9:16,:),1),matNEWVLST(matNEWDVE(9:16,:),2),matNEWVLST(matNEWDVE(9:16,:),3),'r')
% elseif i == 3
% scatter3(matNEWVLST(17:24,1),matNEWVLST(17:24,2),matNEWVLST(17:24,3),'b')
% elseif i == 4
% scatter3(matNEWVLST(25:32,1),matNEWVLST(25:32,2),matNEWVLST(25:32,3),'g')
% end
end
matNEWDVE = [matNEWDVE; matNEWDVE0+max(matNEWDVE(:))];
end

vecNEWAIRFOIL = repmat(vecAIRFOIL,valNUMB,1);
vecNEWM = repmat(vecM, valNUMB,1);
vecNEWN = repmat(vecN, valNUMB,1);

% Define which DVE is on which rotor
vecDVEROTOR = repmat(vecDVEWING,valNUMB,1);

% Create vecDVEVLSTROTOR to define which VLST point is associated to which
% rotor
vecDVEVLSTROTOR = zeros(size(matNEWVLST,1),1);
for i = 1:size(matNEWDVE,1)
    for j = 1:size(matNEWDVE,2)
        tempVLST = matNEWVLST(matNEWDVE(i,j),:);
        idx = (sum(tempVLST(1,:) == matNEWVLST(:,:),2) == 3);
        vecDVEVLSTROTOR(idx) = vecDVEROTOR(i);
    end
end

% New matADJE
if isempty(matADJE) == 1 % If there is only one panel per wing, matADJE will be empty
    tempADDI = (reshape(repmat(1:valNUMB,numel(matADJE)/4,2),[valNUMB*numel(matADJE)/4 2]))*0;
else
    tempADDI = (reshape(repmat(1:valNUMB,numel(matADJE)/4,2),[valNUMB*numel(matADJE)/4 2]))*(max(max(matADJE(:,1),max(matADJE(:,3)))))- (max(max(matADJE(:,1),max(matADJE(:,3)))));
end
tempADJEADD = repmat([matADJE(:,1) matADJE(:,3)],[valNUMB 1]) + tempADDI;
tempADJDUP = repmat([matADJE(:,2) matADJE(:,4)],[valNUMB,1]);
    
matNEWADJE = [tempADJEADD(:,1),tempADJDUP(:,1),tempADJEADD(:,2),tempADJDUP(:,2)];

% New vecDVEWING
tempADDI = (reshape(repmat(1:valNUMB,numel(vecDVEWING),1),[valNUMB*numel(vecDVEWING) 1]))*(max(max(vecDVEWING))) - max(max(vecDVEWING));
vecNEWDVEWING = repmat(vecDVEWING,[valNUMB 1])+tempADDI;

% New vecDVEPANEL
tempADDI = (reshape(repmat(1:valNUMB,numel(vecDVEPANEL),1),[valNUMB*numel(vecDVEPANEL) 1]))*(max(max(vecDVEPANEL))) - max(max(vecDVEPANEL));
vecNEWDVEPANEL = repmat(vecDVEPANEL,[valNUMB 1])+tempADDI;

valNEWNELE = valNELE*valNUMB;

%% Updating roll, pitch, yaw, etc
[~, ~, vecNEWDVEROLL, vecNEWDVEPITCH, vecNEWDVEYAW, ~, ~, ~, ~, matNEWDVENORM, ...
    ~, ~, matNEWCENTER] = fcnVLST2DVEPARAM(matNEWDVE, matNEWVLST);

% Chordwise radial distances
vecQARM = abs(matNEWCENTER(:,2)-matROTAX(vecDVEROTOR,2));


% Calculate a del time using information from first rotor
valDELTIME = 1/((vecRPM(1)/60)*(valAZNUM));

vecAZNUM = 1./((vecRPM./60).*valDELTIME);

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

