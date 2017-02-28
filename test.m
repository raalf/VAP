clc
clear

zp = [0.699743 0.112781 0.216457 0.067837 0.188418 0.309922 7.190567 0.877769 0.051913 -2.729426 0.306655 7.462307 0.868553 0.058878 -4.654932 0.594783 6.951788 0.404617 0.300000 -2.776077 0.595538 7.465613 0.403285 0.058173 -4.346947];
  
[out] = fcnOBJFUNC2(zp)

% z = [...
%     zp(6:10); ... % forward inboard
%     zp(11:15); ... % forward outboard
%     zp(16:20); ... % rear inboard
%     zp(21:25); ... % rear outboard
%     ];
% 
% load('Standard Cirrus Input.mat');
% 
% valMAXTIME = 10;
% seqALPHA = 10;
% 
% out_len = (sqrt(sum(abs(matGEOM(2,1:3,2)-matGEOM(1,1:3,2)).^2)));
% out_vec = (matGEOM(2,1:3,2) - matGEOM(1,1:3,2))./out_len; % Unit vector of outboard leading edge
% 
% lop_loc = out_vec*zp(1); % Where we will cut the wing
% 
% tr = matGEOM(2,4,2)/matGEOM(1,4,2);
% matGEOM(2,4,2) = matGEOM(1,4,2)*tr*(1+(zp(1)/out_len)); % Finding the chord at the cut
% 
% matGEOM(2,5,2) = matGEOM(2,5,2)*(1-(zp(1)/out_len)); % Finding the twist angle at the cut
% matGEOM(2,1:3,2) = matGEOM(2,1:3,2) - lop_loc; % Making the cut
% 
% valPANELS = 9;
% vecAIRFOIL = [1 1 7 6 6 6 6 6 6]';
% vecN = [6 8 3 2 4 4 2 4 4]';
% vecM = [1 1 1 1 1 1 1 1 1]';
% 
% matGEOM(:,:,4) = [matGEOM(2,:,2); [matGEOM(2,1,2) matGEOM(2,2,2) + 0.1 matGEOM(2,3,2) + 0.05 zp(2) zp(3)] ]; % Transition
% matGEOM(:,:,5) = [matGEOM(2,:,4); z(1,:)]; % front inboard
% matGEOM(:,:,6) = [z(1,:); z(2,:)]; % front outboard
% 
% matGEOM(:,:,7) = [matGEOM(2,:,2); [matGEOM(2,1,2) + zp(2) + 0.1 matGEOM(2,2,2) + 0.1 matGEOM(2,3,2) zp(4) zp(5)] ]; % Transition
% matGEOM(:,:,8) = [matGEOM(2,:,7); z(3,:)]; % rear inboard
% matGEOM(:,:,9) = [z(3,:); z(4,:)]; % rear inboard
% 
% [vecCLv, vecCD, vecCDi, vecVINF, vecCLDIST, matXYZDIST, vecAREADIST] = fcnVAP_MAIN(flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF);