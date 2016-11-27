clc
clear

zp = [0.218004 0.100000 0.525675 0.196431 0.834489 0.300000 7.432013 0.541958 0.103197 -0.863533 0.404563 7.482095 0.786142 0.040076 -0.646307 0.52 7.432051 0.393546 0.15 -2.547441 0.60 7.482120 0.393542 0.10 -0.222763];


[out] = fcnOBJFUNC(zp)

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