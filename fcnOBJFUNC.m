function [out] = fcnOBJFUNC(zp)
% clc
% clear

%%
flagRELAX = 1;
flagSTEADY = 1;

valWEIGHT = 7*9.81;
seqALPHA = [4:1:11];
seqBETA = 0;
valKINV = 1.460000e-05;
valDENSITY = 1.2;

valDELTAE = 0;
valDELTIME = 0.05;
valMAXTIME = 30;
valMINTIME = 25;

valINTERF = 20;
%%

valPANELS = 5;

vecSYM = 1;

matGEOM(:,:,1) = [zp(1:5); zp(6:10)];
matGEOM(:,:,2) = [zp(6:10); zp(11:15)];
matGEOM(:,:,3) = [zp(11:15); zp(16:20)];
matGEOM(:,:,4) = [zp(16:20); zp(21:25)];
matGEOM(:,:,5) = [zp(21:25); zp(26:30)];

vecN = [3 3 3 3 3]';
vecM = [1 1 1 1 1]';

vecAIRFOIL = [9 9 9 9 9]';

valVSPANELS = 0;
matVSGEOM = [];
valFPANELS = 0;
matFGEOM = [];
valFTURB = 0;
valFPWIDTH = 0;

%%
valAREA = sum(mean(matGEOM(:,4,:),1).*(matGEOM(2,2,:) - matGEOM(1,2,:)));

valSPAN = matGEOM(2,2,5)*2;

valCMAC = 1;

%% Running VAP2
try
[vecCLv, vecCD, vecCDi, vecVINF, vecCLDIST, matXYZDIST, vecAREADIST, vecDVEAREA] = fcnVAP_MAIN(flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF);
catch
   zp 
end

%% High speed drag coefficient
% Drag coefficient at 51 m/s

highspeed_cd = interp1(vecVINF,vecCD,51,'linear','extrap');

%% Estimating wing weight

rho_cloth = 0.2; % 200 grams per m^2, 6oz carbon fiber

w_cloth = rho_cloth*sum(vecDVEAREA);
w_resin = w_cloth*0.6; % Idealy 60-40 cloth-resin but we know we won't be perfect

w_struct = 1; % 1 kg for misc. structural components? Maybe?

wing_weight = (w_cloth + w_resin)*2 + w_struct; % Top and bottom of the wing in 6oz carbon

%%
drag = vecCD.*0.5.*valDENSITY.*(vecVINF.^2).*valAREA;

preq = drag.*vecVINF;

w_sink = vecVINF.*vecCD./vecCLv;

out = [wing_weight min(preq) min(w_sink) mean(preq(3:7)) mean(w_sink(3:7)) mean(preq) mean(w_sink) highspeed_cd];
%% Writing iteration

fp2 = fopen('optihistory1.txt','at');
fprintf(fp2,'%f %f ', out, zp);
fprintf(fp2,'\r\n');
fclose(fp2);

