clc
clear

   zp = [0.624926261540485   0.100003187623566   0.488206990962351 ...
   0.101765181067386   1.379755063829769   0.448021616045744 ...
   6.926771038990306   0.425596447324849   0.102688960691146 ...
  -0.657720100180761   0.509690005977665   6.988465177690766 ...
   0.421884717910574   0.040307863495002   1.203540165231571 ...
   0.599970714438376   6.925537065588152   0.364500401041310 ...
   0.102583187307354   0.194383972203211   0.599997869472668 ...
   6.981743110780286   0.362682086697931   0.149425757244740 ...
   1.441955913052306]

% zp(1) = lop distance;


z = [...
    zp(6:10); ... % forward inboard
    zp(11:15); ... % forward outboard
    zp(16:20); ... % rear inboard
    zp(21:25); ... % rear outboard
    ];

load('Standard Cirrus Input.mat');

valMAXTIME = 10;
seqALPHA = 10;

out_len = (sqrt(sum(abs(matGEOM(2,1:3,2)-matGEOM(1,1:3,2)).^2)));
out_vec = (matGEOM(2,1:3,2) - matGEOM(1,1:3,2))./out_len; % Unit vector of outboard leading edge

lop_loc = out_vec*zp(1); % Where we will cut the wing

tr = matGEOM(2,4,2)/matGEOM(1,4,2);
matGEOM(2,4,2) = matGEOM(1,4,2)*tr*(1+(zp(1)/out_len)); % Finding the chord at the cut

matGEOM(2,5,2) = matGEOM(2,5,2)*(1-(zp(1)/out_len)); % Finding the twist angle at the cut
matGEOM(2,1:3,2) = matGEOM(2,1:3,2) - lop_loc; % Making the cut

valPANELS = 9;
vecAIRFOIL = [1 1 7 6 6 6 6 6 6]';
vecN = [6 8 3 2 4 4 2 4 4]';
vecM = [1 1 1 1 1 1 1 1 1]';

matGEOM(:,:,4) = [matGEOM(2,:,2); [matGEOM(2,1,2) matGEOM(2,2,2) + 0.1 matGEOM(2,3,2) + 0.05 zp(2) zp(3)] ]; % Transition
matGEOM(:,:,5) = [matGEOM(2,:,4); z(1,:)]; % front inboard
matGEOM(:,:,6) = [z(1,:); z(2,:)]; % front outboard

matGEOM(:,:,7) = [matGEOM(2,:,2); [matGEOM(2,1,2) + zp(2) + 0.1 matGEOM(2,2,2) + 0.1 matGEOM(2,3,2) zp(4) zp(5)] ]; % Transition
matGEOM(:,:,8) = [matGEOM(2,:,7); z(3,:)]; % rear inboard
matGEOM(:,:,9) = [z(3,:); z(4,:)]; % rear inboard

[vecCLv, vecCD, vecCDi, vecVINF, vecCLDIST, matXYZDIST, vecAREADIST] = fcnVAP_MAIN(flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF);