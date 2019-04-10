clc
clear
tic 
warning off

disp('===========================================================================');
disp('+---------------+');
disp('| RYERSON       |       VAPTOR (Based on FreeWake 2015)');
disp('| APPLIED       | ');
disp('| AERODYNAMICS  | ');
disp('| LABORATORY OF | ');
disp('| FLIGHT        |        ---------------+--------------- ');
disp('+---------------+                  ___ /^^[___              _');
disp('                                  /|^+----+   |#___________//');
disp('                                ( -+ |____|    ______-----+/');
disp('                                  ==_________--            \');
disp('                                    ~_|___|__');
disp(' ');
disp('===========================================================================');
disp(' ');

%% User Inputs
% strFILE is a string with the input filename
% strFILE = 'inputs/TMotor.txt';
strFILE = 'inputs/TMotorQuad.txt';

% Optional flags for more options
flagVISCOUS = 1; % Apply viscous effects
flagPRINT   = 1; % Print out results into command window
flagPLOT    = 1; % Plot the rotor at the end of the run
flagVERBOSE = 0; % Add numbers to plot
flagPROGRESS = 1; % Use a progress bar when running performance sweeps
flagHOVERWAKE = 1; % Propogate wake downward in hover simulations
flagSAVE = 1; % Save results
filename = 'Results'; % Save workspace name

% Adding a collective pitch
valCOLLECTIVE = 0;

%% Reading in geometry
[flagRELAX, flagSTEADY, valMAXTIME, valAZNUM, seqALPHAR, seqJ, vecRPM, ...
    valDENSITY, valKINV, valDIA, valNUMB, valNUMRO, matROTAX0, ...
    vecRODIR, valPANELS, vecROTAXLOC, matGEOM, vecAIRFOIL, vecN, vecM, ...
    vecSYM] = fcnVAPTORREADMULTI(strFILE);

matGEOM(:,end,:) = matGEOM(:,end,:) + valCOLLECTIVE;
%% Apply multiple rotors
[valPANELS, vecN, vecM, vecAIRFOIL, matGEOM ] = ...
    fcnGEOMMULTIROTOR(valNUMRO, valPANELS, vecRODIR, vecROTAXLOC, vecN, ...
    vecM, vecAIRFOIL, matGEOM, matROTAX0);

%% Discretize geometry into DVEs
[~, vecDVEHVSPN0, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, ...
    vecDVETESWP, ~, ~, ~, vecDVEAREA, ...
    ~, matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, ...
    vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

%% Apply multiple blades
[valDELTIME, vecAZNUM, vecDVEVLSTROTOR, vecDVEROTOR, valNELE, matNPVLST0, vecAIRFOIL, vecDVELE, vecDVETE, ...
    vecDVEYAW, vecDVEPANEL, vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, ...
    vecDVEROLL, vecDVEAREA, vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, vecDVELESWP, ...
    vecDVEHVCRD, vecDVEHVSPN0, vecSYM, vecQARM, matADJE, matCENTER0, matVLST0, matDVE, ...
    matDVENORM0] = fcnDVEMULTIROTORNEW(valNUMRO, valNELE, valNUMB, ...
    vecDVETIP, vecDVETESWP, vecDVEWING, vecDVEMCSWP, vecM, vecN, ...
    vecDVEPANEL, vecDVELESWP, vecDVEHVCRD, vecDVEHVSPN0, vecDVEAREA, ...
    vecDVESYM, vecDVELE, vecDVETE, vecSYM, matROTAX0, vecAIRFOIL, ...
    matNPVLST0, matDVE, matADJE, matVLST0, vecRPM, valAZNUM);

valWSIZE = length(nonzeros(vecDVETE)); % Amount of wake DVEs shed each timestep
%% Add boundary conditions to D-Matrix
[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN0, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN0);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, ...
    matDVENORM0, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, ...
    vecDVETESWP, vecDVEHVSPN0, vecDVEHVCRD, vecSYM);

%% Performance sweeps
% Preallocating for a turbo-boost in performance
OUTPUT.CT = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CFx = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CFy = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CQ = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CMx = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.CMy = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCTTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCPTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCFxTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCFyTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCQTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCMxTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
OUTPUT.matCMyTIMESTEP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matthrustind = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
matthrustfree = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
matthrustCFfree = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
mattempTi = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
matPthrust = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
matPtorque = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
matdifthrustP = zeros(valNELE,valMAXTIME,length(seqJ),length(seqALPHAR));
OUTPUT.seqJ = seqJ;
OUTPUT.seqALPHAR = seqALPHAR;

if flagPROGRESS == 1 % Apply progress bar
    progressbar('Angle of Attack','Advance Ratio','Timestep')
end

for ai = 1:length(seqALPHAR)
    for ji= 1:length(seqJ)
        valALPHAR = (pi/180)*(seqALPHAR(ai));
        valJ = seqJ(ji);

        % This is done for when we are using a parfor loop
        matCENTER = matCENTER0;
        matVLST = matVLST0;
        matNPVLST = matNPVLST0;
        matROTAX = matROTAX0;
        vecDVEHVSPN = vecDVEHVSPN0;
        matDVENORM = matDVENORM0;
        fprintf('      ANGLE OF ATTACK = %0.3f DEG\n',seqALPHAR(ai));
        fprintf('      ADVANCE RATIO = %0.3f\n',valJ);
        fprintf('      RPM = ')
        fprintf('%0.f ',vecRPM);
        fprintf('\n\n')
        fprintf('-------------------');
        for i = 1:valNUMRO
        fprintf('---------------------');
        end
        fprintf('\n                ');
        fprintf('Rotor %d                 ',1:valNUMRO);
        fprintf('\n');
        fprintf(' Timestep  ');
        for i = 1:valNUMRO
            fprintf(' Thurst   Power         ');
        end
        fprintf('\n');      
        fprintf('-------------------');
        for i = 1:valNUMRO
        fprintf('---------------------');
        end
        fprintf('\n');

        % Calculate inflow velocity at each control point    
         [matUINF, matUINFTE, matTEPTS, vecTHETA] = fcnUINFMULTIROT( ...
             matCENTER, matROTAX, 0, vecRPM, valALPHAR, vecAZNUM, valDIA, ...
             valJ, valNUMB, vecDVEROTOR, vecDVEHVSPN, vecDVETE, matVLST, ...
             matDVE, vecRODIR);

        % Initializing wake parameters
        matWAKEGEOM = [];
        matNPWAKEGEOM = [];
        vecWDVEHVSPN = [];
        vecWDVEHVCRD = [];
        vecWDVEROLL = [];
        vecWDVEPITCH = [];
        vecWDVEYAW = [];
        vecWDVELESWP = [];
        vecWDVEMCSWP = [];
        vecWDVETESWP = [];
        vecWDVEAREA = [];
        matWDVENORM = [];
        matWVLST = [];
        matWDVE = [];
        valWNELE = 0;
        matWCENTER = [];
        matWCOEFF = [];
        vecWK = [];
        matWADJE = [];
        vecWDVEPANEL = [];
        valLENWADJE = 0;
        vecWKGAM = [];
        vecWDVESYM = [];
        vecWDVETIP = [];
        vecWDVEWING = [];
        vecCTCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCFyCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCFxCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCPCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCQCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCMxCONV = nan(round(max(vecAZNUM)),valNUMRO);
        vecCMyCONV = nan(round(max(vecAZNUM)),valNUMRO);
        matDISNORM = [];
        matDISTHRUST = [];
        matDISAXIAL = [];
        matDISSIDE = [];
        vecCDPDIST = [];
        vecCLPDIST = [];
        gamma_old = [];
        dGammadt = [];

        % Building rotor resultant
        [vecR] = fcnRESROTOR(valNELE, 0, matCENTER, matDVENORM, matUINF, ...
            valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, ...
            vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, ...
            vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE);

        % Solving for rotor coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

        for valTIMESTEP = 1:valMAXTIME
            %% Timestep to solution
            %   Move rotor
            %   Generate new wake elements
            %   Create and solve WD-Matrix for new elements
            %   Solve wing D-Matrix with wake-induced velocities
            %   Solve entire WD-Matrix
            %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
            %   Calculate surface normal forces
            %   Calculate DVE normal forces
            %   Calculate induced drag
            %   Calculate force and moment coefficents
            %   Calculate viscous effects

            %% Moving the rotor  
            [matVLST, matCENTER, matROTAX, matNEWWAKE, matNPNEWWAKE, ...
                vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, ...
                vecDVEYAW, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
                vecDVEAREA, matDVENORM] = fcnMOVEMULTIROTOR(vecRPM, ...
                matROTAX, matVLST, matCENTER, matNPVLST, vecDVETE, matDVE, ...
                vecAZNUM, valJ, valDIA, valALPHAR, vecDVEROTOR, ...
                vecDVEVLSTROTOR, vecRODIR);

           % Generate new wake elements
            [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, ...
                vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM,...
                matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, ...
                matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM,...
                vecWDVETIP, vecWKGAM, vecWDVEWING] ...
                = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, ...
                matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, ...
                vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVEMCSWP, ...
                vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, ...
                valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, ...
                matWADJE, matNPVLST, vecDVEPANEL, vecWDVEPANEL, vecSYM, ...
                valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, ...
                vecDVEWING, vecWDVEWING, flagSTEADY, valWSIZE);       

            %% Creating and solving WD-Matrix for latest row of wake elements
            % We need to grab from matWADJE only the values we need for this latest row of wake DVEs
            idx = sparse(sum(ismember(matWADJE,[((valWNELE - valWSIZE) + 1)...
                :valWNELE]'),2)>0 & (matWADJE(:,2) == 4 | matWADJE(:,2) == 2));
             temp_WADJE = [matWADJE(idx,1) - (valTIMESTEP-1)*valWSIZE ...
                 matWADJE(idx,2) matWADJE(idx,3) - (valTIMESTEP-1)*valWSIZE];

              [matWD, vecWR] = fcnWDWAKE([1:valWSIZE]', temp_WADJE, ...
                  vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVESYM(end-valWSIZE+ ...
                  1:end), vecWDVETIP(end-valWSIZE+1:end), vecWKGAM(end- ...
                  valWSIZE+1:end)); 
             [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, vecWR, ...
                 valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end- ...
                 valWSIZE+1:end)); 

            %% Generate rotor resultant
            % Calculate inflow velocity   
            [matUINF, matUINFTE, matTEPTS, vecTHETA] = fcnUINFMULTIROT( ...
                matCENTER, matROTAX, valTIMESTEP, vecRPM, valALPHAR, ...
                vecAZNUM, valDIA, valJ, valNUMB, vecDVEROTOR, vecDVEHVSPN, ...
                vecDVETE, matVLST, matDVE, vecRODIR);

            % Generate rotor resultant
            [vecR] = fcnRESROTOR(valNELE, valTIMESTEP, matCENTER, ...
                matDVENORM, matUINF, valWNELE, matWDVE, matWVLST, matWCOEFF,...
                vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH,...
                vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, ...
                flagSTEADY);

            [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

            % Creating and solving WD-Matrix
            [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, ...
                vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
            [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, ...
                vecWDVEHVSPN);

            %% Relaxing wake
            if valTIMESTEP > 2 && flagRELAX == 1
                % Apply wake relaxation   
                [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, ...
                    vecWDVEYAW, vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, ...
                    vecWDVEAREA, matWCENTER, matWDVENORM, matWVLST, matWDVE,...
                    matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] =  ...
                    fcnRELAXROTORWAKE(matUINF, matCOEFF, matDVE, matVLST, ...
                    matWADJE, matWCOEFF, matWDVE, matWVLST, valDELTIME, ...
                    valNELE, valTIMESTEP, valWNELE, valWSIZE, ...
                    vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEPITCH, ...
                    vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, ...
                    vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
                    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, ...
                    vecWDVEYAW, vecWK, vecWDVEWING, flagSTEADY, ...
                    valAZNUM, flagHOVERWAKE, vecN, vecM, vecDVELE, ...
                    vecDVEPANEL);   
                % Creating and solving WD-Matrix
                [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, ...
                    vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, ...
                    vecWKGAM, vecWDVEHVSPN);
            end
            %% Calculate Forces
            [vecCLPDIST, vecCDPDIST, thrustind, thrustfree, thrustCFfree, ...
                tempTi, Pthrust, difthrustP, Ptorque, valCT, valFy, valFx, valCQ, ...
                valCP, valCMy, valCMx, vecCTCONV, vecCFyCONV, vecCFxCONV, ...
                vecCQCONV, vecCPCONV, vecCMyCONV, vecCMxCONV, vecDISTHRUST, ...
                vecDISNORM, vecDISAXIAL, vecDISSIDE, matDISNORM, ...
                matDISTHRUST, matDISAXIAL, matDISSIDE, matWUINF, gamma_old,...
                dGammadt] = fcnRFORCES(flagVERBOSE, valKINV, valMAXTIME, ...
                vecAZNUM, valDIA, vecRPM, valWSIZE, valTIMESTEP, valNELE, ...
                valWNELE, seqALPHAR, vecDVEAREA, vecAIRFOIL, vecN, vecM, ...
                vecDVEPANEL, vecQARM,vecDVEPITCH, vecDVETE, vecDVEWING, ...
                vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, ...
                vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, ...
                vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, ...
                vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, vecCFxCONV, ...
                vecCFyCONV, vecCPCONV, vecCMxCONV, vecCMyCONV, vecCQCONV, ...
                matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, ...
                matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE, ...
                matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE, ...
                flagSTEADY, gamma_old, dGammadt, flagVISCOUS, valNUMRO, vecDVEROTOR);

            %% Parfor Data saving
             OUTPUT.CT(valTIMESTEP,:,ji,ai) = mean(vecCTCONV,'omitnan');
             OUTPUT.CP(valTIMESTEP,:,ji,ai) = mean(vecCPCONV,'omitnan');
             OUTPUT.CFx(valTIMESTEP,:,ji,ai) = mean(vecCFxCONV,'omitnan');
             OUTPUT.CFy(valTIMESTEP,:,ji,ai) = mean(vecCFyCONV,'omitnan');
             OUTPUT.CQ(valTIMESTEP,:,ji,ai) = mean(vecCQCONV,'omitnan');
             OUTPUT.CMx(valTIMESTEP,:,ji,ai) = mean(vecCMxCONV,'omitnan');
             OUTPUT.CMy(valTIMESTEP,:,ji,ai) = mean(vecCMyCONV,'omitnan');

%              matthrustind(:,valTIMESTEP,ji,ai) = thrustind;
%              matthrustfree(:,valTIMESTEP,ji,ai) = thrustfree;
%              matthrustCFfree(:,valTIMESTEP,ji,ai) = thrustCFfree;
%              mattempTi(:,valTIMESTEP,ji,ai) = tempTi;
%              matPthrust(:,valTIMESTEP,ji,ai) = Pthrust;
%              matPtorque(:,valTIMESTEP,ji,ai) = Ptorque;
%              matdifthrustP(:,valTIMESTEP,ji,ai) = difthrustP;

            temp = valTIMESTEP - ceil((floor((valTIMESTEP-1)./((vecAZNUM)))).*((vecAZNUM)));
            for i = 1:valNUMRO
                 OUTPUT.matCTTIMESTEP(valTIMESTEP,i,ji,ai) = vecCTCONV(temp(i),i);
                 OUTPUT.matCFyTIMESTEP(valTIMESTEP,i,ji,ai) = vecCFyCONV(temp(i),i);
                 OUTPUT.matCFxTIMESTEP(valTIMESTEP,i,ji,ai) = vecCFxCONV(temp(i),i);
                 OUTPUT.matCPTIMESTEP(valTIMESTEP,i,ji,ai) = vecCPCONV(temp(i),i);
                 OUTPUT.matCQTIMESTEP(valTIMESTEP,i,ji,ai) = vecCQCONV(temp(i),i);
                 OUTPUT.matCMxTIMESTEP(valTIMESTEP,i,ji,ai) = vecCMxCONV(temp(i),i);
                 OUTPUT.matCMyTIMESTEP(valTIMESTEP,i,ji,ai) = vecCMyCONV(temp(i),i);
            end
            if flagSAVE ==1
                save(filename)
            end

            if flagPRINT == 1   
                fprintf('    %.0f    ',valTIMESTEP);
                for i = 1:valNUMRO
                    fprintf('CT = %0.3f CP = %0.3f   ',mean(vecCTCONV(:,i),1,'omitnan'), mean(vecCPCONV(:,i),1,'omitnan'));
                end
                fprintf('\n')
            end

            if flagPROGRESS == 1
                progressbar([],[],valTIMESTEP/valMAXTIME)
            end
        end
        fprintf('Completed Advance Ratio: %.1f\n\n',seqJ(ji))
        if flagPROGRESS == 1
            progressbar([],ji/length(seqJ))
        end
    end   
	if flagPROGRESS == 1
        progressbar(ai/length(seqALPHAR),[]);
	end
end

if flagSAVE ==1
    save(filename)
end
fprintf('\nDONE\n');
runtime = toc;
if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
end