clc
clear

warning off

disp('===========================================================================');
disp('+---------------+');
disp('| RYERSON       |       VAPTOR (Based on FreeWake 2015)');
disp('| APPLIED       |       Running Version 2017.2');
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
%% Reading in geometry
%strFILE = 'inputs/ROTORINPUT_MA11by7.txt';
%strFILE = 'inputs/rectangle.txt';
strFILE = 'inputs/TMotor.txt';
%strFILE = 'inputs/standard_cirrus.txt';
%strFILE = 'inputs/TMotorQuad.txt';

[flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, valDELTAE, ...
    seqALPHAR, seqJ, valRPM, valDENSITY, valKINV, valAREA, valDIA, ...
    valNUMB ,vecROTAX0, valPANELS, matGEOM, vecAIRFOIL, vecN, vecM, ...
    vecSYM, valINTERF] = fcnVAPTORREAD(strFILE);

% [flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, ...
%     valDELTAE, seqALPHAR, valJ, valRPM, valDENSITY, valKINV, valAREA, valDIA,...
%     valNUMB, valNUMRO, matROTAX, vecRODIR, valPANELS, vecROTAXLOC, matGEOM, vecAIRFOIL, vecN, vecM, vecSYM, ...
%     valINTERF] = fcnVAPTORREADMULTI(strFILE);

flagPRINT   = 0;
flagPLOT    = 1;
flagPLOTWAKEVEL = 0;
flagVERBOSE = 0;
flagSAVE = 0;
flagPROGRESS = 0;
filename = 'ZeroAlphaSweepJune29'; % Save workspace name

% %% Apply multiple rotors
% [valPANELS, vecN, vecM, vecAIRFOIL, matGEOM ] = fcnGEOMMULTIROTOR(valNUMRO, valPANELS, vecRODIR, vecROTAXLOC, vecN, vecM, vecAIRFOIL, matGEOM, matROTAX);

%% Discretize geometry into DVEs
[matCENTER0, vecDVEHVSPN0, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, ...
    vecDVETESWP, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, ...
    matDVENORM0, matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, ...
    vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

%% Apply multiple blades
[valNELE, matNPVLST0, vecAIRFOIL, vecDVELE, vecDVETE, vecDVEYAW, ...
    vecDVEPANEL, vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, ...
    vecDVEROLL, vecDVEAREA, vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, ...
    vecDVELESWP, vecDVEHVCRD, vecDVEHVSPN0, vecSYM, vecQARM, matADJE, matCENTER0,...
    matVLST0, matDVE, matDVENORM0] = fcnDVEMULTIROTOR(valNELE, valNUMB, ...
    vecDVETIP, vecDVETESWP, vecDVEPITCH, vecDVEWING, vecDVEMCSWP, vecM, ...
    vecN, vecDVEPANEL, vecDVEROLL, vecDVELESWP, vecDVEYAW, vecDVEHVCRD, ...
    vecDVEHVSPN0, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, ...
    vecROTAX0, vecAIRFOIL, matNPVLST0, matDVE, matADJE, matVLST0, ...
    matCENTER0, matDVENORM0);

valWSIZE = length(nonzeros(vecDVETE)); % Amount of wake DVEs shed each timestep
%% Add boundary conditions to D-Matrix
[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN0, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN0);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, ...
    matDVENORM0, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, ...
    vecDVETESWP, vecDVEHVSPN0, vecDVEHVCRD,vecSYM);

%% Performance sweeps
% Preallocating for a turbo-boost in performance
matCTCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCFyCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCFxCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCPCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCQCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCMxCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matCMyCONV = zeros(valMAXTIME,length(seqJ),length(seqALPHAR));
matSWPDISNORM = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISTHRUST = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISAXIAL = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISSIDE = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
convCT = zeros(valMAXTIME,1);
convCP = zeros(valMAXTIME,1);
matthrustind = zeros(valNELE,valMAXTIME,length(seqJ));
matthrustfree = zeros(valNELE,valMAXTIME,length(seqJ));
matthrustCFfree = zeros(valNELE,valMAXTIME,length(seqJ));
mattempTi = zeros(valNELE,valMAXTIME,length(seqJ));
matPthrust = zeros(valNELE,valMAXTIME,length(seqJ));
matdifthrustP = zeros(valNELE,valMAXTIME,length(seqJ));
if flagPROGRESS == 1 % Apply progress bar
    progressbar('Angle of Attack','Advance Ratio')
end

for ai = 1:length(seqALPHAR)
    if flagPROGRESS == 1
        progressbar([],0)
    end
    for ji= 1:length(seqJ)
    if flagPROGRESS == 1
        progressbar([],ji/length(seqJ))
    end
    valALPHAR = deg2rad(seqALPHAR(ai));
    valJ = seqJ(ji);
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNPVLST = matNPVLST0;
    vecROTAX = vecROTAX0;
    vecDVEHVSPN = vecDVEHVSPN0;
    matDVENORM = matDVENORM0;
    
 	fprintf('      ANGLE OF ATTACK = %0.3f DEG\n',seqALPHAR(ai));
    fprintf('      ADVANCE RATIO = %0.3f\n',valJ);
    fprintf('      RPM = %0.3f\n\n',valRPM);
    fprintf('-----------------------------------------\n');
    fprintf('  Timestep    Thurst          Power\n');
    fprintf('-----------------------------------------\n');
    
    % Calculate inflow velocity at each control point
    [matUINF, matUINFTE, matTEPTS, vecTHETA] = ...
        fcnUINFROT(matCENTER, vecROTAX, 0, valRPM, valALPHAR, ...
        valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, vecDVETE, matVLST, ...
        matDVE);
    
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
    vecCTCONV = [];
    vecCFyCONV = [];
    vecCFxCONV = [];
    vecCPCONV = [];
    vecCQCONV = [];
    vecCMxCONV = [];
    vecCMyCONV = [];
    matDISNORM = [];
    matDISTHRUST = [];
    matDISAXIAL = [];
    matDISSIDE = [];
    
    % Building rotor resultant
    [vecR] = fcnRESROTOR(valNELE, 0, matCENTER, matDVENORM, matUINF, ...
        valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, ...
        vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, ...
        vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE);

    % Solving for rotor coefficients
    [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
    for valTIMESTEP = 1:valMAXTIME
%         if valTIMESTEP <= 20
%             valJ = 1;
%             flagRELAX = 0;
%         else
%             valJ = 0;
%             flagRELAX = 1;
%         end
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
        [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE, ... 
            vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, ...
            vecDVEYAW, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
            vecDVEAREA, matDVENORM] = fcnMOVEROTOR(vecROTAX, matVLST, ...
            matCENTER, matNPVLST, vecDVETE, matDVE, valAZNUM, valJ, ...
            valDIA, valALPHAR);
        
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
        [matUINF, matUINFTE, matTEPTS, vecTHETA] = ...
            fcnUINFROT(matCENTER, vecROTAX, valTIMESTEP, valRPM, ...
            valALPHAR, valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, ...
            vecDVETE, matVLST, matDVE);
    
        % Generate rotor resultant
        [vecR] = fcnRESROTOR(valNELE, valTIMESTEP, matCENTER, ...
            matDVENORM, matUINF, valWNELE, matWDVE, matWVLST, matWCOEFF,...
            vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH,...
            vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE);
            
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
                matWADJE, matWCOEFF, matWDVE, matWVLST, valAZNUM, ...
                valRPM, valNELE, valTIMESTEP, valWNELE, valWSIZE, ...
                vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEPITCH, ...
                vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, ...
                vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
                vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, ...
                vecWDVEYAW, vecWK, vecWDVEWING);   
            % Creating and solving WD-Matrix
            [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, ...
                vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
            [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, ...
                vecWKGAM, vecWDVEHVSPN);
        end
%         [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
%         [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
%         [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);

        %% Calculate Forces
        [thrustind, thrustfree, thrustCFfree, tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, valCP, valCMy, valCMx, vecCTCONV, ...
            vecCFyCONV, vecCFxCONV, vecCQCONV, vecCPCONV, vecCMyCONV, ...
            vecCMxCONV, vecDISTHRUST, vecDISNORM, vecDISAXIAL, ...
            vecDISSIDE, matDISNORM, matDISTHRUST, matDISAXIAL, ...
            matDISSIDE, matWUINF] = fcnRFORCES(flagVERBOSE,valKINV, valMAXTIME, valAZNUM, valDIA, valRPM, valWSIZE,...
            valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEAREA, vecAIRFOIL,  vecN, vecM, vecDVEPANEL, vecQARM,...
            vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK,...
            vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, ...
            vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL,...
            vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, ...
            vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, ...
            vecCTCONV, vecCFxCONV, vecCFyCONV, vecCPCONV, vecCMxCONV, ...
            vecCMyCONV, vecCQCONV, matVLST, matDVE, matUINF, matCOEFF, ...
            matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS,...
            matUINFTE, matDISNORM, matDISTHRUST, matDISAXIAL, matDISSIDE);
        

        temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
        
        %% Parfor Data saving
         %convCT(valTIMESTEP) = mean(vecCTCONV);
         %convCP(valTIMESTEP) = mean(vecCPCONV);
         matthrustind(:,valTIMESTEP,ji) = thrustind;
         matthrustfree(:,valTIMESTEP,ji) = thrustfree;
         matthrustCFfree(:,valTIMESTEP,ji) = thrustCFfree;
         mattempTi(:,valTIMESTEP,ji) = tempTi;
         matPthrust(:,valTIMESTEP,ji) = Pthrust;
         matdifthrustP(:,valTIMESTEP,ji) = difthrustP;
         matCTCONV(valTIMESTEP,ji,ai) = vecCTCONV(temp);
         matCFyCONV(valTIMESTEP,ji,ai) = vecCFyCONV(temp);
         matCFxCONV(valTIMESTEP,ji,ai) = vecCFxCONV(temp);
         matCPCONV(valTIMESTEP,ji,ai) = vecCPCONV(temp);
         matCQCONV(valTIMESTEP,ji,ai) = vecCQCONV(temp);
         matCMxCONV(valTIMESTEP,ji,ai) = vecCMxCONV(temp);
         matCMyCONV(valTIMESTEP,ji,ai) = vecCMyCONV(temp);
%         if valTIMESTEP == valMAXTIME
%             matSWPDISNORM(:,:,ji,ai) = matDISNORM;
%             matSWPDISTHRUST(:,:,ji,ai) = matDISTHRUST;
%             matSWPDISAXIAL(:,:,ji,ai) = matDISAXIAL;
%             matSWPDISSIDE(:,:,ji,ai) = matDISSIDE;
%          end
           fprintf('      %.0f      CT = %0.3f      CP = %0.3f\n',valTIMESTEP, mean(vecCTCONV),  mean(vecCPCONV));       
    end
    fprintf('Completed Advance Ratio: %.1f\n\n',seqJ(ji))    
    end
    if flagPROGRESS == 1
        progressbar(ai/length(seqALPHAR))
    end
end
% if flagSAVE ==1
%     save(filename)
% end
fprintf('\nDONE\n');
% save(filename)
if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
end