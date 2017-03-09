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

[flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, valDELTAE, ...
    seqALPHAR, seqJ, valRPM, valDENSITY, valKINV, valAREA, valDIA, ...
    valNUMB ,vecROTAX0, valPANELS, matGEOM, vecAIRFOIL, vecN, vecM, ....
    vecSYM, valINTERF] = fcnVAPTORREAD(strFILE);

flagPRINT   = 1;
flagPLOT    = 1;
flagPLOTWAKEVEL = 0;
flagVERBOSE = 0;
flagSAVE = 1;
filename = 'TMotor Data'; % Save workspace name


%% Discretize geometry into DVEs
[matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, ...
    vecDVETESWP, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, ...
    matDVENORM0, matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, ...
    vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);


%% Apply multiple blades
if valNUMB > 1
[valNELE, matNPVLST0, vecAIRFOIL, vecDVELE, vecDVETE, vecDVEYAW, vecDVEPANEL, ...
    vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, vecDVEROLL, vecDVEAREA,...
    vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, vecDVELESWP, vecDVEHVCRD, ...
    vecDVEHVSPN0, vecSYM, matADJE, matCENTER0, matVLST0, matDVE, ...
    matDVENORM0] = fcnDVEMULTIROTOR(valNELE, valNUMB, vecDVETIP, ...
    vecDVETESWP, vecDVEPITCH, vecDVEWING, vecDVEMCSWP, vecM, vecN, ...
    vecDVEPANEL, vecDVEROLL, vecDVELESWP, vecDVEYAW, vecDVEHVCRD, ...
    vecDVEHVSPN, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, ...
    vecROTAX0, vecAIRFOIL, matNPVLST0, matDVE, matADJE, matVLST0, ...
    matCENTER0, matDVENORM0);
end

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
matCTCONV = zeros(valMAXTIME, length(seqJ),length(seqALPHAR));
matCQ = zeros(valMAXTIME, length(seqJ),length(seqALPHAR));
progress = zeros(length(seqJ),length(seqALPHAR));

for ai = 1:length(seqALPHAR)
    for ji= 1:length(seqJ)
    
    valALPHAR = deg2rad(seqALPHAR(ai));
    valJ = seqJ(ji);
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNPVLST = matNPVLST0;
    vecROTAX = vecROTAX0;
    vecDVEHVSPN = vecDVEHVSPN0;
    matDVENORM = matDVENORM0;
    
% 	fprintf('      ANGLE OF ATTACK = %0.3f DEG\n',seqALPHAR(ai));
%     fprintf('      ADVANCE RATIO = %0.3f\n',valJ);
%     fprintf('      RPM = %0.3f\n\',valRPM);
%     fprintf('\n');
    
    % Calculate inflow velocity at each control point
    [matUINF, matUINFTE, matTEPTS, vecTHETA, vecCPRADI] = fcnUINFROT(matCENTER, vecROTAX, 0, valRPM, valALPHAR, ...
        valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, vecDVETE, matVLST, matDVE);
    
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
        [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE, vecDVEHVSPN, ...
            vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
            vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM] = ...
            fcnMOVEROTOR(vecROTAX, matVLST, matCENTER, matNPVLST, ...
            vecDVETE, matDVE, valAZNUM, valJ, valDIA, valALPHAR);
        
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
        [matUINF, matUINFTE, matTEPTS, vecTHETA,vecCPRADI] = fcnUINFROT(matCENTER, vecROTAX, valTIMESTEP, ...
        valRPM, valALPHAR, valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, ...
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

        % Calculate Forces
        [valCT, valCQ, vecCTCONV] = fcnRFORCES(valAZNUM, valDIA, valRPM, valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecCPRADI,vecDVEPITCH, vecDVETE, vecDVEWING, vecWDVEWING, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, vecCTCONV, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF, matTEPTS, matUINFTE);
        fprintf('      CT = %0.3f\n',mean(vecCTCONV));
        matCTCONV(valTIMESTEP,ji,ai) = mean(vecCTCONV);
        matCQ(valTIMESTEP,ji,ai) = valCQ;
        
    end
fprintf('Completed Advance Ratio: %.1f\n\n',seqJ(ji))
    end
end

if flagSAVE ==1
    save(filename)
end

fprintf('DONE\n');
if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
end