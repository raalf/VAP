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
%strFILE = 'inputs/TMotor.txt';
%strFILE = 'inputs/standard_cirrus.txt';
strFILE = 'inputs/TMotorQuad.txt';
%strFILE = 'inputs/TenzinHover.txt';
tic  
% [flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, valDELTAE, ...
%     seqALPHAR, seqJ, valRPM, valDENSITY, valKINV, valAREA, valDIA, ...
%     valNUMB ,vecROTAX0, valPANELS, matGEOM, vecAIRFOIL, vecN, vecM, ...
%     vecSYM, valINTERF] = fcnVAPTORREAD(strFILE);

[flagRELAX, flagSTEADY, valMAXTIME, valMINTIME, valAZNUM, ...
    valDELTAE, seqALPHAR, seqJ, vecRPM, valDENSITY, valKINV, valAREA, valDIA,...
    valNUMB, valNUMRO, matROTAX0, vecRODIR, valPANELS, vecROTAXLOC, matGEOM, vecAIRFOIL, vecN, vecM, vecSYM, ...
    valINTERF] = fcnVAPTORREADMULTI(strFILE);
%valRPM = vecRPM(1);
flagVISCOUS = 1;
flagPRINT   = 1;
flagPLOT    = 1;
flagPLOTWAKEVEL = 0;
flagVERBOSE = 0;
flagSAVE = 0;
flagPROGRESS = 0;
filename = 'TMotorSweepViscousApparentMass'; % Save workspace name

%% Apply multiple rotors
[valPANELS, vecN, vecM, vecAIRFOIL, matGEOM ] = fcnGEOMMULTIROTOR(valNUMRO, valPANELS, vecRODIR, vecROTAXLOC, vecN, vecM, vecAIRFOIL, matGEOM, matROTAX0);

%% Discretize geometry into DVEs
[~, vecDVEHVSPN0, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, ...
    vecDVETESWP, ~, ~, ~, vecDVEAREA, ...
    ~, matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, ...
    vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

%% Apply multiple blades
% [valNELE, matNPVLST0, vecAIRFOIL, vecDVELE, vecDVETE, vecDVEYAW, ...
%     vecDVEPANEL, vecDVETIP, vecDVEWING, vecDVESYM, vecM, vecN, ...
%     vecDVEROLL, vecDVEAREA, vecDVEPITCH, vecDVEMCSWP, vecDVETESWP, ...
%     vecDVELESWP, vecDVEHVCRD, vecDVEHVSPN0, vecSYM, vecQARM, matADJE, matCENTER0,...
%     matVLST0, matDVE, matDVENORM0] = fcnDVEMULTIROTOR(valNELE, valNUMB, ...
%     vecDVETIP, vecDVETESWP, vecDVEPITCH, vecDVEWING, vecDVEMCSWP, vecM, ...
%     vecN, vecDVEPANEL, vecDVEROLL, vecDVELESWP, vecDVEYAW, vecDVEHVCRD, ...
%     vecDVEHVSPN0, vecDVEAREA, vecDVESYM, vecDVELE, vecDVETE, vecSYM, ...
%     vecROTAX0, vecAIRFOIL, matNPVLST0, matDVE, matADJE, matVLST0, ...
%     matCENTER0, matDVENORM0);

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
matCTCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCFyCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCFxCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCPCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCQCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCMxCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matCMyCONV = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matSWPDISNORM = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISTHRUST = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISAXIAL = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
matSWPDISSIDE = zeros(valNELE,valAZNUM,length(seqJ),length(seqALPHAR));
convCT = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
convCP = zeros(valMAXTIME,valNUMRO,length(seqJ),length(seqALPHAR));
matthrustind = zeros(valNELE,valMAXTIME,length(seqJ));
matthrustfree = zeros(valNELE,valMAXTIME,length(seqJ));
matthrustCFfree = zeros(valNELE,valMAXTIME,length(seqJ));
mattempTi = zeros(valNELE,valMAXTIME,length(seqJ));
matPthrust = zeros(valNELE,valMAXTIME,length(seqJ));
matdifthrustP = zeros(valNELE,valMAXTIME,length(seqJ));
matINDUCED_VEL = zeros(valNELE,3,valAZNUM);
matCDPDIST = zeros(valNELE,length(seqJ));
matCLPDIST = zeros(valNELE,length(seqJ));
induced_vel1 = zeros(625,3,valAZNUM);
b = 0;


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
    
    if valJ == 10000
        flagHOVER = 1;
        hoverSCALING = 0.001;
        valJ = -(1)*valRPM/(valAZNUM)*(hoverSCALING) + (valMAXTIME-valAZNUM)*valRPM/(valAZNUM)*(hoverSCALING);
        valALPHAR = deg2rad(90);
    else
        flagHOVER = 0;
    end
    
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
%     [matUINF, matUINFTE, matTEPTS, vecTHETA] = ...
%         fcnUINFROT(matCENTER, vecROTAX, 0, valRPM, valALPHAR, ...
%         valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, vecDVETE, matVLST, ...
%         matDVE);
    
     [matUINF, matUINFTE, matTEPTS, vecTHETA] = fcnUINFMULTIROT(matCENTER, matROTAX, 0, vecRPM, valALPHAR, vecAZNUM, valDIA, valJ, valNUMB, vecDVEROTOR, vecDVEHVSPN, vecDVETE, matVLST, matDVE, vecRODIR);


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
    vecCTCONV = nan(max(vecAZNUM),valNUMRO);
    vecCFyCONV = nan(max(vecAZNUM),valNUMRO);
    vecCFxCONV = nan(max(vecAZNUM),valNUMRO);
    vecCPCONV = nan(max(vecAZNUM),valNUMRO);
    vecCQCONV = nan(max(vecAZNUM),valNUMRO);
    vecCMxCONV = nan(max(vecAZNUM),valNUMRO);
    vecCMyCONV = nan(max(vecAZNUM),valNUMRO);
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
            
        if flagHOVER == 1 && valJ ~= 0
            valJ = -(valTIMESTEP)*valRPM/(valAZNUM)*(hoverSCALING) + (valMAXTIME-valAZNUM)*valRPM/(valAZNUM)*(hoverSCALING);
        end
        
        %% Moving the rotor
%         [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE, ... 
%             vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, ...
%             vecDVEYAW, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
%             vecDVEAREA, matDVENORM] = fcnMOVEROTOR(vecROTAX, matVLST, ...
%             matCENTER, matNPVLST, vecDVETE, matDVE, valAZNUM, valJ, ...
%             valDIA, valALPHAR);
% [matVLST, matCENTER, vecROTAX, matNEWWAKE, matNPNEWWAKE, vecDVEHVSPN, ...
%     vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%     vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM] = ...
%     fcnMOVEROTOR(valDENSITY, valRPM, vecROTAX, vecCTCONV, matUINF, matVLST, matCENTER, matNPVLST, ...
%     vecDVETE, matDVE, valAZNUM, valJ, valDIA, valALPHAR);        

[matVLST, matCENTER, matROTAX, matNEWWAKE, matNPNEWWAKE, vecDVEHVSPN, ...
    vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM] = ...
    fcnMOVEMULTIROTOR(vecRPM, matROTAX, matVLST, matCENTER, matNPVLST, ...
    vecDVETE, matDVE, vecAZNUM, valJ, valDIA, valALPHAR, vecDVEROTOR, vecDVEVLSTROTOR, vecRODIR);


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
%         [matUINF, matUINFTE, matTEPTS, vecTHETA] = ...
%             fcnUINFROT(matCENTER, vecROTAX, valTIMESTEP, valRPM, ...
%             valALPHAR, valAZNUM, valDIA, valJ, valNUMB,vecDVEHVSPN, ...
%             vecDVETE, matVLST, matDVE);
%     

        [matUINF, matUINFTE, matTEPTS, vecTHETA] = fcnUINFMULTIROT(matCENTER, matROTAX, valTIMESTEP, vecRPM, valALPHAR, vecAZNUM, valDIA, valJ, valNUMB, vecDVEROTOR, vecDVEHVSPN, vecDVETE, matVLST, matDVE, vecRODIR);

        % Generate rotor resultant
        [vecR] = fcnRESROTOR(valNELE, valTIMESTEP, matCENTER, ...
            matDVENORM, matUINF, valWNELE, matWDVE, matWVLST, matWCOEFF,...
            vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH,...
            vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, flagSTEADY);
            
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
                vecWDVEYAW, vecWK, vecWDVEWING, flagSTEADY);   
            % Creating and solving WD-Matrix
            [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, ...
                vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
            [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, ...
                vecWKGAM, vecWDVEHVSPN);
        end
%         [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
%         [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
%         [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
%         %quiver3(matCENTER(:,1,1),matCENTER(:,2,1),matCENTER(:,3,1),matUINF(:,1,1),matUINF(:,2,1),matUINF(:,3,1))

        %% Calculate Forces
%HERE
        [vecCLPDIST, vecCDPDIST, thrustind, thrustfree, thrustCFfree, ...
            tempTi, Pthrust, difthrustP, valCT, valFy, valFx, valCQ, ...
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
        
        temp = valTIMESTEP - (floor((valTIMESTEP-1)/valAZNUM))*(valAZNUM);
        %% Parfor Data saving
         convCT(valTIMESTEP,:,ji,ai) = nanmean(vecCTCONV);
         convCP(valTIMESTEP,:,ji,ai) = nanmean(vecCPCONV);
         matthrustind(:,valTIMESTEP,ji) = thrustind;
         matthrustfree(:,valTIMESTEP,ji) = thrustfree;
         matthrustCFfree(:,valTIMESTEP,ji) = thrustCFfree;
         mattempTi(:,valTIMESTEP,ji) = tempTi;
         matPthrust(:,valTIMESTEP,ji) = Pthrust;
         matdifthrustP(:,valTIMESTEP,ji) = difthrustP;
         
           temp = valTIMESTEP - ceil((floor((valTIMESTEP-1)./((vecAZNUM)))).*((vecAZNUM)));
        for i = 1:valNUMRO
         matCTCONV(valTIMESTEP,i,ji,ai) = vecCTCONV(temp(i),i);
         matCFyCONV(valTIMESTEP,i,ji,ai) = vecCFyCONV(temp(i),i);
         matCFxCONV(valTIMESTEP,i,ji,ai) = vecCFxCONV(temp(i),i);
         matCPCONV(valTIMESTEP,i,ji,ai) = vecCPCONV(temp(i),i);
         matCQCONV(valTIMESTEP,i,ji,ai) = vecCQCONV(temp(i),i);
         matCMxCONV(valTIMESTEP,i,ji,ai) = vecCMxCONV(temp(i),i);
         matCMyCONV(valTIMESTEP,i,ji,ai) = vecCMyCONV(temp(i),i);
        end
%HERE
%         if valTIMESTEP == valMAXTIME
%             matSWPDISNORM(:,:,ji,ai) = matDISNORM;
%             matSWPDISTHRUST(:,:,ji,ai) = matDISTHRUST;
%             matSWPDISAXIAL(:,:,ji,ai) = matDISAXIAL;
%             matSWPDISSIDE(:,:,ji,ai) = matDISSIDE;
%          end
        if flagPRINT == 1   
            fprintf('    %.0f    ',valTIMESTEP);
            for i = 1:valNUMRO
                fprintf('CT = %0.3f CP = %0.3f   ',nanmean(vecCTCONV(:,i),1), nanmean(vecCPCONV(:,i),1));
            end
            fprintf('\n')
        end
%     if valTIMESTEP > valMAXTIME - valAZNUM
%         b = b+1;
%         x_offset = vecROTAX(1);
%         z_offset = vecROTAX(3);
%         x = -0.6:0.01:0.6;
%         x = x + x_offset;
%         fpg1 = zeros(size(x,2),3);
%         fpg1(:,1) = x';
%         fpg1(:,3) = (fpg1(:,3)+1)*z_offset;
% 
%         % Along y line when z=x=0
%         y = -0.6:0.01:0.6;
% 
%         fpg2 = zeros(size(y,2),3);
%         fpg2(:,2) = y';
%         fpg2(:,1) = (fpg2(:,1)+1)*x_offset;
%         fpg2(:,3) = (fpg2(:,3)+1)*z_offset;
% 
%         % Calculate induced velocities
%         induced_vel1(:,:,b) = fcnINDVEL(fpg1,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%         vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%         vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%         matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%         vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
% 
%         induced_vel2(:,:,b) = fcnINDVEL(fpg2,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%         vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%         vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%         matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%         vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
% %     end
%     if valTIMESTEP > valMAXTIME - valAZNUM
%         b = b+1;
%         x_offset = vecROTAX(1);
%         z_offset = vecROTAX(3);
% %         x = -0.6:0.01:0.6;
% %         x = x + x_offset;
%         [x,y] = meshgrid(-0.6:0.05:0.6, -0.6:0.05:0.6);
%         x = reshape(x,1,numel(x))' + x_offset;
%         y = reshape(y,1,numel(y))';
%         fpg = zeros(size(x,1),3);
%         fpg(:,1) = x;
%         fpg(:,2) = y;
%         fpg(:,3) = (fpg(:,3)+1)*z_offset;
% 
% %         % Along y line when z=x=0
% %         y = -0.6:0.01:0.6;
% % 
% %         fpg2 = zeros(size(y,2),3);
% %         fpg2(:,2) = y';
% %         fpg2(:,1) = (fpg2(:,1)+1)*x_offset;
% %         fpg2(:,3) = (fpg2(:,3)+1)*z_offset;
% 
%         % Calculate induced velocities
%         induced_vel1(:,:,b) = fcnINDVEL(fpg,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%         vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%         vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%         matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%         vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
% 
% %         induced_vel2(:,:,b) = fcnINDVEL(fpg2,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
% %         vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
% %         vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
% %         matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
% %         vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
%     end
   
    %fprintf('Completed Advance Ratio: %.1f\n\n',seqJ(ji))    
    end
    if flagPROGRESS == 1
        progressbar(ai/length(seqALPHAR));
    end
    end
end
if flagSAVE ==1
    save(filename)
end
fprintf('\nDONE\n');
runtime = toc;
% save(filename)
if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
end