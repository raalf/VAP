clc
clear

warning off

tic
disp('===========================================================================');
disp('+---------------+');
disp('| RYERSON       |       WinDySIM (Based on FreeWake 2015)');
disp('| APPLIED       |       Running Version 2017.03');
disp('| AERODYNAMICS  |       Includes stall model');
disp('| LABORATORY OF |       No trim solution');
disp('| FLIGHT        |       o                         o');
disp('+---------------+        \                       /');
disp('                          \                     /');
disp('                           \                   /');
disp('                            \       -^-       /');
disp('                             \    _([|])_    /');
disp('                             _\__/ \   / \__/_');
disp('   +X+````````\\\\\\RCAF\\\\\___/\  \ /  /\___/////RCAF//////''''''''+X+');
disp('              /             ||  \ \_(o)_/ /  ||             \');
disp('            +X+     "```````\\__//\__^__/\\__//```````"     +X+');
disp('                              |H|   |H|   |H|');
disp('   ___________________________/______Y______\___________________________');
disp('                            {}+      |      +{}');
disp('                                   {}+{}');
disp('===========================================================================');
disp(' ');

%% Best Practices
% 1. Define wing from one wingtip to another in one direction
% 2. When using symmetry, define from symmetry plane outward

%% Reading in geometry

% strFILE = 'inputs/VAP christmas.txt';
strFILE = 'inputs/VAP_HALE.txt';
strSTRUCT_INPUT = 'inputs/Struct_Input_HALE.txt';

[flagRELAX, flagSTEADY, flagSTIFFWING, valAREA, valSPAN, valCMAC, valWEIGHT, ...
    valCM, seqALPHA, seqBETA, valKINV, valUINF, valDENSITY, valPANELS, matGEOM, vecSYM, ...
    vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
    valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
    valINTERF] = fcnVAPREAD(strFILE);

<<<<<<< HEAD:WinDySIM_MAIN_STATIC.m
[valSDELTIME, valSTIFFSTEPS, flagSTATIC, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF] = fcnSTRUCTREAD(strSTRUCT_INPUT);

% seqALPHA = [2];
=======
flagSTEADY = 2;
>>>>>>> refs/remotes/origin/master:VAP_MAIN.m

% strFILE = 'inputs/input.txt';
% strFILE = 'inputs/Config 1.txt';
% strFILE = 'inputs/Config 2.txt';

% [flagRELAX, flagSTEADY, valAREA, valSPAN, valCMAC, valWEIGHT, ...
%     seqALPHA, seqBETA, valKINV, valDENSITY, valPANELS, matGEOM, vecSYM, ...
%     vecAIRFOIL, vecN, vecM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, ...
%     valFTURB, valFPWIDTH, valDELTAE, valDELTIME, valMAXTIME, valMINTIME, ...
%     valINTERF] = fcnFWREAD(strFILE);

% flagRELAX = 0;
% valMAXTIME = 68;

flagPRINT   = 1;
flagPLOT    = 1;
flagPLOTWAKEVEL = 0;
flagVERBOSE = 0;

%% Discretize geometry into DVEs

[matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST0, matNPVLST0, matNTVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVEWING, vecDVELE, vecDVETE, vecDVEPANEL, vecLEDVES] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

valWSIZE = length(nonzeros(vecDVETE)); % Amount of wake DVEs shed each timestep

matNPDVE = matDVE;

%% Discretize geometry into structural parameters
% [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
%     vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matVLST0, matDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW);

%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVEWING, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecDVEHVCRD,vecSYM);

%% Alpha Loop

% Preallocating for a turbo-boost in performance
vecCL = zeros(valMAXTIME, length(seqALPHA));
vecCDI = zeros(valMAXTIME, length(seqALPHA));
vecE = zeros(valMAXTIME, length(seqALPHA));

for ai = 1:length(seqALPHA)
    
    valALPHA = deg2rad(seqALPHA(ai));
    
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNTVLST = matNTVLST0;
    matNPVLST = matNPVLST0;
    
    for bi = 1:length(seqBETA)
        
        fprintf('      ANGLE OF ATTACK = %0.3f DEG\n',seqALPHA(ai));
        fprintf('    ANGLE OF SIDESLIP = %0.3f DEG\n',seqBETA(bi));
        fprintf('\n');
        
        valBETA = deg2rad(seqBETA(bi));
        
        % Determining freestream vector
        vecUINF = fcnUINFWING(valALPHA, valBETA, valUINF);
        
        matUINF = repmat(vecUINF,size(matCENTER,1),1);
        
        % Initializing wake parameters
        matWAKEGEOM = [];
        matNPWAKEGEOM = [];
        matDEFGLOB = [];
        matTWISTGLOB = [];
        matDEF = zeros(2,valNELE+1); % <-------------- FIX THIS TO WORK WITH m > 1
        matTWIST = zeros(2,valNELE+1); % <-------------- FIX THIS TO WORK WITH m > 1
        matSLOPE = [];
        vecLIFTSTATIC = [];
        vecMOMSTATIC = [];
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
        gamma_old = [];
        dGammadt = [];
        
        % Building wing resultant
        [vecR] = fcnRWING(valNELE, 0, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
            matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            vecWDVETESWP, vecSYM, valWSIZE, flagSTEADY);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
        [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
            vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
        
        [ledves, ~, ~] = find(vecDVELE > 0);
        
        tempSPANDIST = matCENTER(ledves,2); % Y coordinate of DVE mid-point (point where DVEPITCH is applied) --> used as "x" term for linear interpolation

        vecEDGEPITCH = vecDVEPITCH(ledves); % DVE pitch along span --> used as "y" term for linear interpolation

        % Some setup of work to be able to perform linear interpolation without a
        % for loop
        tempSPANDIST = repmat(tempSPANDIST', size(tempSPANDIST,1),1);

        tempSPANDIST = triu(tempSPANDIST);

        vecEDGEPITCH = repmat(vecEDGEPITCH', size(vecEDGEPITCH,1),1);

        vecEDGEPITCH = triu(vecEDGEPITCH);

        vecEDGEPITCH = ((vecSPANDIST(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
            tempSPANDIST(1,1:(end-1)))).*(vecEDGEPITCH(2,2:end)-vecEDGEPITCH(1,1:(end-1))) + vecEDGEPITCH(1,1:(end-1)); % Linear interpolation

        % Adding in root and tip values using a linear extrapolation
        pitch_root = vecEDGEPITCH(1,1) - (tempSPANDIST(1,1) - vecSPANDIST(1)).*(vecEDGEPITCH(1,2)-vecEDGEPITCH(1,1))./(tempSPANDIST(1,2)-tempSPANDIST(1,1));

        pitch_tip = vecEDGEPITCH(1,end) + (vecSPANDIST(end) - tempSPANDIST(1,end)).*(vecEDGEPITCH(1,end)-vecEDGEPITCH(1,end-1))./(tempSPANDIST(1,end)-tempSPANDIST(1,end-1));

        vecEDGEPITCH = [pitch_root, vecEDGEPITCH, pitch_tip];
        
        n = 1;
        for valTIMESTEP = 1:valMAXTIME
            %% Timestep to solution
            %   Move wing
            %   Generate new wake elements
            %   Create and solve WD-Matrix for new elements
            %   Solve wing D-Matrix with wake-induced velocities
            %   Solve entire WD-Matrix
            %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
            %   Calculate surface normal forces
            %   Calculate DVE normal forces
            %   Calculate induced drag
            %   Calculate cn, cl, cy, cdi
            %   Calculate viscous effects
            
            %% Moving the wing and structure
                    
            if valTIMESTEP <= (n*valSTIFFSTEPS)
                [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matUINF, matDEFGLOB, matTWISTGLOB, valUINF] = fcnSTIFFWING_STATIC(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE,...
                    matNTVLST, matNPVLST, vecN, valTIMESTEP, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF, matNPDVE, matDEFGLOB, matTWISTGLOB);
                
                [matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST, vecSPANDIST, matSC, vecMAC] = fcnSTRUCTDIST(vecDVEHVSPN, vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF,...
                    vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL, vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW);
                              
            % Remaining timesteps compute wing deflection and translate the
            % wing accordingly

            else
                
                [valDELTIME, matEIx, matGJt, vecEA, vecCG, vecJT, vecLM, vecLSM, vecLSAC, matAEROCNTR, matSCLST,...
                    vecSPANDIST, matSC, vecMAC, vecDEF, vecTWIST, matDEFGLOB, matTWISTGLOB, matDEF, matTWIST, matSLOPE,...
                    matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
                    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, matVLST, matDVE, matCENTER, vecEDGEPITCH] = fcnFLEXWING_STATIC(vecDVEHVSPN,...
                    vecDVELE, vecDVETE, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF, matNPVLST, matNPDVE, vecDVEPANEL,...
                    vecN, vecM, vecDVEWING, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecLIFTDIST, vecMOMDIST, valSPAN, valTIMESTEP, matDEFGLOB, matTWISTGLOB,...
                    matSLOPE, vecLIFTSTATIC, vecMOMSTATIC, valALPHA, valBETA, matVLST, matCENTER, matDVE, vecCL, valWEIGHT, valAREA, valDENSITY, valUINF,...
                    flagSTATIC, valSDELTIME, valDELTIME, matDEF, matTWIST, valSTIFFSTEPS, vecCLDIST, vecEDGEPITCH, vecLEDVES);
                
                n = n + 1;

            end
            
            % Update structure location after moving wing
            [vecSPNWSECRD, vecSPNWSEAREA, matQTRCRD, vecQTRCRD] = fcnWINGSTRUCTGEOM(vecDVEWING, vecDVELE, vecDVEPANEL, vecM, vecN, vecDVEHVCRD, matDVE, matVLST, vecDVEAREA);
            
            %% Generating new wake elements
            [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
                = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, ...
                vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVEWING, vecWDVEWING, flagSTEADY, valWSIZE);
            
            %% Creating and solving WD-Matrix for latest row of wake elements
            % We need to grab from matWADJE only the values we need for this latest row of wake DVEs
            idx = sparse(sum(ismember(matWADJE,[((valWNELE - valWSIZE) + 1):valWNELE]'),2)>0 & (matWADJE(:,2) == 4 | matWADJE(:,2) == 2));
            temp_WADJE = [matWADJE(idx,1) - (valTIMESTEP-1)*valWSIZE matWADJE(idx,2) matWADJE(idx,3) - (valTIMESTEP-1)*valWSIZE];
            
            [matWD, vecWR] = fcnWDWAKE([1:valWSIZE]', temp_WADJE, vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVESYM(end-valWSIZE+1:end), vecWDVETIP(end-valWSIZE+1:end), vecWKGAM(end-valWSIZE+1:end));
            [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, vecWR, valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end-valWSIZE+1:end));
            
            %% Rebuilding and solving wing resultant
            [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
                matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVETESWP, vecSYM, valWSIZE, flagSTEADY);
            
            [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
            
            %% Creating and solving WD-Matrix
            [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
            [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
            
            %% Relaxing wake
            if valTIMESTEP > 2 && flagRELAX == 1
                
                [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
                    vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
                    matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKE(vecUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
                    matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
                    vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
                    vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING, flagSTEADY);
                
                % Creating and solving WD-Matrix
                [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
            end
            
            %% Timing
            %             eltime(valTIMESTEP) = toc;
            %             ttime(valTIMESTEP) = sum(eltime);
            
            %% Forces
            
            [vecCL(valTIMESTEP,ai), vecCLF(valTIMESTEP,ai),vecCLI(valTIMESTEP,ai),vecCDI(valTIMESTEP,ai), vecE(valTIMESTEP,ai), vecDVENFREE, vecDVENIND, ...
<<<<<<< HEAD:WinDySIM_MAIN_STATIC.m
                vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecLIFTDIST, vecMOMDIST, vecCLDIST] = ...
                fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP,...
                vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD,vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE,...
                valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, ...
                vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, vecSYM, vecDVETESWP, valAREA, valSPAN, valBETA, ...
                vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEAREA, vecSPNWSECRD, vecSPNWSEAREA, matQTRCRD, valDENSITY, valWEIGHT,...
                vecLEDVES, vecUINF, matSCLST, vecSPANDIST, matNPVLST, matNPDVE, matSC, vecMAC, valCM, valUINF, matAEROCNTR);
            
            if valTIMESTEP == valSTIFFSTEPS && flagSTATIC == 0
                
                vecLIFTSTATIC = vecLIFTDIST;
                vecMOMSTATIC = vecMOMDIST;
                
            end
=======
                vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, gamma_old, dGammadt] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, vecUINF, vecDVELESWP, ...
                vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE, valWNELE, matWDVE, matWVLST, ...
                matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, ...
                vecSYM, vecDVETESWP, valAREA, valSPAN, valBETA, vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, flagSTEADY, gamma_old, dGammadt, valDELTIME);
>>>>>>> refs/remotes/origin/master:VAP_MAIN.m
            
            if flagPRINT == 1 && valTIMESTEP == 1
                fprintf(' TIMESTEP    CL          CDI          Tip Def.       Twist (deg)\n'); %header
                fprintf('----------------------------------------------------------------\n'); 
            end
            if flagPRINT == 1 && flagSTIFFWING == 2
                fprintf('  %4d     %0.5f     %0.5f         %0.5f          %0.5f\n',valTIMESTEP,vecCL(valTIMESTEP,ai),vecCDI(valTIMESTEP,ai),...
                    matDEFGLOB(valTIMESTEP,end),(180/pi)*matTWISTGLOB(valTIMESTEP,end)); %valTIMESTEP
            else
                fprintf('  %4d     %0.5f     %0.5f\n',valTIMESTEP,vecCL(valTIMESTEP,ai),vecCDI(valTIMESTEP,ai)); %valTIMESTEP               
            end
            
%             fprintf('\n\tTimestep = %0.0f', valTIMESTEP);
%             fprintf('\tCL = %0.5f',vecCL(valTIMESTEP,ai));
%             fprintf('\tCDi = %0.5f',vecCDI(valTIMESTEP,ai));

        end
        
        %% Viscous wrapper
        
        [vecCLv(1,ai), vecCD(1,ai), vecPREQ(1,ai), valVINF(1,ai), valLD(1,ai)] = fcnVISCOUS(vecCL(end,ai), vecCDI(end,ai), ...
            valWEIGHT, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
            vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
            matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
            valFPWIDTH, valINTERF, vecDVEROLL);
                
    end
end

fprintf('\n');

%% Plotting

if flagPLOT == 1
    [hFig2] = fcnPLOTBODY(flagVERBOSE, valNELE, matDVE, matVLST, matCENTER);
    [hFig2] = fcnPLOTWAKE(flagVERBOSE, hFig2, valWNELE, matWDVE, matWVLST, matWCENTER);
    [hLogo] = fcnPLOTLOGO(0.97,0.03,14,'k','none');
    
    if flagPLOTWAKEVEL == 1
        try
        quiver3(matWDVEMP(:,1),matWDVEMP(:,2),matWDVEMP(:,3),matWDVEMPIND(:,1),matWDVEMPIND(:,2),matWDVEMPIND(:,3));
        end
    end

end

if flagSTIFFWING ~= 1
figure(3)
clf
plot(vecSPANDIST, matDEFGLOB(valTIMESTEP,(1:size(matDEFGLOB,2))));
ylabel('Deflection (m)')
xlabel('Span Location (m)')
hold on
yyaxis right
plot(vecSPANDIST, (180/pi)*matTWISTGLOB(valTIMESTEP,(1:size(matDEFGLOB,2))));
ylabel('Twist (deg)')
hold off

figure(4)
clf
plot(valDELTIME*(1:valTIMESTEP),(180/pi)*matTWISTGLOB(:,end))
xlabel('Time (s)')
ylabel('Tip Twist (deg)')
grid on
box on

figure(5)
clf
plot(valDELTIME*(1:valTIMESTEP),matDEFGLOB(:,end))
xlabel('Time (s)')
ylabel('Tip Deflection (m)')
grid on
box on
end

save('Structural_Dynamics.mat');

toc
%% Viscous wrapper

% whos