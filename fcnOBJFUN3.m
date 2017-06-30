function out = fcnOBJFUN3(orig2)


%% Reading in geometry

filename = 'inputs/twoVehicles.vap';

[flagRELAX, flagSTEADY, flagTRI, matGEOM, valMAXTIME, valMINTIME, valDELTIME, valDELTAE, ...
    valDENSITY, valKINV, valVEHICLES, matVEHORIG, vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHROLL, ...
    vecVEHFPA, vecVEHTRK, ~, ~, vecAREA, vecSPAN, ~, ~, ...
    ~, vecSYM, vecN, vecM, ~, ~, ~, ~, ...
    vecSURFACEVEHICLE, valPANELS, ~, vecROTORRPM, vecROTDIAM, matROTORHUB, matROTORAXIS, vecROTORBLADES, ~, vecROTOR,...
    vecFTURB, vecFUSESECTIONS, matFGEOM, matSECTIONFUSELAGE, vecFUSEVEHICLE, matFUSEAXIS, matFUSEORIG...
    ] = fcnXMLREAD(filename);

matVEHORIG(2,:) = orig2;

% flagRELAX = 0
% valMAXTIME = 3

seqALPHA = 0;
seqBETA = 0;

flagPREVIEW = 0;

% tranlsate matGEOM to vehicle origin
matGEOM(:,1:3,:) = matGEOM(:,1:3,:)+permute(reshape(matVEHORIG(matGEOM(:,6,:),:)',3,2,[]),[2,1,3]);

[matCENTER0, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, vecDVEMCSWP, vecDVETESWP, ...
    vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVEAREA, matDVENORM, ...
    matVLST0, matNPVLST0, matDVE, valNELE, matADJE, ...
    vecDVESYM, vecDVETIP, vecDVESURFACE, vecDVELE, vecDVETE, vecDVEPANEL] = fcnGENERATEDVES(valPANELS, matGEOM, vecSYM, vecN, vecM);

% Identifying which DVEs belong to which vehicle, as well as which type of lifting surface they belong to (wing or rotor)
vecDVEVEHICLE = vecSURFACEVEHICLE(vecDVESURFACE);
vecDVEWING = vecDVESURFACE;

vecDVEROTOR = vecROTOR(vecDVEPANEL); % Alton-Y
vecDVEROTORBLADE = vecDVEROTOR; % Current rotor DVEs are for Blade 1 (they are duplicated to Blade 2, 3, etc etc below)
idx_rotor = vecDVEROTOR>0; % Alton-Y
vecDVEWING(idx_rotor) = 0;

matSURFACETYPE = zeros(size(unique(vecDVESURFACE),1),2);
matSURFACETYPE(nonzeros(unique(vecDVEWING)),1) = nonzeros(unique(vecDVEWING));
matSURFACETYPE(nonzeros(unique(vecDVESURFACE(idx_rotor))),2) = nonzeros(unique(vecDVEROTOR));

% Identifying which ROTOR belongs to which vehicle.
vecROTORVEH = vecSURFACEVEHICLE(matSURFACETYPE(:,2)~=0);

% Duplicate Blades in a Rotor
[ matVLST0, matNPVLST0, matCENTER0, matDVE, matADJE, vecDVEVEHICLE, ...
    vecDVEWING, vecDVEROTOR, matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, ...
    vecDVETIP, vecDVELE, vecDVETE, vecDVEROTORBLADE, vecDVESYM, ...
    valNELE ] = fcnDUPBLADE( vecROTORVEH, vecDVEROTOR, ...
    matVLST0, matCENTER0, matNPVLST0, matDVE, matADJE, vecROTORBLADES, ...
    valNELE, matROTORHUB, matVEHORIG, vecDVEVEHICLE, vecDVEWING, ...
    matSURFACETYPE, vecDVESURFACE, vecDVEPANEL, vecDVETIP, vecDVELE, ...
    vecDVETE, vecDVEROTORBLADE, vecDVESYM, matROTORAXIS );

matFUSEGEOM = fcnCREATEFUSE(matSECTIONFUSELAGE, vecFUSESECTIONS, matFGEOM, matFUSEAXIS, matFUSEORIG, vecFUSEVEHICLE);

[ matVEHUVW, matVEHROT, vecVEHPITCH, vecVEHYAW ] = fcnINITVEHICLE( vecVEHVINF, vecVEHALPHA, vecVEHBETA, vecVEHFPA, vecVEHROLL, vecVEHTRK );
[ matVLST0, matCENTER0, matFUSEGEOM, matROTORHUBGLOB] = fcnROTVEHICLE( matDVE, matVLST0, matCENTER0, valVEHICLES, vecDVEVEHICLE, matVEHORIG, matVEHROT, matFUSEGEOM, vecFUSEVEHICLE, matFUSEAXIS, matROTORHUB, matROTORAXIS, vecROTORVEH);

[ matUINF ] = fcnINITUINF( matCENTER0, matVEHUVW, matVEHROT, vecDVEVEHICLE, ...
    vecDVEROTOR, vecROTORVEH, matVEHORIG, matROTORHUBGLOB, matROTORAXIS, vecROTORRPM );

% update DVE params after vehicle rotation
[ vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW,...
    vecDVELESWP, vecDVEMCSWP, vecDVETESWP, vecDVEAREA, matDVENORM, ~, ~, ~ ] ...
    = fcnVLST2DVEPARAM(matDVE, matVLST0);

% [hFig2] = fcnPLOTBODY(0, valNELE, matDVE, matVLST0, matCENTER0,matFUSEGEOM);

%%
% valWSIZE = length(nonzeros(vecDVETE.*(vecDVEWING > 0))); % Amount of wake DVEs shed each timestep
valWSIZE = length(nonzeros(vecDVETE));
%% Add boundary conditions to D-Matrix

[matD] = fcnDWING(valNELE, matADJE, vecDVEHVSPN, vecDVESYM, vecDVETIP);

%% Add kinematic conditions to D-Matrix

[vecK] = fcnSINGFCT(valNELE, vecDVESURFACE, vecDVETIP, vecDVEHVSPN);
[matD] = fcnKINCON(matD, valNELE, matDVE, matCENTER0, matVLST0, matDVENORM, vecK, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELESWP, vecDVETESWP, vecDVEHVSPN, vecDVEHVCRD,vecSYM);

%% Alpha Loop

% Preallocating for a turbo-boost in performance
vecCL = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCLF = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCLI = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecCDI = nan(valMAXTIME,valVEHICLES,length(seqALPHA));
vecE = nan(valMAXTIME,valVEHICLES,length(seqALPHA));

for ai = 1:length(seqALPHA)
    
    valALPHA = deg2rad(seqALPHA(ai));
    
    % This is done for when we are using a parfor loop
    matCENTER = matCENTER0;
    matVLST = matVLST0;
    matNPVLST = matNPVLST0;
    
    for bi = 1:length(seqBETA)
                
        valBETA = deg2rad(seqBETA(bi));
        
        % Determining freestream vector
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
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
        
        % Building wing resultant
        [vecR] = fcnRWING(valNELE, 0, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
            matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            vecWDVETESWP, vecSYM, valWSIZE);
        
        % Solving for wing coefficients
        [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
        
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
            
            %% Moving the wing
            %             [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE.*(vecDVEWING > 0), matNPVLST);
            %             [matUINF, matVEHORIG, matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, vecDVETE.*(vecDVEWING > 0), matNPVLST, matFUSEGEOM, vecFUSEVEHICLE);
            [matUINF, matUINFTE, matVEHORIG, ...
                matVLST, matCENTER, ...
                matNEWWAKE, matNPNEWWAKE, ...
                matFUSEGEOM] = fcnMOVESURFACE(matVEHORIG, matVEHUVW, ...
                valDELTIME, matVLST, matCENTER, matDVE, vecDVEVEHICLE, ...
                vecDVETE, matNPVLST, matFUSEGEOM, vecFUSEVEHICLE, ...
                matVEHROT, vecROTORVEH, matROTORHUBGLOB, ...
                matROTORHUB, matROTORAXIS, vecDVEROTOR, vecROTORRPM );
            %% Generating new wake elements
            %             [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            %                 vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
            %                 = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
            %                 vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE.*(vecDVEWING > 0), matWADJE, matNPVLST, vecDVEPANEL, ...
            %                 vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVESURFACE, vecWDVEWING, flagSTEADY, valWSIZE);
            
            [matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matWADJE, matNPVLST, vecWDVEPANEL, valLENWADJE, vecWDVESYM, vecWDVETIP, vecWKGAM, vecWDVEWING] ...
                = fcnCREATEWAKEROW(matNEWWAKE, matNPNEWWAKE, matWAKEGEOM, matNPWAKEGEOM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                vecWDVEMCSWP, vecWDVETESWP, vecWDVEAREA, matWDVENORM, matWVLST, matWDVE, valWNELE, matWCENTER, matWCOEFF, vecWK, matCOEFF, vecDVETE, matWADJE, matNPVLST, vecDVEPANEL, ...
                vecWDVEPANEL, vecSYM, valLENWADJE, vecWKGAM, vecWDVESYM, vecWDVETIP, vecK, vecDVESURFACE, vecWDVEWING, flagSTEADY, valWSIZE);
            
            if flagPREVIEW ~= 1
                %% Creating and solving WD-Matrix for latest row of wake elements
                % We need to grab from matWADJE only the values we need for this latest row of wake DVEs
                idx = sparse(sum(ismember(matWADJE,[((valWNELE - valWSIZE) + 1):valWNELE]'),2)>0 & (matWADJE(:,2) == 4 | matWADJE(:,2) == 2));
                temp_WADJE = [matWADJE(idx,1) - (valTIMESTEP-1)*valWSIZE matWADJE(idx,2) matWADJE(idx,3) - (valTIMESTEP-1)*valWSIZE];
                
                [matWD, vecWR] = fcnWDWAKE([1:valWSIZE]', temp_WADJE, vecWDVEHVSPN(end-valWSIZE+1:end), vecWDVESYM(end-valWSIZE+1:end), vecWDVETIP(end-valWSIZE+1:end), vecWKGAM(end-valWSIZE+1:end));
                [matWCOEFF(end-valWSIZE+1:end,:)] = fcnSOLVEWD(matWD, vecWR, valWSIZE, vecWKGAM(end-valWSIZE+1:end), vecWDVEHVSPN(end-valWSIZE+1:end));
                
                %% Rebuilding and solving wing resultant
                [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, matUINF, valWNELE, matWDVE, ...
                    matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
                    vecWDVETESWP, vecSYM, valWSIZE);
                
                [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
                
                %% Creating and solving WD-Matrix
                [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
                
                %% Relaxing wake
                if valTIMESTEP > 2 && flagRELAX == 1
                    
                    [vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW,...
                        vecWDVELESWP, vecDVEWMCSWP, vecDVEWTESWP, vecWDVEAREA, matWCENTER, matWDVENORM, ...
                        matWVLST, matWDVE, matWDVEMP, matWDVEMPIND, idxWVLST, vecWK] = fcnRELAXWAKE(matUINF, matCOEFF, matDVE, matVLST, matWADJE, matWCOEFF, ...
                        matWDVE, matWVLST, valDELTIME, valNELE, valTIMESTEP, valWNELE, valWSIZE, vecDVETE, vecDVEHVSPN, vecDVEHVCRD, vecDVELESWP, ...
                        vecDVEPITCH, vecDVEROLL, vecDVETESWP, vecDVEYAW, vecK, vecSYM, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVELESWP, vecWDVEPITCH, ...
                        vecWDVEROLL, vecWDVESYM, vecWDVETESWP, vecWDVETIP, vecWDVEYAW, vecWK, vecWDVEWING);
                    
                    % Creating and solving WD-Matrix
                    [matWD, vecWR] = fcnWDWAKE([1:valWNELE]', matWADJE, vecWDVEHVSPN, vecWDVESYM, vecWDVETIP, vecWKGAM);
                    [matWCOEFF] = fcnSOLVEWD(matWD, vecWR, valWNELE, vecWKGAM, vecWDVEHVSPN);
                end
                
                %% Forces
                [vecCL(valTIMESTEP,:,ai), vecCLF(valTIMESTEP,:,ai), vecCLI(valTIMESTEP,:,ai), vecCDI(valTIMESTEP,:,ai), vecE(valTIMESTEP,:,ai), vecDVENFREE, vecDVENIND, ...
                    vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND] = fcnFORCES(matCOEFF, vecK, matDVE, valNELE, matCENTER, matVLST, matUINF, vecDVELESWP, ...
                    vecDVEMCSWP, vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, vecDVELE, vecDVETE, matADJE, valWNELE, matWDVE, matWVLST, ...
                    matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP, ...
                    vecSYM, vecDVETESWP, vecAREA, vecSPAN, [], vecDVEWING, vecWDVEWING, vecN, vecM, vecDVEPANEL, vecDVEVEHICLE, valVEHICLES);
                
            end
                       
        end
        
        %% Viscous wrapper
        
        %         [vecCLv(1,ai), vecCD(1,ai), vecPREQ(1,ai), valVINF(1,ai), valLD(1,ai)] = fcnVISCOUS(vecCL(end,ai), vecCDI(end,ai), ...
        %             valWEIGHT, valAREA, valDENSITY, valKINV, vecDVENFREE, vecDVENIND, ...
        %             vecDVELFREE, vecDVELIND, vecDVESFREE, vecDVESIND, vecDVEPANEL, vecDVELE, vecDVEWING, vecN, vecM, vecDVEAREA, ...
        %             matCENTER, vecDVEHVCRD, vecAIRFOIL, flagVERBOSE, vecSYM, valVSPANELS, matVSGEOM, valFPANELS, matFGEOM, valFTURB, ...
        %             valFPWIDTH, valINTERF, vecDVEROLL);
        
    end
end
    
out = vecE(end,1,end)./vecE(end,2,end);

fp2 = fopen('optihistory3.txt','at');
fprintf(fp2,'%f %f ', out, orig2);
fprintf(fp2,'\r\n');
fclose(fp2);

end

