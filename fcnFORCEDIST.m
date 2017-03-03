function [vecLIFTDIST, vecMOMDIST] = fcnFORCEDIST(liftfree, liftind, matSCLST, valDENSITY, valWEIGHT, valCL, vecDVEHVSPN, vecLEDVES, vecN, vecM,...
    vecDVEWING, vecDVEPANEL, matCENTER, vecSPANDIST, matNPVLST, matNPDVE, matSC, matLIFTDIR, vecMAC, valCM, valAREA, vecSPNWSEAREA,vecLSAC)

% This function computes the dimensional force and moment distribution
% across the wing, resolved to the shear center line. Moment is taken
% about the elastic axis.
%
% INPUT:
%
% OUTPUT:
% vecLIFTDIST - 1 x sum(vecN) matrix of the total lift at each spanwise
% station
% vecMOMDIST - 1 x sum(vecN) matrix of the total pitching moment at each
% spanwise station (+ve nose up)

% Calculate qinf required for steady level flight. This will be used for
% load calculations
q_inf = valWEIGHT/(valCL*valAREA);
% valVINF = sqrt(2*q_inf/valDENSITY);

[matROWS] = fcnDVEROW(vecLEDVES, vecDVEPANEL, vecDVEWING, vecM, vecN);

% ======================================================================= % 
%                        Start Lift Calculation                           %
% ======================================================================= %

% Convert DVE force from force/density to force/unit length
vecLIFTDIST = ((sum(liftfree(matROWS),2) + sum(liftind(matROWS),2))*valDENSITY)./(2*vecDVEHVSPN(vecLEDVES));

% Interpolate lift distribution to find lift at DVE edges (i.e. Structural
% grid points)
tempSPANDIST = matCENTER(vecLEDVES,2); % Y coordinate of DVE midpoints

% Some setup work to be able to perform linear interpolation without a
% for loop
tempSPANDIST = repmat(tempSPANDIST', size(tempSPANDIST,1),1);

tempSPANDIST = triu(tempSPANDIST);

vecLIFTDIST = repmat(vecLIFTDIST', size(vecLIFTDIST,1),1);

vecLIFTDIST = triu(vecLIFTDIST);

vecLIFTDIST = ((vecSPANDIST(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
    tempSPANDIST(1,1:(end-1)))).*(vecLIFTDIST(2,2:end)-vecLIFTDIST(1,1:(end-1))) + vecLIFTDIST(1,1:(end-1)); % Linear interpolation of lift

vecLIFTDIST = [vecLIFTDIST(1), vecLIFTDIST, 0]; % Add zero lift to tip. Use the same lift at root as at the next neighbouring node. This doesn't really have any impact since the structure solver doesn't use this value

% ======================================================================= % 
%                        Start Moment Calculation                         %
% ======================================================================= %
lift_chord = liftfree(matROWS) + liftind(matROWS); % Lift distribution along each chordwise location

% Interpolate lift force at each chordwise station to use for moment
% calculations. Also determine the LE coordinates of each chordwise DVE to
% be used to calculate a moment arm between DVE and elastic axis
for i = 1:vecM(1)
    
    row_ledves(:,:,i) = [matNPVLST(matNPDVE(matROWS(:,i),1),:); matNPVLST(matNPDVE(matROWS(end,i),2),:)]; % LE coordinates of each set of chordwise DVEs
    
    temp_liftdist = repmat(lift_chord(:,i)', size(lift_chord,1),1);

    temp_liftdist = triu(temp_liftdist);

    temp_liftdist = ((vecSPANDIST(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
        tempSPANDIST(1,1:(end-1)))).*(temp_liftdist(2,2:end)-temp_liftdist(1,1:(end-1))) + temp_liftdist(1,1:(end-1)); % Linear interpolation of lift
    
    temp_lift = temp_liftdist;
   
    temp_lift = [temp_lift(1); temp_liftdist']; % Use the same lift at root as at the next neighbouring node
    
    % Lift force in X, Y, Z components to be used in cross product to
    % calculate moment about elastic axis. Each step into the 3rd dimension
    % is the lift distribution at each chordwise station
    lift_moment(:,:,i) = [temp_lift.*matLIFTDIR(matROWS(:,i),:); zeros(1,3)]; 
    
end

matSC = repmat(matSC,1,1,vecM(1));

% Compute moment arm for cross product
% delX = (row_ledves - matSC);
delX = [vecLSAC, zeros(size(vecLSAC,1),2)];
force = [zeros(size(vecLIFTDIST,2),2),vecLIFTDIST'];
% tempMOMDIST = cross(-delX,force); % M' = delx X lift
tempMOMDIST = vecLIFTDIST'.*vecLSAC;

% Compute magnitude of moment at each chordwise location
% for i = 1:vecM(1)
%    
%     matMOMDIST(:,i) = sign(tempMOMDIST(:,2,i)).*sqrt(sum(abs(tempMOMDIST(:,:,i)).^2,2));
%     
% end
% 
% vecMOMDIST = sum(matMOMDIST,2);
% vecMOMDIST = [vecMOMDIST(1:(end-1)).*valDENSITY./(2*vecDVEHVSPN(vecLEDVES)); 0] + [q_inf*vecMAC.*vecMAC.*valCM;0];

vecMOMDIST = tempMOMDIST + [q_inf*vecMAC.*vecMAC.*valCM;0];

end