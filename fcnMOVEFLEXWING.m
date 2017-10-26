function [matNPVLST, matNPNEWWAKE, matNEWWAKE, valUINF] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, vecDVEHVSPN, vecDVELE, matNPVLST, matDEFGLOB,...
    matTWISTGLOB, matSLOPE, valTIMESTEP, vecN, vecM, vecDVEWING, vecDVEPANEL, matSCLST, vecDVEPITCH, matNPDVE, vecSPANDIST, vecCL, valWEIGHT, valAREA, valDENSITY,valUINF)
% matNEWWAKE, matNPNEWWAKE,
% This function determines the velocities with which the DVEs are moved
% based on the deflection and twist of the wing. The corresponding
% translations are then computed of the DVE vertices and control points.

% q_inf = valWEIGHT/(vecCL(valTIMESTEP-1)*valAREA);
% valUINF = sqrt(2*q_inf/valDENSITY);

[ledves, ~, ~] = find(vecDVELE > 0);

% Span of each spanwise set of DVEs
vecDVESPAN = 2*vecDVEHVSPN(ledves)';

% Deflection velocity after first timestep (referenced to zero initial
% deflection and twist)

% Calculate cartesian velocity of DVE edges
del_twist = ((matTWISTGLOB(valTIMESTEP,:) - matTWISTGLOB(valTIMESTEP-1,:)));
omega = ((matTWISTGLOB(valTIMESTEP,:) - matTWISTGLOB(valTIMESTEP-1,:)))./valDELTIME;
vecXVEL = repmat(valUINF*cos(valALPHA)*cos(valBETA),1,sum(vecN,1)+1);
% vecYVEL = repmat(valUINF*sin(valBETA),1,size(matSLOPE,2)) + [0,((matSLOPE(valTIMESTEP,2:end) - matSLOPE(valTIMESTEP-1,2:end))./valDELTIME).*vecDVESPAN.*cos(repmat(pi/2,1,size(matSLOPE,2)-1) - matSLOPE(valTIMESTEP,2:end))];
% vecYVEL = repmat(valUINF*sin(valBETA),1,size(matSLOPE,2)) + [0, (vecDVESPAN - vecDVESPAN.*cos(matSLOPE(valTIMESTEP,2:end)-matSLOPE(valTIMESTEP-1,2:end)))./valDELTIME];
vecYVEL = repmat(valUINF*sin(valBETA),1,size(matSLOPE,2)+1) + [0, (matDEFGLOB(valTIMESTEP,2:end)-matDEFGLOB(valTIMESTEP-1,2:end))./...
    (valDELTIME.*tan(repmat(pi/2,1,size(matSLOPE,2))-(matSLOPE(valTIMESTEP,:)-matSLOPE(valTIMESTEP-1,:))./2))];
vecZVEL = repmat(valUINF*sin(valALPHA)*cos(valBETA),1,sum(vecN,1)+1) + ((matDEFGLOB(valTIMESTEP,:) - matDEFGLOB(valTIMESTEP-1,:))./valDELTIME);

% Determine DVEs in each spanwise station
[matROWS] = fcnDVEROW(ledves, vecDVEPANEL, vecDVEWING, vecM, vecN);

%% Determining displacement distances

% Determine vertices that need to be moved at each spanwise station

% All left LE and TE points to move
temp_leftV = [matNPDVE(matROWS,1),matNPDVE(matROWS,4)];
temp_leftV = reshape(temp_leftV,sum(vecN,1),[]);

[move_row,~] = find(temp_leftV); % Vector correspond to which index of deflection velocity matrix should be used for each element

% Allocate space for translation matrices
translateNPVLST = zeros(size(matNPVLST,1),3);
temp_translate = zeros(size(matNPVLST,1),3);

temp_r = [sqrt(sum(matSCLST.^2,2)), zeros(length(matSCLST(:,1)),2)]; % Distance between vertex and shear center

xz_sign = sign(matSCLST(:,1)); % Determines whether positive or negative contribution to X and Z velocity for each vertex

% Perform a linear interpolation/extrapolation to determine the pitch of
% the DVE's at their respective left and right edges. This is used to
% determine the initial orientation of the DVE before applying the twist
% caused by elastic deformation
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

% ======================== Left Edge Displacements ========================
% Translate left edge vertices due to twist
twistXDIST = -temp_r(temp_leftV,1).*cos(del_twist(move_row)+vecEDGEPITCH(move_row))' + ...
    temp_r(temp_leftV,1).*cos(vecEDGEPITCH(move_row))'; % X component of twist 
twistZDIST = temp_r(temp_leftV,1).*sin(del_twist(move_row)+vecEDGEPITCH(move_row))' -  ...
    temp_r(temp_leftV,1).*sin(vecEDGEPITCH(move_row))'; % Z component of twist

v_rot = temp_r(temp_leftV,1).*omega(move_row)';

% Assign twist displacement to translation matrix
temp_translate(temp_leftV,1) = twistXDIST;
temp_translate(temp_leftV,3) = twistZDIST;

test = zeros(size(matNPVLST,1),3);
test(temp_leftV,1) = v_rot.*sin(matTWISTGLOB(valTIMESTEP,move_row)+vecEDGEPITCH(move_row))'.*valDELTIME;
test(temp_leftV,3) = v_rot.*cos(matTWISTGLOB(valTIMESTEP,move_row)+vecEDGEPITCH(move_row))'.*valDELTIME;

% Translate left edge vertices due to freestream and bending
translateNPVLST(temp_leftV,1) = valDELTIME.*vecXVEL(move_row);
translateNPVLST(temp_leftV,2) = valDELTIME.*vecYVEL(move_row);
translateNPVLST(temp_leftV,3) = -1*valDELTIME.*vecZVEL(move_row);

% ======================== Right Edge Displacements =======================
% All right LE and TE points to move
temp_rightV = [matNPDVE(matROWS,2), matNPDVE(matROWS,3)];
temp_rightV = reshape(temp_rightV,sum(vecN,1),[]);

[move_row,~] = find(temp_rightV); % Vector correspond to which index of deflection velocity matrix should be used for each element

v_rot = temp_r(temp_leftV,1).*omega(move_row+1)';
% Translate right edge vertices due to twist
twistXDIST = -temp_r(temp_rightV,1).*cos(del_twist(move_row+1)+vecEDGEPITCH(move_row+1))' + ...
    temp_r(temp_rightV,1).*cos(vecEDGEPITCH(move_row+1))'; % X component of twist
twistZDIST = temp_r(temp_rightV,1).*sin(del_twist(move_row+1)+vecEDGEPITCH(move_row+1))' - ...
    temp_r(temp_rightV,1).*sin(vecEDGEPITCH(move_row+1))'; % Z component of twist 

temp_translate(temp_rightV,1) = twistXDIST';
temp_translate(temp_rightV,3) = twistZDIST';

test(temp_rightV,1) = v_rot.*sin(matTWISTGLOB(valTIMESTEP,move_row+1)+vecEDGEPITCH(move_row+1))'.*valDELTIME;
test(temp_rightV,3) = v_rot.*cos(matTWISTGLOB(valTIMESTEP,move_row+1)+vecEDGEPITCH(move_row+1))'.*valDELTIME;

% Assign appropriate sign to twist movement
temp_translate(:,1) = xz_sign.*temp_translate(:,1);
temp_translate(:,3) = xz_sign.*temp_translate(:,3);

% ======================================================================= %
% =============== TRYING NEW TRANSLATION FROM TWIST ===================== %
% ======================================================================= %
% temp_translate(:,1) = xz_sign.*test(:,1);
% temp_translate(:,3) = xz_sign.*test(:,3);
% ======================================================================= %
% ======================================================================= %

% Translate right edge vertices due to freestream and bending
translateNPVLST(temp_rightV,1) = valDELTIME.*vecXVEL(move_row+1);
translateNPVLST(temp_rightV,2) = valDELTIME.*vecYVEL(move_row+1);
translateNPVLST(temp_rightV,3) = -1*valDELTIME.*vecZVEL(move_row+1);

%% Move wing and generate new wake elements

% Old trailing edge vertices
matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,4) = matNPVLST(matNPDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matNPDVE(vecDVETE>0,3),:);

% Update matVLST and matNPVLST
matNPVLST = matNPVLST - (translateNPVLST - temp_translate);

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNPVLST(matNPDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNPVLST(matNPDVE(vecDVETE>0,3),:);

