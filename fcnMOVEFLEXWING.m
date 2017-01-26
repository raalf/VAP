function [matNPVLST, matNPNEWWAKE] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, vecDVEHVSPN, vecDVELE, matNPVLST, matDEFGLOB,...
    matTWISTGLOB, valVINF, matSLOPE, valTIMESTEP, vecN, vecM, vecDVEWING, vecDVEPANEL, vecLSAC, matSCLST, vecDVEPITCH)
% matNEWWAKE, matNPNEWWAKE, 
% This function determines the velocities with which the DVEs are moved
% based on the deflection and twist of the wing. The corresponding
% translations are then computed of the DVE vertices and control points.

uinf = 1;

[ledves, ~, ~] = find(vecDVELE > 0);
[tedves, ~, ~] = find(vecDVETE > 0);

% Span of each spanwise set of DVEs
vecDVESPAN = 2*vecDVEHVSPN(ledves)';

% Deflection velocity after first timestep (referenced to zero initial
% deflection and twist)
if valTIMESTEP > 2

    % Calculate cartesian velocity of DVE edges
    del_twist = ((matTWISTGLOB(valTIMESTEP-1,3:sum(vecN,1)+3) - matTWISTGLOB(valTIMESTEP-2,3:sum(vecN,1)+3)));
    vecXVEL = repmat(uinf*cos(valALPHA)*cos(valBETA),1,sum(vecN,1)+1);
    vecYVEL = repmat(uinf*sin(valBETA),1,size(matSLOPE,2)) + [0,((matSLOPE(valTIMESTEP-1,2:end) - matSLOPE(valTIMESTEP-2,2:end))./valDELTIME).*vecDVESPAN.*cos(repmat(pi/2,1,size(matSLOPE,2)-1) - matSLOPE(valTIMESTEP-1,2:end))];
    vecZVEL = repmat(uinf*sin(valALPHA)*cos(valBETA),1,sum(vecN,1)+1) + ((matDEFGLOB(valTIMESTEP-1,3:sum(vecN,1)+3) - matDEFGLOB(valTIMESTEP-2,3:sum(vecN,1)+3))./valDELTIME);

    % Calculate cartesian velocity of DVE control points
%     for yy = 2:size(vecXVEL,2)
%         vecCENXVEL(yy-1) = (vecXVEL(yy)+vecXVEL(yy-1))/2;
%         vecCENZVEL(yy-1) = (vecZVEL(yy)+vecZVEL(yy-1))/2;
%     end
%     vecCENYVEL = repmat(uinf*sin(valBETA),1,size(matSLOPE,2)-1) + ((matSLOPE(valTIMESTEP-1,2:end) - matSLOPE(valTIMESTEP-2,2:end))./valDELTIME).*vecDVEHVSPN(ledves)'.*cos(repmat(pi/2,1,size(matSLOPE,2)-1) - matSLOPE(valTIMESTEP-1,2:end));
    
else
    
    % Calculate cartesian velocity of DVE edges
    del_twist = 5000*((matTWISTGLOB(valTIMESTEP-1,3:sum(vecN,1)+3)));
%     del_twist = (pi/2).*ones(1,sum(vecN,1)+1);
    vecXVEL = repmat(uinf*cos(valALPHA)*cos(valBETA),1,sum(vecN,1)+1);
    vecYVEL = repmat(uinf*sin(valBETA),1,size(matSLOPE,2)) + [0,((matSLOPE(valTIMESTEP-1,2:end))./valDELTIME).*vecDVESPAN.*cos(repmat(pi/2,1,size(matSLOPE,2)-1) - matSLOPE(valTIMESTEP-1,2:end))];
    vecZVEL = repmat(uinf*sin(valALPHA)*cos(valBETA),1,sum(vecN,1)+1) + ((matDEFGLOB(valTIMESTEP-1,3:sum(vecN,1)+3))./valDELTIME);

    % Calculate cartesian velocity of DVE control points
%     for yy = 2:size(vecXVEL,2)
%         vecCENXVEL(yy-1) = (vecXVEL(yy)+vecXVEL(yy-1))/2;
%         vecCENZVEL(yy-1) = (vecZVEL(yy)+vecZVEL(yy-1))/2;
%     end
%     vecCENYVEL = repmat(uinf*sin(valBETA),1,size(matSLOPE,2)-1) + ((matSLOPE(valTIMESTEP-1,2:end))./valDELTIME).*vecDVEHVSPN(ledves)'.*cos(repmat(pi/2,1,size(matSLOPE,2)-1) - matSLOPE(valTIMESTEP-1,2:end));

end
% uinf = repmat([uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)],sum(vecN,1)+1,1);
lepanels = vecDVEPANEL(ledves);

% Determine DVEs in each spanwise station
for i = 1:max(vecDVEWING)

	idxdve = ledves(vecDVEWING(ledves) == i);
	idxpanel = lepanels(vecDVEWING(ledves) == i);

    m = vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
    
    rows = repmat(idxdve,1,m) + tempm;

end

% Determine vertices that need to be moved at each spanwise station

% All left LE and TE points to move
temp_leftV = [matDVE(rows,1),matDVE(rows,4)];
temp_leftV = reshape(temp_leftV,sum(vecN,1),[]);

[move_row,~] = find(temp_leftV); % Vector correspond to which deflection velocity should be used for each element

translateNPVLST = zeros(size(matNPVLST,1),3);
temp_translate = zeros(size(matNPVLST,1),3);

temp_r = [sqrt(sum(matSCLST.^2,2)), zeros(length(matSCLST(:,1)),2)]; % Distance between vertex and shear center

xz_sign = sign(matSCLST(:,1)); % Determines whether positive or negative contribution to X and Z velocity for each vertex

test = temp_r(temp_leftV,1);
test2 = sin(del_twist(move_row))';

test_theta = vecDVEPITCH(ledves);

test_theta = [test_theta(1);test_theta]';

% Translate left edge vertices due to twist
twistXDIST = -temp_r(temp_leftV,1).*cos(del_twist(move_row)+test_theta(move_row))' +  temp_r(temp_leftV,1).*cos(test_theta(move_row))';%.*del_twist(move_row)'; % X component of twist (X = omega*r*sin(theta)*dt)
twistZDIST = temp_r(temp_leftV,1).*sin(del_twist(move_row)+test_theta(move_row))' -  temp_r(temp_leftV,1).*sin(test_theta(move_row))';%.*del_twist(move_row)'; % Z component of twist (Z = omega*r*cos(theta)*dt)

temp_translate(temp_leftV,1) = twistXDIST';
temp_translate(temp_leftV,3) = twistZDIST';

% Translate left edge vertices due to freestream and bending
translateNPVLST(temp_leftV,1) = zeros(length(twistXDIST),1);%valDELTIME.*vecXVEL(move_row);
translateNPVLST(temp_leftV,2) = zeros(length(twistXDIST),1);%100*valDELTIME.*vecYVEL(move_row);
translateNPVLST(temp_leftV,3) = zeros(length(twistXDIST),1);%-100*valDELTIME.*vecZVEL(move_row);

% All right LE and TE points to move
temp_rightV = [matDVE(rows,2), matDVE(rows,3)];
temp_rightV = reshape(temp_rightV,sum(vecN,1),[]);

[move_row,~] = find(temp_rightV); % Vector correspond to which deflection velocity should be used for each element

% Translate right edge vertices due to twist
twistXDIST = -temp_r(temp_rightV,1).*cos(del_twist(move_row+1)+test_theta(move_row+1))' + temp_r(temp_rightV,1).*cos(test_theta(move_row+1))';%.*del_twist(move_row+1)'; % X component of twist (X = omega*r*sin(theta)*dt)
twistZDIST = temp_r(temp_rightV,1).*sin(del_twist(move_row+1)+test_theta(move_row+1))' - temp_r(temp_rightV,1).*sin(test_theta(move_row+1))';%.*del_twist(move_row+1)'; % Z component of twist (Z = omega*r*cos(theta)*dt)

temp_translate(temp_rightV,1) = twistXDIST';
temp_translate(temp_rightV,3) = twistZDIST';

% Assign appropriate sign to twist movement
temp_translate(:,1) = xz_sign.*temp_translate(:,1);
temp_translate(:,3) = xz_sign.*temp_translate(:,3);

% Translate right edge vertices due to freestream and bending
translateNPVLST(temp_rightV,1) = zeros(length(twistXDIST),1);%valDELTIME.*vecXVEL(move_row+1);
translateNPVLST(temp_rightV,2) = zeros(length(twistXDIST),1);%100*valDELTIME.*vecYVEL(move_row+1);
translateNPVLST(temp_rightV,3) = zeros(length(twistXDIST),1);%-100*valDELTIME.*vecZVEL(move_row+1);

figure(3)
clf
plot(matNPVLST(temp_rightV(14,:),1), matNPVLST(temp_rightV(14,:),3),'-o')
hold on
% patch('Faces',matDVE,'Vertices',matNPVLST,'FaceColor','r')
% hold on

% Old trailing edge vertices
% matNEWWAKE(:,:,4) = matVLST(matDVE(vecDVETE>0,4),:);
% matNEWWAKE(:,:,3) = matVLST(matDVE(vecDVETE>0,3),:);

% Old non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,4) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,3) = matNPVLST(matDVE(vecDVETE>0,3),:);

% Update matVLST and matNPVLST
% matVLST = matVLST - translateNPVLST;
matNPVLST = matNPVLST - (translateNPVLST - temp_translate);

plot(matNPVLST(temp_rightV(14,:),1), matNPVLST(temp_rightV(14,:),3),'-o')

clf
patch('Faces',matDVE,'Vertices',matNPVLST,'FaceColor','b')
axis equal
axis tight
hold off
 
% % New trailing edge vertices
% matNEWWAKE(:,:,1) = matVLST(matDVE(vecDVETE>0,4),:);
% matNEWWAKE(:,:,2) = matVLST(matDVE(vecDVETE>0,3),:);
% 

% New non-planar trailing edge vertices (used to calculate matWADJE)
matNPNEWWAKE(:,:,1) = matNPVLST(matDVE(vecDVETE>0,4),:);
matNPNEWWAKE(:,:,2) = matNPVLST(matDVE(vecDVETE>0,3),:);



