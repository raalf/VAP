function [nfreecs,thrustCFfree, axialCFfree, sideCFfree] = fcnRCROSSFLOWFORCE(valNELE, vecTHETA, vecDVETESWP,vecDVELESWP,vecDVEHVSPN,vecDVEHVCRD,matDVE,matUINF,matVLST,B,C)
% This function calculated the force due to cross flow. This occurs in
% forward flight at angles other than 90 deg. This is done but integrating
% vorticity across the a DVE, calculating the lift across a vortex sheet.

% OUTPUTS
% nfreecs - normal force due the freestream in the crossflow direction
% thrustCFfree - thrust due to the crossflow at each control point


% Calculate UxXi
% Calculate leading edge and trailing control points 
tempCPLE = matVLST(matDVE(:,1),:)+(matVLST(matDVE(:,2),:)-matVLST(matDVE(:,1),:))/2;
tempCPTE = matVLST(matDVE(:,3),:)+(matVLST(matDVE(:,4),:)-matVLST(matDVE(:,3),:))/2;
% Create Xi unit vectors at control points
matXi = (tempCPLE-tempCPTE)./(2*vecDVEHVCRD);

% UxXi
tempb = cross(matUINF,matXi);
UxXi = sqrt(sum(abs(tempb).^2,2));

% Normal Directions
en = tempb.*repmat((1./UxXi),1,3);
% Side direction
es = [-sin(vecTHETA) cos(vecTHETA) zeros(valNELE,1)];
% Axial direction
ea = [cos(vecTHETA) sin(vecTHETA) zeros(valNELE,1)];

% Cross-flow force:
% nfreecs = 2*(U x Xi)*(2*Xi*B*eta-0.5*(tan(Zeta_LE)+tan(Zeta_TE))*(4/3)*C*eta*eta*eta)
nfreecs = (2*UxXi'.*(2*vecDVEHVCRD'.*B.*vecDVEHVSPN'-(2/3)*(tan(vecDVELESWP')+tan(vecDVETESWP')).*C.*vecDVEHVSPN'.*vecDVEHVSPN'.*vecDVEHVSPN'))';
thrustCFfree = nfreecs.*(en(:,3));
sideCFfree = nfreecs.*(dot(es,en,2));
axialCFfree = -1*nfreecs.*(dot(ea,en,2));

% Test Plotting
%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),nfreecs.*en(:,1),nfreecs.*en(:,2),nfreecs.*en(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),thrustCWfree(:,1),thrustCWfree(:,2),thrustCWfree(:,3))

% THESE COMMENTS ARE FOR LATER REFERENCE
% THIS IS NO LONGER NEEDED USING A,B & C from fcnRDVENFORCE
%
% idx1 = vecDVETIP == 2;
% idx1 = vecDVELE == 1;
%
% % Find B, C values required
% A = zeros(1,valNELE);
% B = zeros(1,valNELE);
% C = zeros(1,valNELE);
% 
% % Calculate the A,B,C values
% %A(idx1) = matCOEFF(idx1,1); % A is not used
% B(idx1) = matCOEFF(idx1,2);
% C(idx1) = matCOEFF(idx1,3);
% 
% % % If any other row, A= A-Aupstream, B= B-Bupstream, C= C-Cupstream (cross
% % % flow dir.) This subtraction is done
% % dvenum = find(idx1==0); %dvenum in question
% % idxf = matADJE((ismember(matADJE(:,1), dvenum) & matADJE(:,2) == 2),3); 
% % A(idx1 ==0) = (matCOEFF(idx1==0,1)-matCOEFF(idxf,1));
% % B(idx1 ==0) = (matCOEFF(idx1==0,2)-matCOEFF(idxf,2));
% % C(idx1 ==0) = (matCOEFF(idx1==0,3)-matCOEFF(idxf,3));
% 
% % if any other row, A= A-Aupstream, B= B-Bupstream, C= C-Cupstream
% dvenum = find(idx1==0); %dvenum in question
% idxf = matADJE((ismember(matADJE(:,1), dvenum) & matADJE(:,2) == 1),3); %upstream dve num
% %A(idx1 ==0) = (matCOEFF(idx1==0,1)-matCOEFF(idxf,1)); A is not used
% B(idx1 ==0) = (matCOEFF(idx1==0,2)-matCOEFF(idxf,2));
% C(idx1 ==0) = (matCOEFF(idx1==0,3)-matCOEFF(idxf,3));
