function [nind, nfree, thrustind, thrustfree, sideind, sidefree, axialind, axialfree, A, B, C] = fcnRDVENFORCE(valWSIZE, valTIMESTEP, valNELE, valWNELE, seqALPHAR, vecDVEPITCH, vecK, vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecDVEYAW, vecDVEMCSWP, vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecDVEROLL,  vecDVEHVCRD, vecDVELE, vecDVEHVSPN, vecWDVEPITCH, vecDVELESWP, vecDVETESWP, vecSYM, vecTHETA, matVLST, matDVE, matUINF, matCOEFF, matADJE, matWDVE, matWVLST, matCENTER, matWCOEFF)

% A modified DVENFORCE function that has been tailored to calculate thrust
% and axial force. Mostly transfered directly from fcnDVENFORCE.
%
% Modifications include:
%  - Using the unqiue velocities at each DVE
%  - Calculating thrust, side and axial forces from the original nromal
%  force calculation.
%  - Thrust is define in the global z-direction
%  - Side is defined in the global y-direction
%  - Axial is defined in the global x-direction


idx1 = vecDVELE == 1; %index of LE vectors

% Calculate s matrix
s = zeros(valNELE,3);
s(idx1,:) =( matVLST(matDVE(idx1,2),:) -matVLST(matDVE(idx1,1),:) )  ./ repmat((vecDVEHVSPN(idx1).*2),1,3);

% vector along midchord for remaining elements
if any(idx1 == 0)
    tempa = zeros(valNELE,1);
    tempa(idx1==0,1) = tan(vecDVEMCSWP(idx1==0));
    tempa= [tempa ones(valNELE,1)  zeros(valNELE,1)];
    
    % global vectors to find ind. velocities
    s(idx1 ==0,:)= fcnSTARGLOB(tempa(idx1==0,:), vecDVEROLL(idx1==0), vecDVEPITCH(idx1==0), vecDVEYAW(idx1==0));
end

% 80% of halfspan
eta8 = vecDVEHVSPN*0.8;
% matUINF(:,1) = abs(matUINF(:,1));
% matUINF(:,2) = abs(matUINF(:,2));
% UxS
tempb = cross(matUINF,s,2);

% norm(UxS)
uxs = sqrt(sum(abs(tempb).^2,2));

% eN = tempa.*(1/UxS);
en = tempb.*repmat((1./uxs),1,3);

% Thrust direction
et = repmat([0 0 1],[valNELE,1]);

% Side direction
es = [-sin(vecTHETA) cos(vecTHETA) zeros(valNELE,1)];

% Axial direction
ea = [cos(vecTHETA) sin(vecTHETA) zeros(valNELE,1)];

%% Thrust due to freestream

% N_free = (A*2*eta + C/3*2*eta*eta*eta)*UxS;
% if first row, A=A, B=B, C=C
A = zeros(1,valNELE);
B = zeros(1,valNELE);
C = zeros(1,valNELE);

A(idx1) = matCOEFF(idx1,1);
B(idx1) = matCOEFF(idx1,2);
C(idx1) = matCOEFF(idx1,3);
% if any other row, A= A-Aupstream, B= B-Bupstream, C= C-Cupstream
dvenum = find(idx1==0); %dvenum in question
idxf = matADJE((ismember(matADJE(:,1), dvenum) & matADJE(:,2) == 1),3); %upstream dve num
A(idx1 ==0) = (matCOEFF(idx1==0,1)-matCOEFF(idxf,1));
B(idx1 ==0) = (matCOEFF(idx1==0,2)-matCOEFF(idxf,2));
C(idx1 ==0) = (matCOEFF(idx1==0,3)-matCOEFF(idxf,3));


nfree = ((A .*2 .* vecDVEHVSPN'+  C./3.*2.*vecDVEHVSPN'.*vecDVEHVSPN'.*vecDVEHVSPN') .*uxs')';

%% Induced force

% for first row (m=1):
%	compute 3 velocities along LE of DVE

% for remaining rows (m>1):
%	3 velocities are average of element center and upstream DVE center

% element leading edge midpoint of LE elements only:
xle = (matVLST(matDVE(idx1,1),:) + matVLST(matDVE(idx1,2),:))/2;

% fpg will be all points that we grab velocity at. It will be
% valNELE x XYZ x 3 for now, then we will reshape after
fpg = zeros(valNELE,3,3);

% leading edge row:
fpg(idx1,:,1) = (xle + s(idx1==1,:).*repmat(-eta8(idx1==1),1,3)); %left side
fpg(idx1,:,2) = xle ; %middle
fpg(idx1,:,3) = (xle + s(idx1==1,:).*repmat(eta8(idx1==1),1,3)); %right ride

% Remaining rows:
if any(idx1 == 0)
    fpg(idx1==0,:,1) = (matCENTER(idx1==0,:) + s(idx1==0,:).*repmat(-eta8(idx1==0),1,3)); %left side
    fpg(idx1==0,:,2) = matCENTER(idx1==0,:) ; %middle
    fpg(idx1==0,:,3) = (matCENTER(idx1==0,:) + s(idx1==0,:).*repmat(eta8(idx1==0),1,3)); %right ride
end

% if there are any elements in not the first row, we will append all the LE center
% points to the bottom, necessary for averaging. the if statement will be ignored if all m=1.
% need to remove tan, we should already have the vector
if any(idx1 == 0)
    fpg(1:(valNELE+sum(idx1)),:,1) = [fpg(1:valNELE,:,1) ; (matCENTER(idx1==1,:) + repmat(tan(vecDVEMCSWP(idx1==1)).*-eta8(idx1==1),1,3))]; %left side
    fpg(1:(valNELE+sum(idx1)),:,2) = [fpg(1:valNELE,:,2) ; matCENTER(idx1==1,:)] ; %middle
    fpg(1:(valNELE+sum(idx1)),:,3) = [fpg(1:valNELE,:,3) ; (matCENTER(idx1==1,:) + repmat(tan(vecDVEMCSWP(idx1==1)).*eta8(idx1==1),1,3))]; %right ride
end

% take second dimension, move to bottom. then take third dimension and move
% to bottom
fpg = reshape(permute(fpg,[1 3 2]),[],3);

% get velocities
[w_total] = fcnINDVEL(fpg,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
    vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
    matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
    vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);

% undo reshape and permute
% matrix is now LE vels for all LE elements, center vels for remaining DVES,
% and the center vels for the LE elements are appended to the bottom
w_total = permute(reshape(w_total,[],3,3),[1 3 2]);

%grab LE DVE values into final w matrix
w = zeros(valNELE,3,3);
w(idx1,:,:) = w_total(idx1,:,:);

%make a matrix with all DVE center velocities
if any(idx1 ==0)
    w_center = zeros(valNELE,3,3);
    w_center(idx1 ==0,:,:) = w_total(idx1 ==0,:,:);
    w_center(idx1,:,:) = w_total(valNELE+1:end,:,:); %add center vels from LE DVES
    
    %//case of multiple lifting lines along the span
    %//the induced velocity at the lifting line is averaged with the
    %//velocities at mid chord locations of the DVES upstream and
    %//downstream of the bound vortex. Otherwise, the singularity of
    %//the bound vortex and the discontinuity of the bound vortex sheet
    %//of a blade with twist causes trouble.
    w(idx1 ==0,:,:) = (w_center(idx1 ==0,:,:)+w_center(idxf,:,:))./2;
end

% perform integration
tempd = cross(w,repmat(s,[1,1,3]),2);
gamma(:,1) = A - B.*eta8' + C.*eta8'.*eta8';
gamma(:,2) = A;
gamma(:,3) = A + B.*eta8' + C.*eta8'.*eta8';
tempr = tempd .* repmat(permute(gamma,[1 3 2]),1,3,1);

%//The resulting induced force is
%//determined by numerically integrating forces across element
%//using Simpson's Rule with overhaning parts
%  R = (R1 + 4*Ro+R2).*eta8/3;
r = (tempr(:,:,1) + 4.*tempr(:,:,2) + tempr(:,:,3)) .* repmat(eta8,1,3) ./3;
%  R = R + ((7.*R1 - 8.*Ro + 7.*R2).*(eta-eta8)./3);
r = r + ((7.*tempr(:,:,1) - 8.*tempr(:,:,2) + 7.*tempr(:,:,3)).*repmat((vecDVEHVSPN-eta8),1,3)./3);

% Induced normal force
nind = dot(r,en,2);

% Induced forces
thrustind = dot(r,et,2);
sideind = dot(r,es,2);
axialind = dot(r,ea,2);

% Freestream forces
thrustfree = nfree.*(en(:,3));
sidefree = nfree.*(dot(es,en,2));
axialfree = nfree.*(dot(ea,en,2));

% Test plotting
%hold on
%quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),es(:,1),es(:,2),es(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),matUINF(:,1),matUINF(:,2),matUINF(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),s(:,1),s(:,2),s(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),A.*en(:,1),A.*en(:,2),A.*en(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),(thrustfree+thrustind).*et(:,1),(thrustfree + thrustind).*et(:,2),(thrustind+thrustfree).*et(:,3),'k')
%quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3),A(:,1),A(:,2),A(:,3))
%quiver3(matCENTER(:,1),matCENTER(:,2), matCENTER(:,3),ea(:,1),ea(:,2),ea(:,3))
%clf
%plot((matCENTER(:,1)/abs(matCENTER(:,1)))*sqrt(matCENTER(:,1).^2+matCENTER(:,1).^2),thrustfree+thrustind,'*')
end