function [matDEF, matTWIST] = fcnIMPLICIT(matEIx, matGJt, matDEF, matTWIST, vecJT, valDELTIME, vecDVEHVSPN, vecLSM,...
    vecLM, vecLIFTDIST, vecMOMDIST, vecSPANDIST, valTIMESTEP)

% vecLM(end) = 0;

nele = size(vecLM,1)-1;

matDAMP = [0*ones(size(vecLM,1),1), 0*ones(size(vecLM,1),1)];

vecJT = 0.00037078.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
    - 0.01102270.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
    + 0.12838255.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST - 0.73708913.*vecSPANDIST.*vecSPANDIST.*vecSPANDIST...
    + 2.15067037.*vecSPANDIST.*vecSPANDIST - 2.99312818.*vecSPANDIST + 1.84576176;

delx = sum(2*vecDVEHVSPN,1)/length(vecDVEHVSPN);
%% Bending matrix terms

a = 6.*matEIx(2:end,1)./(delx^4) + 6.*matEIx(2:end,2)./(delx^3) - 2.*matEIx(2:end,3)./(delx^2); % Main diagonal elements
b = -4.*matEIx(2:end-1,1)./(delx^4) - 6.*matEIx(2:end-1,2)./(delx^3) + matEIx(2:end-1,3)./(delx^2); % Elements 1 above main diagonal
c = matEIx(2:end-2,1)./(delx^4) + 2.*matEIx(2:end-2,2)./(delx^3); % Elements 2 above main diagonal
d = -4.*matEIx(3:end,1)./(delx^4) - 2.*matEIx(3:end,2)./(delx^3) + matEIx(3:end,3)./(delx^2); % Elements 1 below main diagonal
e = matEIx(4:end,1)./(delx^4); % Elements 2 below main diagonal

a(1) = a(1) + matEIx(2,1)./(delx^4); % Accounts for boundary condition at root of zero slope and zero displacement

% Accounting for boundary conditions at the tip of zero shear and zero
% moment
% a(end-1) = 6.*matEIx(end-1,1)./(delx^4) + 6.*matEIx(end-1,2)./(delx^3) - 2.*matEIx(end-1,3)./(delx^2) + vecLM(end-1)./(valDELTIME^2)...
%     - (matEIx(end-1,1)/(delx^4) + 2*matEIx(end-1,2)/(delx^3));
% a(end) = 6.*matEIx(end,1)./(delx^4) + 6.*matEIx(end,2)./(delx^3) - 2.*matEIx(end,3)./(delx^2) + vecLM(end)./(valDELTIME^2)...
%     + 3*(matEIx(end,1)/(delx^4) + 2*matEIx(end,2)/(delx^3)) + 2*(-4*matEIx(end,1)/(delx^4) - 6*matEIx(end,2)/(delx^3) + matEIx(end,3)/(delx^2));
% 
% b(end) = -4.*matEIx(end-1,1)./(delx^4) - 6.*matEIx(end-1,2)./(delx^3) + matEIx(end-1,3)./(delx^2) + 2*(matEIx(end-1,1)/(delx^4) + 2*matEIx(end-1,2)/(delx^3));
% 
% d(end) = -4.*matEIx(end,1)./(delx^4) - 2.*matEIx(end,2)./(delx^3) + matEIx(end,3)./(delx^2) - 2*(matEIx(end,1)/(delx^4) + 2*matEIx(end,2)/(delx^3))...
%     - (-4.*matEIx(end,1)./(delx^4) - 6.*matEIx(end,2)./(delx^3) + matEIx(end,3)./(delx^2));

A0 = diag(a,0);
A1 = diag(b,1);
A2 = diag(c,2);
A3 = diag(d,-1);
A4 = diag(e,-2);

matA = A0 + A1 + A2 + A3 + A4;

%% Torsion matrix terms

aprime = -2.*matGJt(2:end,1)./(delx^2) - matGJt(2:end,2)./delx; % Main diagonal elements
bprime = matGJt(2:end-1,1)./(delx^2); % Elements above main diagonal
cprime = matGJt(3:end,2)./delx + matGJt(3:end,1)./(delx^2); % Elements below main diagonal

cprime(end) = matGJt(end,2)./delx + 2*matGJt(end,1)/(delx^2); % Accounts for boundary of zero slope at wing tip

B0 = diag(aprime,0);
B1 = diag(bprime,1);
B2 = diag(cprime,-1);

matB = B0 + B1 + B2;

%% Set diagonal matrices for acceleration terms

S0 = [diag(-2*vecLM(2:end)./(valDELTIME^2) + 11*matDAMP(2:end,1)./(6*valDELTIME),0), diag(2*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0); diag(2*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0), diag(-2*vecJT(2:end)./(valDELTIME^2) + 11*matDAMP(2:end,2)./(6*valDELTIME),0)];
S1 = [diag(5*vecLM(2:end)./(valDELTIME^2),0) - 18*matDAMP(2:end,1)./(6*valDELTIME), diag(-5*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0); diag(-5*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0), diag(5*vecJT(2:end)./(valDELTIME^2) - 18*matDAMP(2:end,2)./(6*valDELTIME),0)];
S2 = [diag(-4*vecLM(2:end)./(valDELTIME^2),0) + 9*matDAMP(2:end,1)./(6*valDELTIME), diag(4*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0); diag(4*vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0), diag(4*vecJT(2:end)./(valDELTIME^2) + 9*matDAMP(2:end,2)./(6*valDELTIME),0)];
S3 = [diag(vecLM(2:end)./(valDELTIME^2),0) - 2*matDAMP(2:end,1)./(6*valDELTIME), diag(-vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0); diag(-vecLM(2:end).*vecLSM(2:end)./(valDELTIME^2),0), diag(-vecJT(2:end)./(valDELTIME^2) - 2*matDAMP(2:end,2)./(6*valDELTIME),0)];


% Assemble global matrix to solve for displacements
matC = [matA, zeros(size(vecLM,1)-1); zeros(size(vecLM,1)-1), matB];

matD = matC - S0;

matD(nele-1,nele-1) = 6.*matEIx(end-1,1)./(delx^4) + 6.*matEIx(end-1,2)./(delx^3) - 2.*matEIx(end-1,3)./(delx^2) + vecLM(end-1)./(valDELTIME^2)...
    - (matEIx(end-1,1)/(delx^4) + 2*matEIx(end-1,2)/(delx^3));

matD(nele,nele) = 6.*matEIx(end,1)./(delx^4) + 6.*matEIx(end,2)./(delx^3) - 2.*matEIx(end,3)./(delx^2) + vecLM(end)./(valDELTIME^2)...
    + 3*(matEIx(end,1)/(delx^4) + 2*matEIx(end,2)/(delx^3)) + 2*(-4*matEIx(end,1)/(delx^4) - 6*matEIx(end,2)/(delx^3) + matEIx(end,3)/(delx^2));

matD(nele-1,nele) = -4.*matEIx(end-1,1)./(delx^4) - 6.*matEIx(end-1,2)./(delx^3) + matEIx(end-1,3)./(delx^2) + 2*(matEIx(end-1,1)/(delx^4) + 2*matEIx(end-1,2)/(delx^3));

matD(nele,nele-1) = -4.*matEIx(end,1)./(delx^4) - 2.*matEIx(end,2)./(delx^3) + matEIx(end,3)./(delx^2) - 2*(matEIx(end,1)/(delx^4) + 2*matEIx(end,2)/(delx^3))...
    - (-4.*matEIx(end,1)./(delx^4) - 6.*matEIx(end,2)./(delx^3) + matEIx(end,3)./(delx^2));

%% Solve for RHS of equation that includes external forces

if valTIMESTEP == 4

    vecRES = (matD + S2 + 8*S3)\([vecLIFTDIST(2:end)' - vecLM(2:end).*9.81; vecMOMDIST(2:end) - vecLM(2:end).*vecLSM(2:end).*9.81]);

    matDEF(valTIMESTEP,:) = vecRES(1:nele);
    matTWIST(valTIMESTEP,:) = vecRES(nele+1:end);

    matDEF(1,:) = -8*vecRES(1:nele);
    matTWIST(1,:) = -8*vecRES(nele+1:end);

    matDEF(2,:) = -vecRES(1:nele);
    matTWIST(2,:) = -vecRES(nele+1:end);

    matDEF(3,:) = zeros(1,nele);
    matTWIST(3,:) = zeros(1,nele);

else

    vecR = S1*[matDEF(end,:)'; matTWIST(end,:)'] + S2*[matDEF(end-1,:)'; matTWIST(end-1,:)'] + S3*[matDEF(end-2,:)'; matTWIST(end-2,:)'] +...
        [vecLIFTDIST(2:end)' - vecLM(2:end).*9.81; vecMOMDIST(2:end) - vecLM(2:end).*vecLSM(2:end).*9.81];

    % Solve structural dynamics matrix. Resulting vector will have deflection
    % results stack on top of the twist results
    vecRES = matD\vecR;

    % Extract individual results from vecRES
    matDEF(valTIMESTEP,:) = vecRES(1:nele);
    matTWIST(valTIMESTEP,:) = vecRES(nele+1:end);

end

end
