function [matDEF, matTWIST] = fcnWINGTWISTBEND(valDENSITY,valDELTIME,valSPAN,valAREA,valTIMESTEP,vecDVEHVSPN,vecDVEHVCRD,...
    vecLEDVES,vecLSAC,vecJT,vecLSM,vecLAMBDA,vecLIFTDIST,vecMOMDIST,valUINF,matEIx,matGJt,matDEF,matTWIST,vecLM)

%--- IMPLICT RECURRENCE MATRIX SOLUTION FOR STRUCTURAL DYNAMIC RESPONSE ---%
% vecLAMBDA = [0.09; 0.18;  0.17; 0.16; 0.16; 0.16];
% matEIx(:,1) = [28976640003; 20069006224; 11805297779; 5580686223; 2414720000; 724416000.1];
% vecDVEHVSPN = [101; 101; 90; 90; 90; 90]./2;
% vecCRD = [154; 136; 118; 102; 85; 68];
% valSPAN = 560*2;
% vecMASS = [10769; 6050; 1433; 382; 201; 118];
% valDENSITY = 0.0765/12^3;
% valUINF = 3700;
% valAR = 10;
% num_el = 6;

%% Determining coefficients using eqn C8
num_el = size(vecLEDVES,1);
vecCRD = 2*vecDVEHVCRD;
valAR = valSPAN*valSPAN/valAREA;
vecMASS = (2*vecDVEHVSPN).*vecLM;

valBETA = (valAR/(2+valAR))*pi*valDENSITY*valUINF; % Forward speed and aspect ratio factor for wing

i = 2:num_el; % Spatial index

% Pre-allocate coefficients used in H1 matrix eqns
A = zeros(length(i)+1,1);
B = zeros(length(i)+1,1);
C = zeros(length(i)+1,1);
D = ones(length(i)+1,1);

% Calculate coefficients using moments of inertia at spanwise locations
A(i,1) = (1/4)*matEIx(1,1)./matEIx(i-1,1) + (1/12)*matEIx(1,1)./matEIx(i,1);
B(i,1) = (1/3)*matEIx(1,1)./matEIx(i-1,1) + (1/6)*matEIx(1,1)./matEIx(i,1);
C(i,1) = (1/12)*matEIx(1,1)./matEIx(i-1,1) + (1/12)*matEIx(1,1)./matEIx(i,1);
D(i,1) = (1/6)*matEIx(1,1)./matEIx(i-1,1) + (1/3)*matEIx(1,1)./matEIx(i,1);

%% Determine H2 matrix elements from eqn C10

% Bunch of crap to make an upper triangular matrix
temp1 = sum(triu(repmat(vecLAMBDA(2:length(vecLAMBDA)),1,length(vecLAMBDA)-1)),1);
temp1 = repmat(temp1,length(vecLAMBDA)-1,1);
temp2 = zeros(length(vecLAMBDA));
temp2(2:end,2:end) = triu(repmat(vecLAMBDA(2:length(vecLAMBDA)),1,length(vecLAMBDA)-1));
temp2(1,:) = [];
temp2(:,end) = [];
temp3 = sum(temp2,1);
temp4 = temp1 - repmat(temp3',1,length(vecLAMBDA)-1);

% Final H2 matrix from eqn C10
H2 = triu(temp4);

% Mirror the H2 matrix across the diagonal for use later on
[n,~]=size(H2);
temp_sum=H2'+H2;
temp_sum(1:n+1:end)=diag(H2);
temp_sum = tril(temp_sum);

temp_sum_2 = zeros(num_el-1);
temp_sum_2(2:end,1:end-1) = temp_sum(2:end,2:end);

%% Determine H1 matrix elements from eqn C7

m = 2:num_el;

% General eqn to compute elements of H1 matrix
a = repmat((vecLAMBDA(m-1,1)'.^2.*C(m-1,1)'),num_el-1,1) + repmat(vecLAMBDA(m-1,1)'.*D(m-1,1)',num_el-1,1).*temp_sum +...
    tril(repmat(vecLAMBDA(m,1)'.^2.*A(m,1)',num_el-1,1)) + repmat(vecLAMBDA(m,1)'.*B(m,1)',num_el-1,1).*temp_sum_2;

% Create lower triangular matrix with one upper diagonal row of ones. This
% is used to filter out the unnecessary values from the "a" matrix. The
% unnecessary elements are multiplied by the zeros while the elements to
% keep are multiplied by 1.
xx = tril(ones(num_el-2));
yy = ones(num_el-1);
yy(1:end-1,2:end) = xx;

% Final H1 matrix from eqn C9
H1 = yy.*a;

%% Calculate little b matrix from C16

b_inner = (matEIx(1,1)/(valSPAN/2)^3).*inv(H1*H2); % Inner elements of A matrix

% Calculate outer most row and column to form A matrix
b_01_0n = -sum(b_inner,1);
b_00 = -sum(b_01_0n,2);

% Outer matrix elements defined by C18 and C21
b_outer = [b_00, b_01_0n];

% Concatenate outer elements to inner matrix elements
matA = horzcat(b_01_0n',b_inner);
matA = vertcat(b_outer, matA);

%% Calculating B matrix from C30 for torsion response

j = (2./((valSPAN/2).*vecLAMBDA(i)')).*(1./((1./matGJt(i-1,1)'+1./matGJt(i,1)')));

matB = zeros(num_el);

% Lower and upper diagonal elements of the B matrix
matB_lower = diag(-j,-1);
matB_upper = diag(-j,1);

% Temp vectors to create diagonal elements
temp_j1 = [j(1:end),0];
temp_j2 = [0, j(1:end)];

% Diagonal elements of B matrix
matB_diag = diag(temp_j1+temp_j2,0);

% Final B matrix from C30
matB = matB_diag + matB_lower + matB_upper + matB;

% C matrix of eqn 43
matC = [matA, zeros(num_el); zeros(num_el), matB];

%% Calculating load vectors p and q
vecMBAR = vecMASS + (pi.*valDENSITY.*(2*vecDVEHVSPN).*vecCRD.^2)./4; % Mass including apparent mass effects
vecMEBAR = vecMASS.*vecLSM + ((pi.*valDENSITY.*(2*vecDVEHVSPN).*vecCRD.^3)./4).*(0.5 - (0.25.*vecCRD + vecLSAC)./vecCRD); % Mass moment including apparent mass effects
vecK2 = vecJT./vecMASS; % Radiation of gyration squared
vecMKBAR = vecMASS.*vecK2 +((pi.*valDENSITY.*(2*vecDVEHVSPN).*vecCRD.^4)./4).*(0.5 - (0.25.*vecCRD + vecLSAC)./vecCRD).^2 ...
    + (pi*valDENSITY.*(2*vecDVEHVSPN).*vecCRD.^4)./128;

% Eta terms from eqn A7
eta0 = -2*vecMBAR./(valDELTIME*valDELTIME);
eta1 = 5*vecMBAR./(valDELTIME*valDELTIME);
eta2 = -4*vecMBAR./(valDELTIME*valDELTIME);
eta3 = vecMBAR./(valDELTIME.*valDELTIME);

eta0prime = 2*vecMEBAR./(valDELTIME*valDELTIME) + (11/(24*valDELTIME))*valBETA.*vecCRD.*vecCRD.*(2*vecDVEHVSPN);
eta1prime = -5*vecMEBAR./(valDELTIME*valDELTIME) - (3/(4*valDELTIME))*valBETA.*vecCRD.*vecCRD.*(2*vecDVEHVSPN);
eta2prime = 4*vecMEBAR./(valDELTIME*valDELTIME) + (9/(24*valDELTIME))*valBETA.*vecCRD.*vecCRD.*(2*vecDVEHVSPN);
eta3prime = -vecMEBAR./(valDELTIME*valDELTIME) - (1/(12*valDELTIME))*valBETA.*vecCRD.*vecCRD.*(2*vecDVEHVSPN);

% Form diagonal eta matrices from Appendix A
matETA0 = diag(eta0,0);
matETA1 = diag(eta1,0);
matETA2 = diag(eta2,0);
matETA3 = diag(eta3,0);

matETA0PRIME = diag(eta0prime,0);
matETA1PRIME = diag(eta1prime,0);
matETA2PRIME = diag(eta2prime,0);
matETA3PRIME = diag(eta3prime,0);

% Nu terms from eqn A4
nu0 = 2*vecMEBAR./(valDELTIME*valDELTIME);
nu1 = -5*vecMEBAR./(valDELTIME*valDELTIME);
nu2 = 4*vecMEBAR./(valDELTIME*valDELTIME);
nu3 = -vecMEBAR./(valDELTIME*valDELTIME);

nu0prime = 2*vecMKBAR./(valDELTIME*valDELTIME) + (11/(24*valDELTIME))*valBETA.*vecCRD.*vecCRD.*vecCRD.*(2*vecDVEHVSPN).*(0.75 - (0.25.*vecCRD + vecLSAC)./vecCRD);
nu1prime = -5*vecMKBAR./(valDELTIME*valDELTIME) - (3/(4*valDELTIME))*valBETA.*vecCRD.*vecCRD.*vecCRD.*(2*vecDVEHVSPN).*(0.75 - (0.25.*vecCRD + vecLSAC)./vecCRD);
nu2prime = 4*vecMKBAR./(valDELTIME*valDELTIME) + (9/(24*valDELTIME))*valBETA.*vecCRD.*vecCRD.*vecCRD.*(2*vecDVEHVSPN).*(0.75 - (0.25.*vecCRD + vecLSAC)./vecCRD);
nu3prime = -vecMKBAR./(valDELTIME*valDELTIME) - (1/(12*valDELTIME))*valBETA.*vecCRD.*vecCRD.*vecCRD.*(2*vecDVEHVSPN).*(0.75 - (0.25.*vecCRD + vecLSAC)./vecCRD);

% Form diagonal nu matrices from Appendix A
matNU0 = diag(nu0,0);
matNU1 = diag(nu1,0);
matNU2 = diag(nu2,0);
matNU3 = diag(nu3,0);

matNU0PRIME = diag(nu0prime,0);
matNU1PRIME = diag(nu1prime,0);
matNU2PRIME = diag(nu2prime,0);
matNU3PRIME = diag(nu3prime,0);

% Form Si matrices from eqn 62
S0 = [matETA0, matETA0PRIME; matNU0, matNU0PRIME];
S1 = [matETA1, matETA1PRIME; matNU1, matNU1PRIME];
S2 = [matETA2, matETA2PRIME; matNU2, matNU2PRIME];
S3 = [matETA3, matETA3PRIME; matNU3, matNU3PRIME];

%% Assemble D matrix and Q vector and solve for deflection and twist according to eqn 65 and 66
matD = matC - S0;

% Determining initial response at first 3 timesteps
if valTIMESTEP == 1
    
    % Response at n
    matRES = (matD + S2 + 8*S3)\[vecLIFTDIST;vecMOMDIST];
    matDEF(valTIMESTEP+3,:) = matRES(1:num_el,1);
    matTWIST(valTIMESTEP+3,:) = matRES(num_el+1:end,1);
    
    % Response at n-1
    matDEF(valTIMESTEP+2,:) = zeros(1,num_el);
    matTWIST(valTIMESTEP+2,:) = zeros(1,num_el);

    % Response at n-2
    matDEF(valTIMESTEP+1,:) = -matRES(1:num_el,1);
    matTWIST(valTIMESTEP+1,:) = -matRES(num_el+1:end,1);    
    
    %Response at n-3
    matDEF(valTIMESTEP,:) = -8*matRES(1:num_el,1);
    matTWIST(valTIMESTEP,:) = -8*matRES(num_el+1:end,1);
    
else
    
    vecQ = S1*[matDEF(end,:)'; matTWIST(end,:)'] + S2*[matDEF(end-1,:)'; matTWIST(end-1,:)'] + S3*[matDEF(end-2,:)'; matTWIST(end-2,:)'] + [vecLIFTDIST; vecMOMDIST];

    matRES = matD\vecQ;

    matDEF(valTIMESTEP+3,:) = matRES(1:num_el,1);
    matTWIST(valTIMESTEP+3,:) = matRES(num_el+1:end,1);

end

end