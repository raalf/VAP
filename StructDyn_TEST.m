clear,clc

valNSELE = 30;
valMAXTIME = 600000;
Len = 6;
omega = 1*2*pi;

valDY = Len/(valNSELE-1);

valSTRUCTDELTIME = 0.00001;

matEIx(:,1) = 9772210*ones(valNSELE,1);
matEIx(:,2) = zeros(valNSELE,1);
matEIx(:,3) = zeros(valNSELE,1);

matGJt(:,1) = 985810*ones(valNSELE,1);
matGJt(:,2) = zeros(valNSELE,1);

vecJT = 8.6505832*ones(valNSELE,1);

vecLM = 35.75121*ones(valNSELE,1);
vecLSM = 0.18*ones(valNSELE,1);

vecDEF = zeros(1,valNSELE+3);
vecTWIST = zeros(1,valNSELE+3);
vecSLOPE = zeros(1,valNSELE-1);

matDEF = zeros(valMAXTIME,valNSELE+3);
matTWIST = zeros(valMAXTIME,valNSELE+3);
%% Beam boundary conditions

matDEF(1:2,:) = zeros(2,valNSELE+3);
matTWIST(1:2,:) = zeros(2,valNSELE+3);

for valSTRUCTTIME = 3:valMAXTIME

vecDEF(2) = 0; % Zero deflection at root BC
vecTWIST(2) = 0; % Zero twist at root BC

vecLIFTDIST = 100*sin(omega*valSTRUCTTIME*valSTRUCTDELTIME)*ones(1,valNSELE);
vecMOMDIST = 100*sin(omega*valSTRUCTTIME*valSTRUCTDELTIME)*ones(1,valNSELE);

% Assemble load matrix
matLOAD = [vecLIFTDIST' - vecLM.*9.81, vecMOMDIST' - vecLM.*vecLSM.*9.81];
% matLOAD(end,:) = [0,0];

for yy = 3:(valNSELE+1)

    %% Geometric property assembly

    % Assemble mass matrix
    matMASS = [vecLM(yy-1), -vecLM(yy-1).*vecLSM(yy-1); -vecLM(yy-1).*vecLSM(yy-1), vecJT(yy-1)];

    % Assemble stiffness matrices
    matK_1 = [matEIx(yy-1,3), 0; 0, 0];
    matK_2 = [matEIx(yy-1,2), 0; 0, -matGJt(yy-1,2)]; 
    matK_3 = [matEIx(yy-1,1), 0; 0, -matGJt(yy-1,1)];
    matB = [0 0; 0 0];

    %% Finite difference relations for partial derivatives

    % Finite difference relations for partial derivatives w.r.t
    % time
    valUDOT = (matDEF(valSTRUCTTIME-1,yy) - matDEF(valSTRUCTTIME - 2, yy))./valSTRUCTDELTIME;
    valTDOT = (matTWIST(valSTRUCTTIME-1,yy) - matTWIST(valSTRUCTTIME - 2,yy))./valSTRUCTDELTIME;
        
    % Finite difference relations for partial derivative of deflection w.r.t Y
    valU_yy = (matDEF(valSTRUCTTIME-1,yy+1) - 2*matDEF(valSTRUCTTIME-1,yy) + matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^2;
    valU_yyy = (matDEF(valSTRUCTTIME-1,yy+2) - 3*matDEF(valSTRUCTTIME-1,yy+1) + 3*matDEF(valSTRUCTTIME-1,yy)- ...
        matDEF(valSTRUCTTIME-1,yy-1))/(valDY)^3;
    valU_yyyy = (matDEF(valSTRUCTTIME-1,yy+2) - 4*matDEF(valSTRUCTTIME-1,yy+1) + 6*matDEF(valSTRUCTTIME-1,yy) - ...
        4*matDEF(valSTRUCTTIME-1,yy-1) + matDEF(valSTRUCTTIME-1,yy-2))/(valDY)^4;

    % Finite difference relations for partial derivative of twist w.r.t Y
    valTHETA_y = (matTWIST(valSTRUCTTIME-1,yy+1) - matTWIST(valSTRUCTTIME-1,yy-1))/(2*valDY);
    valTHETA_yy = (matTWIST(valSTRUCTTIME-1,yy+1) - 2*matTWIST(valSTRUCTTIME-1,yy) + matTWIST(valSTRUCTTIME-1,yy-1))/(valDY^2);

    %% Solve matrix equation

    % Temp variable with the wing deflection and twist stored as a matrix. The
    % first row is the deflection, w/ each column as a spanwise station. The
    % second row is the twist, w/ each column as a spanwise station.

    tempTWISTBEND = 2.*[matDEF(valSTRUCTTIME-1,yy); matTWIST(valSTRUCTTIME-1,yy)] - [matDEF(valSTRUCTTIME-2,yy); matTWIST(valSTRUCTTIME-2,yy)] ...
        + (valSTRUCTDELTIME^2).*inv(matMASS)*([matLOAD(yy-1,1); matLOAD(yy-1,2)] - matK_1*[valU_yy; 0] - matK_2*[valU_yyy; valTHETA_y] -...
        matK_3*[valU_yyyy; valTHETA_yy] - matB*[valUDOT; valTDOT]);

    % Output result of deflection and twist to separate vectors
    vecDEF(yy) = tempTWISTBEND(1,:);
    vecTWIST(yy) = tempTWISTBEND(2,:);

    % Calculate angle between DVE and horizontal based on
    % deflection
    vecSLOPE(yy-2) = asin((vecDEF(yy)-vecDEF(yy-1))/(valDY));

end

vecDEF(valNSELE+2) = 2*vecDEF(valNSELE+1)...
    -vecDEF(valNSELE); % BC for deflection one element beyond wing (positive span direction)

vecDEF(valNSELE+3) = 3*vecDEF(valNSELE+1)...
    -2*vecDEF(valNSELE); % BC for deflection two elements beyond wing (positive span direction)

vecDEF(1) = vecDEF(3); % BC for deflection one element beyond root (negative span direction)

vecTWIST(valNSELE+2) = vecTWIST(valNSELE); % BC for twist one element beyond wing tip (positive span direction)

matDEF(valSTRUCTTIME,:) = vecDEF;
matTWIST(valSTRUCTTIME,:) = vecTWIST;

% Spanwise deflection and twist wrt structural timestep
vecDEF = matDEF(end,:);
vecTWIST = matTWIST(end,:);

end

res = fft(matTWIST(:,end-1));
Fs = 1/valSTRUCTDELTIME;
L = valMAXTIME;
P2 = abs(res/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure;
plot(f,P1)
xlim([0 0.5e4])
xlim([0 100])

fn1 = (1.875^2/(2*pi*Len^2))*sqrt(matEIx(1,1)/(vecLM(1)*Len))
fn2 = (4.694^2/(2*pi*Len^2))*sqrt(matEIx(1,1)/(vecLM(1)*Len))
fn3 = (7.8539^2/(2*pi*Len^2))*sqrt(matEIx(1,1)/(vecLM(1)*Len))

fnt1 = (pi/(2*Len)).*sqrt(matGJt(1,1)/vecJT(1))/(2*pi)
fnt2 = (3*pi/(2*Len)).*sqrt(matGJt(1,1)/vecJT(1))/(2*pi)
fnt3 = (5*pi/(2*Len)).*sqrt(matGJt(1,1)/vecJT(1))/(2*pi)

figure;
plot((1:valMAXTIME).*valSTRUCTDELTIME,matTWIST(:,end-1))
% valNELE = size(vecLIFTDIST,2);
% 
% valDY = 0.01;
% valNELE = 500;
% 
% vecLOAD = -5*ones(1,valNELE)*9.81;
% matEIx = [];
% matEIx(:,1) = 20000*ones(valNELE,1);
% 
% S(valNELE) = 0;
% M(valNELE) = 0;
% 
% for yy = (valNELE-1):-1:1
%     S(yy) = S(yy+1) - ((vecLOAD(yy+1)+vecLOAD(yy))/2)*valDY;
%     M(yy) = M(yy+1) - ((S(yy+1)+S(yy))/2)*valDY;
% end
% 
% theta(1) = 0;
% w(1) = 0;
% 
% for yy = 2:valNELE
%     theta(yy) = theta(yy-1) + 0.5*((M(yy)+M(yy-1))/matEIx(yy,1))*valDY;
%     w(yy) = w(yy-1) + ((theta(yy)+theta(yy-1))/2)*valDY;
%     vecSLOPE(yy) = asin((w(yy)-w(yy-1))/(valDY));
% end