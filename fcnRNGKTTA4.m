function [Z, Zdot] = fcnRNGKTTA4(valDELTIME,vecLEDVES,matEIx,vecLIFTDIST,Z,Zdot)

Faero = vecLIFTDIST;

% Number of spanwise elements
nsele = size(vecLEDVES,1);

num_DOF = nsele+1; % Number of degrees of freedom
h = valDELTIME; % Set timestep size
damp = 0.01; % Structural damping is 1% of stiffness

M = vecLM;

K = [matEIx, 0];

if valTIMESTEP == 1
    
    Z(1,1:num_DOF) = 0;
    Zdot(1,1:num_DOF) = 0;
    
else
     
    % Computing X1,Y1
    X(1,1:num_DOF) = Z(valTIMESTEP-1,1:num_DOF);
    Y(1,1:num_DOF) = Zdot(valTIMESTEP-1,1:num_DOF);
    
    % Computing F1
    tempS = 0;
    
    tempS = tempS + K(1:num_DOF,:).*X(1,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*X(1,(nsele+1):(2*nsele));
    tempS = tempS + (K(1:num_DOF,:).*Y(1,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*Y(1,(nsele+1):(2*nsele))).*damp;
    
    F(1,1:num_DOF) = (Faero(1:num_DOF) - tempS)./M(1:num_DOF);
    
    % Computing X2,Y2
    X(2,1:num_DOF) = X(1,:) + Y(1,:)*0.5*h;
    Y(2,1:num_DOF) = Y(1,:) + F(1,:)*0.5*h;
    
    % Computing F2
    tempS = 0;
    
    tempS = tempS + K(1:num_DOF,:).*X(2,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*X(2,(nsele+1):(2*nsele));
    tempS = tempS + (K(1:num_DOF,:).*Y(2,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*Y(2,(nsele+1):(2*nsele))).*damp;
    
    F(2,1:num_DOF) = (Faero(1:num_DOF) - tempS)./M(1:num_DOF);
    
    % Computing X3,Y3
    X(3,1:num_DOF) = X(1,:) + Y(2,:)*0.5*h;
    Y(3,1:num_DOF) = Y(1,:) + F(2,:)*0.5*h;
    
    % Computing F3
    tempS = 0;
    
    tempS = tempS + K(1:num_DOF,:).*X(3,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*X(3,(nsele+1):(2*nsele));
    tempS = tempS + (K(1:num_DOF,:).*Y(3,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*Y(3,(nsele+1):(2*nsele))).*damp;
    
    F(3,1:num_DOF) = (Faero(1:num_DOF) - tempS)./M(1:num_DOF);
    
    % Computing X4,Y4
    X(4,1:num_DOF) = X(1,:) + Y(3,:)*0.5*h;
    Y(4,1:num_DOF) = Y(1,:) + F(3,:)*0.5*h;
    
    % Computing F4
    tempS = 0;
    
    tempS = tempS + K(1:num_DOF,:).*X(4,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*X(1:num_DOF,(nsele+1):(2*nsele));
    tempS = tempS + (K(1:num_DOF,:).*Y(1,1:nsele) + K(1:num_DOF,(nsele+1):(2*nsele)).*Y(1:num_DOF,(nsele+1):(2*nsele))).*damp;
    
    F(4,1:num_DOF) = (Faero(1:num_DOF) - tempS)./M(1:num_DOF);
    
    % Computing Z and Zdot at current timestep
    Z(valTIMESTEP,:) = Z(valTIMESTEP - 1,:) + (h/6)*(Y(1,:) + 2*Y(2,:) + 2*Y(3,:) + Y(4,:));
    Zdot(valTIMESTEP,:) = Zdot(valTIMESTEP-1,:) + (h/6)*(F(1,:) + 2*F(2,:) + 2*F(3,:) + F(4,:));

end

end