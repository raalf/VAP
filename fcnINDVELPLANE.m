function [ induced_vel ] = fcnINDVELPLANE(vecROTAX, valMAXTIME, vecAZNUM, valNELE, matDVE, matVLST, matCOEFF, vecK, ...
    vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
    matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
    vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP)
% This function calculated the induced drag in a desired plane. This was
% created to validate Julia's wake interference model

% if valTIMESTEP > valMAXTIME - valAZNUM
%     b = b+1;
%     x_offset = vecROTAX(1);
%     z_offset = vecROTAX(3);
%     x = -0.6:0.01:0.6;
%     x = x + x_offset;
%     fpg1 = zeros(size(x,2),3);
%     fpg1(:,1) = x';
%     fpg1(:,3) = (fpg1(:,3)+1)*z_offset;
% 
%     % Along y line when z=x=0
%     y = -0.6:0.01:0.6;
% 
%     fpg2 = zeros(size(y,2),3);
%     fpg2(:,2) = y';
%     fpg2(:,1) = (fpg2(:,1)+1)*x_offset;
%     fpg2(:,3) = (fpg2(:,3)+1)*z_offset;
% 
%     % Calculate induced velocities
%     induced_vel1(:,:,b) = fcnINDVEL(fpg1,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%     vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%     vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%     matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%     vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
% 
%     induced_vel2(:,:,b) = fcnINDVEL(fpg2,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%     vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%     vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%     matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%     vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
% end
if valTIMESTEP > valMAXTIME - max(vecAZNUM)
    x_offset = vecROTAX(1,1);
    z_offset = vecROTAX(1,3);
%         x = -0.6:0.01:0.6;
%         x = x + x_offset;
    [x,y] = meshgrid(-0.6:0.05:0.6, -0.6:0.05:0.6);
    x = reshape(x,1,numel(x))' + x_offset;
    y = reshape(y,1,numel(y))';
    fpg = zeros(size(x,1),3);
    fpg(:,1) = x;
    fpg(:,2) = y;
    fpg(:,3) = (fpg(:,3)+1)*z_offset;

%         % Along y line when z=x=0
%         y = -0.6:0.01:0.6;
% 
%         fpg2 = zeros(size(y,2),3);
%         fpg2(:,2) = y';
%         fpg2(:,1) = (fpg2(:,1)+1)*x_offset;
%         fpg2(:,3) = (fpg2(:,3)+1)*z_offset;

    % Calculate induced velocities
    induced_vel(:,:,v-(valMAXTIME-valTIMESTEP)) = fcnINDVEL(fpg,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
    vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
    vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
    matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
    vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);

%         induced_vel2(:,:,b) = fcnINDVEL(fpg2,valNELE, matDVE, matVLST, matCOEFF, vecK, ...
%         vecDVEHVSPN, vecDVEHVCRD, vecDVEROLL, vecDVEPITCH, vecDVEYAW, ...
%         vecDVELESWP, vecDVETESWP, vecSYM, valWNELE, matWDVE, matWVLST, ...
%         matWCOEFF,vecWDVEPITCH , vecWDVEHVSPN,vecWDVEHVCRD,vecWDVEROLL, ...
%         vecWK, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, valWSIZE, valTIMESTEP);
end


end

