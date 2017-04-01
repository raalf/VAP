function [theta_dist] = fcnSTATICTWIST(vecLIFTDIST, vecDVEHVSPN, matGJt, vecEDGEPITCH, vecCLDIST, matTWISTGLOB,...
    valTIMESTEP, valDENSITY, valUINF, valSPAN, vecLEDVES, matCENTER, vecSPANDIST, vecDVELE, vecMOMDIST)

s = 0.5*valSPAN;
% 
% matGJt = 20000*(1-vecSPANDIST/(2*s)).^2;
% [ledves, ~, ~] = find(vecDVELE > 0);
% 
% tempSPANDIST = matCENTER(ledves,2);
% 
% tempSPANDIST = repmat(tempSPANDIST', size(tempSPANDIST,1),1);
% 
% tempSPANDIST = triu(tempSPANDIST);
% 
% vecCLDIST = repmat(vecCLDIST', size(vecCLDIST,1),1);
% 
% vecCLDIST = triu(vecCLDIST);
% 
% vecCLDIST = ((vecSPANDIST(2:(end-1))' - tempSPANDIST(1,1:(end-1)))./(tempSPANDIST(2,2:end)-...
%     tempSPANDIST(1,1:(end-1)))).*(vecCLDIST(2,2:end)-vecCLDIST(1,1:(end-1))) + vecCLDIST(1,1:(end-1)); % Linear interpolation
% 
% CL_root = vecCLDIST(1,1) - (tempSPANDIST(1,1) - vecSPANDIST(1)).*(vecCLDIST(1,2)-vecCLDIST(1,1))./(tempSPANDIST(1,2)-tempSPANDIST(1,1));
% vecCLDIST = [CL_root, vecCLDIST, 0];
% 
% % vecCL_alpha = vecCLDIST./(vecEDGEPITCH+matTWISTGLOB(valTIMESTEP-1,:));
% vecCL_alpha = (0.07*(1-vecSPANDIST./s).^2)./sqrt(1-0.5*0.5);
% 
% % q = 0.5*valDENSITY*valUINF*valUINF;
% q = 0.5*1.4*50000*0.5*0.5;
% 
% % f = -0.0000809.*(vecSPANDIST).^2 + 0.0028988.*(vecSPANDIST);
% f = vecSPANDIST./s;
% 
% % f_prime = 2*-0.0000809.*(vecSPANDIST) + 0.0028988;
% f_prime = 1/s;
% 
% A = (1/s).*trapz(linspace(0,s,18),matGJt(:,1).*f_prime.*f_prime);
% 
% % vecLIFT = vecLIFTDIST(1:(end-1))'.*(2*vecDVEHVSPN(vecLEDVES));
% % vecLIFT = [vecLIFT; 0];
% 
% B = ((1.5*1.5)*0.1/s).*trapz(linspace(0,s,18),vecCL_alpha.*f.*f);
% 
% C = ((1.5*1.5)*0.1/s).*trapz(linspace(0,s,18),(vecCL_alpha*5.*f));
% 
% theta_r = q*C/(A - B);
% 
% qd = A/(s*s*B);
% 
% theta_dist = (theta_r.*f);

coeff = coeffvalues(fit(vecSPANDIST,vecMOMDIST,'poly9'));

n = 1;
for i = 10:-1:1
    
    temp_int(i) = coeff(n)*s^(i)/(i);
    coeff_new(n) = coeff(n)/i;
    n = n + 1;
    
end

C1 = sum(temp_int);

coeff_new = [coeff_new, C1];

coeff_new = coeff_new./((11:-1:1)*matGJt(1,1));

n = 1;
for i = 11:-1:1
    
    temp_theta(n,:) = coeff_new(n)*vecSPANDIST'.^(i);
    n = n + 1;
    
end

theta_dist = sum(temp_theta,1);

end