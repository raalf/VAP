function [vecR] = fcnRWING(valNELE, valTIMESTEP, matCENTER, matDVENORM, vecUINF, valWNELE, matWDVE, ...
    matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
    vecWDVETESWP, vecSYM, valWSIZE, flagSTEADY)
% Resultant
% Kinematic resultant is the freestream (and wake-induced velocities summed) dotted with the
% norm of the point we are influencing on, multiplied by 4*pi

vecR = zeros(valNELE*3,1);

    len = length(matCENTER(:,1));

if valTIMESTEP < 1;
    % Flow tangency at control points goes at the bottom of the resultant
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1), matDVENORM,2);    
else
    [w_wake] = fcnWDVEVEL(matCENTER, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEHVCRD,vecWDVEROLL, ...
        vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, valTIMESTEP, flagSTEADY);

    % Including the wake-induced velocities,
    vecR(end-(len-1):end) = (4*pi).*dot(repmat(vecUINF,len,1)+w_wake, matDVENORM,2);  

end

end

