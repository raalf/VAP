function [matUINF, gust_vel_old] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valGUSTSTART,matCENTER,gust_vel_old,test)

% This function modifies matUINF to model a sinusoidal gust.

% 2017/06/04 - I-80W, Omaha, NE

% Gust period
valPER = valGUSTL/valUINF;

start_loc = repmat([-valGUSTSTART*valDELTIME*valUINF,0,0],size(matCENTER,1),1); % Location (in meters) in global frame where gust starts

delx = start_loc - matCENTER; % Distance between DVE points and gust starting point

idx1 = delx(:,1) >= 0 & delx(:,1) <= valGUSTL;
idx2 = find(idx1 > 0);

tau = delx(idx2,1)./valUINF;

% matUINF = repmat([valUINF*cos(valALPHA)*cos(valBETA) valUINF*sin(valBETA) valUINF*sin(valALPHA)*cos(valBETA)],size(matUINF,1),1);

% Create gust velocity for sine wave gust
if flagGUSTMODE == 1
    
    if valPER >=  valGUSTTIME*valDELTIME
        valGUSTVEL = valGUSTAMP*sin((2*pi/valPER)*valDELTIME*valGUSTTIME);
    else
        valGUSTVEL = 0;
    end

% Create gust velocity for 1-cosine gust
elseif flagGUSTMODE == 2
    
    if any(idx1) > 0
        gust_vel = 0.5*valGUSTAMP*(1 - cos((2*pi*tau/(valGUSTL/valUINF))));
        matUINF(idx2,3) = matUINF(idx2,3) + (gust_vel - gust_vel_old(idx2,1));
        gust_vel_old(idx2) = gust_vel;
    end
    
% Create gust velocity for sharp edge gust
elseif flagGUSTMODE == 3
    
    if any(idx1) > 0
        gust_vel = valGUSTAMP;
        matUINF(idx2,3) = matUINF(idx2,3) + (gust_vel);
        gust_vel_old(idx2) = gust_vel;
    end
    
elseif flagGUSTMODE == 0
    
    gust_vel = 0;
    
else
    
    disp('No gust mode exists for entered value')
    
end

end
