function [matUINF,test] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valALPHA,valBETA,valGUSTSTART,matDVE,matCENTER,vecDVEHVSPN,test)

% This function modifies matUINF to model a sinusoidal gust.

% 2017/06/04 - I-80W, Omaha, NE

% Gust period
valPER = valGUSTL/valUINF;
x0 = valGUSTTIME*valDELTIME*valUINF;
x = 0;

start_loc = repmat([-valGUSTSTART*valDELTIME*valUINF,0,0],size(matCENTER,1),1); % Location (in meters) in global frame where gust starts

delx = start_loc - matCENTER; % Distance between DVE points and gust starting point

idx1 = delx(:,1) >= 0 & delx(:,1) <= valGUSTL;
idx2 = find(idx1 > 0);

tau = delx(idx2,1)./valUINF;

matUINF = repmat([valUINF*cos(valALPHA)*cos(valBETA) valUINF*sin(valBETA) valUINF*sin(valALPHA)*cos(valBETA)],size(matUINF,1),1);

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
        matUINF(idx2,3) = matUINF(idx2,3) + 0.5*valGUSTAMP*(1 - cos((2*pi*tau/(valGUSTL/valUINF))));
        test(valGUSTTIME,:) = matUINF(4,3);
    end
    
% Create gust velocity for sharp edge gust
elseif flagGUSTMODE == 3
    
    if any(idx1) > 0
        matUINF(idx2,3) = matUINF(idx2,3) + valGUSTAMP;
    end
    
else
    
    disp('No gust mode exists for entered value');
    
end

end
