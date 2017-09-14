function [matUINF] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF,valALPHA,valBETA)

% This function modifies matUINF to model a sinusoidal gust.

% 2017/06/04 - I-80W, Omaha, NE

% Gust period
valPER = valGUSTL/valUINF;

matUINF = repmat([valUINF*cos(valALPHA)*cos(valBETA) valUINF*sin(valBETA) valUINF*sin(valALPHA)*cos(valBETA)],size(matUINF,1),1);

% Create gust velocity for sine wave gust
if flagGUSTMODE == 1
    
    if valPER >=  valGUSTTIME*valDELTIME
        valGUSTVEL = valGUSTAMP*sin((2*pi/valPER)*valDELTIME*valGUSTTIME);
    else
        valGUSTVEL= 0;
    end

% Create gust velocity for 1-cosine gust
elseif flagGUSTMODE == 2
    
    if valPER >= valGUSTTIME*valDELTIME
        valGUSTVEL = 0.5*valGUSTAMP*(1 - cos((2*pi*valUINF*valDELTIME*valGUSTTIME)/valGUSTL));
    else
        valGUSTVEL = 0;
    end
    
% Create gust velocity for sharp edge gust
elseif flagGUSTMODE == 3
    
    if valPER >= valGUSTTIME*valDELTIME
        valGUSTVEL = valGUSTAMP;
    else
        valGUSTVEL = 0;
    end
    
else
    
    valGUSTVEL = 0;
    
end

% Add gust velocity to matUINF
matUINF = matUINF + repmat([0,0,valGUSTVEL],size(matUINF,1),1);

