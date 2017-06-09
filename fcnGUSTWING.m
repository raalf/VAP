function [matUINF] = fcnGUSTWING(matUINF,valGUSTAMP,valGUSTL,flagGUSTMODE,valDELTIME,valGUSTTIME,valUINF)

% This function modifies matUINF to model a sinusoidal gust.

% 2017/06/04 - I-80W, Omaha, NE

% Gust period
valPER = valGUSTL/valUINF;

% Create gust velocity for sine wave gust
if flagGUSTMODE == 1
    
    valGUSTVEL = valGUSTAMP*sin((2*pi/valPER)*valDELTIME*valGUSTTIME);

% Create gust velocity for 1-cosine gust
elseif flagGUSTMODE == 2
    
    valGUSTVEL = 0.5*valGUSTAMP*(1 - cos((2*pi*valUINF*valDELTIME*valGUSTTIME)/valGUSTL));
    
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

