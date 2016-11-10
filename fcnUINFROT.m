function [vecUINF] = fcnUINFROT(valALPHAR)

% This function defines the direction of the inflow velocity

% INPUT:
%   valALPHA - Angle of attack for this move (radians)
% OUTPUT:
%   vecUINF - vector of freestream velocity

uinf = 1;

vecUINF = [uinf*cos(valALPHAR) 0 uinf*sin(valALPHAR)];

end

