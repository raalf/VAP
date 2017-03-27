function [matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME)

matUINF = ((matCENTER_old - matCENTER)./valDELTIME);
% matUINF = [(matCENTER_old(:,1) - matCENTER(:,1))./valDELTIME,(matCENTER(:,2) - matCENTER_old(:,2))./valDELTIME,(matCENTER(:,3) - matCENTER_old(:,3))./valDELTIME,];

% matUINF(:,3) = 0.01*matUINF(:,3);

end