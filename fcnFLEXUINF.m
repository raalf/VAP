function [matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME)

matUINF = ((matCENTER_old - matCENTER)./valDELTIME);

end