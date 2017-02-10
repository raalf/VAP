function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST, matNPVLST, matUINF] = fcnSTIFFWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST)

% This function moves the wing in the freestream direction and calculates
% the new wake elements, asssuming no bending of the wing.

matCENTER_old = matCENTER;

[matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE, matNTVLST] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNTVLST);

matNPVLST = matNTVLST;

[matUINF] = fcnFLEXUINF(matCENTER_old, matCENTER, valDELTIME);

end