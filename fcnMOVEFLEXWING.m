function [matVLST, matCENTER, matNEWWAKE, matNPNEWWAKE] = fcnMOVEFLEXWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matDVE, vecDVETE, matNPVLST, matDEFGLOB, matTWISTGLOB, valVINF)

matDEFVEL = (matDEFGLOB./valDELTIME)./valVINF;

uinf = 1;

uinf = repmat([uinf*cos(valALPHA)*cos(valBETA) uinf*sin(valBETA) uinf*sin(valALPHA)*cos(valBETA)],sum(vecN,1)+1,1);

translation = valDELTIME.*uinf;

end