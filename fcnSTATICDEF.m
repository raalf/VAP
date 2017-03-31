function [w,vecSLOPE] = fcnSTATICDEF(vecLIFTDIST,vecLM,matEIx,vecDVEHVSPN)

vecLOAD = vecLIFTDIST' - vecLM*9.81;
vecLOAD(end) = 0;
valDY = sum(2*vecDVEHVSPN,1)/length(vecDVEHVSPN);
valNELE = size(vecLIFTDIST,2);

S(valNELE) = 0;
M(valNELE) = 0;

for yy = (valNELE-1):-1:1
    S(yy) = S(yy+1) - ((vecLOAD(yy+1)+vecLOAD(yy))/2)*valDY;
    M(yy) = M(yy+1) - ((S(yy+1)+S(yy))/2)*valDY;
end

theta(1) = 0;
w(1) = 0;

for yy = 2:valNELE
    theta(yy) = theta(yy-1) + 0.5*((M(yy)+M(yy-1))/matEIx(yy))*valDY;
    w(yy) = w(yy-1) + ((theta(yy)+theta(yy-1))/2)*valDY;
    vecSLOPE(yy) = asin((w(yy)-w(yy-1))/(valDY));
end

end