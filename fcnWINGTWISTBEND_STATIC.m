function [vecDEF_res] = fcnWINGTWISTBEND_STATIC(matEIx, matGJt, vecLIFTDIST, vecMOMDIST, vecLM,vecDVEHVSPN)

valSDELTIME = 0.0015;
valDY = sum(2*vecDVEHVSPN,1)/length(vecDVEHVSPN);

vecLOAD = vecLIFTDIST' - vecLM*9.81;
valNELE = size(vecLIFTDIST,2)+3;

vecDEF = zeros(2,valNELE);

for valSTIME = 3:20000

vecDEF(valSTIME,2) = 0;

for yy = 3:(valNELE-2)
    vecDEF(valSTIME,yy) = 2*vecDEF(valSTIME-1,yy) - vecDEF(valSTIME-2,yy) -...
        (matEIx(yy-1,1)*valSDELTIME^2/(vecLM(yy-1)*valDY^4))*(vecDEF(valSTIME-1,yy+2) -...
        4*vecDEF(valSTIME-1,yy+1) + 6*vecDEF(valSTIME-1,yy) - 4*vecDEF(valSTIME-1,yy-1) + vecDEF(valSTIME-1,yy-2))+(vecLOAD(yy-1)*valSDELTIME^2/vecLM(yy-1));
end
vecDEF(valSTIME,valNELE-1) = 2*vecDEF(valSTIME,valNELE-2) - vecDEF(valSTIME,valNELE-3);
vecDEF(valSTIME, valNELE) = 3*vecDEF(valSTIME,valNELE-2) - 2*vecDEF(valSTIME,valNELE-3);
vecDEF(valSTIME,1) = vecDEF(valSTIME,3);
end


vecDEF_res = vecDEF(end,:);

end

