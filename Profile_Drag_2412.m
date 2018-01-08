clear cd cd_t sectionCL
load naca2412Mike.dat

end_time = valMAXTIME-valGUSTSTART;

for t=1:end_time
    
    sectionCL(t,1) = 0.75*vecCL(t-1+valGUSTSTART)./12;
    cd_t(t,1) = linterp(naca2412Mike(:,2),naca2412Mike(:,3),sectionCL(t,1));
    
end

cd = mean(cd_t);
cdtotal0 = cd + vecCDI(valGUSTSTART);
cdi = mean(vecCDI(valGUSTSTART:end));
cdtotal = cd + cdi;

cdtotal_t = cd_t + vecCDI(valGUSTSTART:end-1);
cdtotal0_t = repmat(cdtotal0,size(cdtotal_t,1),1);
integrand_t = (cdtotal_t-cdtotal0_t)./cdtotal0_t;

reduction = (cdtotal-cdtotal0)/cdtotal0