function [w_wake] = fcnWDVEVEL(fpg, valWNELE, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM, valWSIZE, valTIMESTEP)

% vecWDVEHVSPN, vecWDVEHVCRD, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, ...
% vecWDVEMCSWP, vecWDVETESWP, matWDVENORM, matWVLST, matWDVE, valWNELE, ...
% matWCENTER, matWCOEFF, vecWK

% List of wDVEs we are influencing from (each one for each of the fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valWNELE,len,1),[],1);

fpg = repmat(fpg,valWNELE,1);

% DVE type 1 is a regular wake DVE
dvetype = ones(length(dvenum),1);

% Newest row of wake DVEs have a filament at the leading edge
newest_row = [((valWNELE-valWSIZE)+1):1:valWNELE]';
dvetype(ismember(dvenum, newest_row)) = 2;

% Oldest row of wake DVEs are semi-infinite
oldest_row = [1:valWSIZE]';

if valTIMESTEP == 1
    dvetype(ismember(dvenum, oldest_row)) = -3;
else
    dvetype(ismember(dvenum, oldest_row)) = 3;
end


[w_ind] = fcnDVEVEL(dvenum, fpg, dvetype, matWDVE, matWVLST, matWCOEFF, vecWK, vecWDVEHVSPN, vecWDVEROLL, vecWDVEPITCH, vecWDVEYAW, vecWDVELESWP, vecWDVETESWP, vecSYM);

w_wake = reshape(sum(reshape(w_ind', len*3, [])',1),3,[])';
end
