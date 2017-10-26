function [valSDELTIME, valNSELE, valSTIFFSTEPS, flagSTATIC, vecEIxCOEFF, vecGJtCOEFF, vecEACOEFF, vecCGCOEFF, vecJTCOEFF, vecLMCOEFF] = fcnSTRUCTREAD(strSTRUCT_INPUT)

% This function reads in the structure input file and stores each parameter
% in a vector. The vectors are all 1 x 3 with col 1 = A, col 2 = B, and col
% 3 = C, where A, B, and C make up the function y  = A*x^2 + B*x + C
%
% INPUT
% strSTRUCT_INPUT - String pointing to path of structure input file
%
% OUTPUT
% vec_EIx - 1 x 3 vector containing A, B, and C coefficients of quadratic
% bending stiffness distribution
% vec_GJt - 1 x 3 vector containing A, B, and C coefficients of quadratic
% torsional stiffness distribution
% vec_EA - 1 x 3 vector containing A, B, and C coefficients of quadratic
% elastic axis location distribution
% vec_CG - 1 x 3 vector containing A, B, and C coefficients of quadratic
% center of gravity location distribution
% vec_LM - 1 x 3 vector containing A, B, and C coefficients of quadratic
% linear mass distribution

% Preallocate variables
vecEIxCOEFF = zeros(1,3);
vecGJtCOEFF = zeros(1,3);
vecEACOEFF = zeros(1,3);
vecCGCOEFF = zeros(1,3);
vecJTCOEFF = zeros(1,3);
vecLMCOEFF = zeros(1,3);

% Open structure input file
fp = fopen(strSTRUCT_INPUT);

% Read structural timestep size 
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valSDELTIME = fscanf(fp,'%lf',1);

% Read number of structural grid points
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valNSELE = fscanf(fp,'%d',1);

% Read number of timesteps to run without deflecting the wing
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
valSTIFFSTEPS = fscanf(fp,'%d',1);

% Read in flag to denote static or dynamic aeroelasticity
ch = fscanf(fp,'%c',1);
while(ch~='=');
    ch = fscanf(fp,'%c',1);
end
flagSTATIC = fscanf(fp,'%d',1);

% Read A,B,C values for bending stiffness distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecEIxCOEFF(i) = fscanf(fp,'%lf',1);

end

% Read A,B,C values for torsional stiffness distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecGJtCOEFF(i) = fscanf(fp,'%lf',1);

end

% Read A,B,C values for elastic axis distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecEACOEFF(i) = fscanf(fp,'%lf',1);

end

% Read A,B,C values for center of gravity distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecCGCOEFF(i) = fscanf(fp,'%lf',1);

end

% Read A,B,C values for torsion constant distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecJTCOEFF(i) = fscanf(fp,'%lf',1);

end

% Read A,B,C values for linear mass distribution
for i = 1:3
    ch = fscanf(fp,'%c',1);
    while(ch~='=')
        ch = fscanf(fp,'%c',1);
    end
    vecLMCOEFF(i) = fscanf(fp,'%lf',1);

end

fclose(fp);
