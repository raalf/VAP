function [ Polars ] = XFOIL_Polars( filename, ReRange, alphaRange )
%XFOIL_ Summary of this function goes here
%   Detailed explanation goes here
%%
%
% clear
% clc
% filename = 'e1230';
% ReRange = [3e6;6e6;1e7];
% alphaRange = -5:1:5;


filename = sprintf('airfoil/%s.dat',filename); % add .dat


% exePath = '/Applications/Xfoil.app/Contents/Resources/xfoil';
exePath = 'xfoil.exe';
% n2 = 0;

for j = 1:length(ReRange)
    try
        fclose all;
        delete('XFOIL_Polar_Output.dat');
    end
    
    Re = ReRange(j)
    
    fileID = fopen('XFOIL_Commands.txt','w');
    fprintf(fileID,'plop\n');
    fprintf(fileID,'g\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'load\n');
    fprintf(fileID,'%s\n', filename);
    fprintf(fileID,'oper\n');
    fprintf(fileID,'iter\n');
    fprintf(fileID,'100\n'); %any more it will probably never converge
    fprintf(fileID,'visc\n');
    fprintf(fileID,'%i\n', Re);
    fprintf(fileID,'pacc\n');
    fprintf(fileID,'%s\n', 'XFOIL_Polar_Output.dat');
    fprintf(fileID,'\n');
    for alpha = alphaRange
         
        if alpha < -5 || alpha > 14
           
            fprintf(fileID,'init\n');
        end
        fprintf(fileID,'alfa\n');
        fprintf(fileID,'%.3f\n', alpha);
    end
    fprintf(fileID,'pacc\n');
    fclose(fileID);
    
    %% RUN XFOIL
    [~,~] = system(sprintf('%s <XFOIL_Commands.txt> output.txt',exePath));
    
    %% Read XFOIL_Polar_Output.dat
    
    
    fileID = fopen('XFOIL_Polar_Output.dat','r');
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,12, 'ReturnOnError', false);
    fclose(fileID);
    matArray = [dataArray{:,1:7}];
    heading = {'alpha','CL','CD','CDp','CM','Top_Xtr','Bot_Xtr'};
    
    %%
    
%         for n = 1:length(matArray(:,1))
%             for i = 1:length(matArray(1,:))
%                 Polars.(heading{i})(n+((j-1)*length(matArray(:,1))),:) = matArray(n,i);
%                 Polars.('Re')(n+((j-1)*length(matArray(:,1))),:) = Re;
%             end
%         end
%         if any((Polars.Re) == 0) == 1
%             warning('Some runs are not converged, check Polars.Re for zeros')
%         end
            for n = 1:length(matArray(:,1))
                for i = 1:length(matArray(1,:))
                    %             Polars(n+n2).(heading{i}) = matArray(n,i);
                    %             Polars(n+n2).Re = Re;
        
                    %Polars2
                    Polars.(sprintf('Re_%i',Re'))(n).(heading{i}) = matArray(n,i);
                    Polars.(sprintf('Re_%i',Re'))(n).('Re') = Re;
                end
            end
        
        
            if length(ReRange) == 1
                Polars = Polars.(sprintf('Re_%i',Re'));
            end
        
        %     n2 = length(Polars);
        
        
        
    
    
end













