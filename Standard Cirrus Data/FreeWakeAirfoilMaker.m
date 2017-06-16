clear
clc

% AirfoilName = 'fx6617a2';
AirfoilName = 'NACA63012A';

Polars = XFOIL_Polars(AirfoilName, [20000,50000,100000,150000,200000,250000,500000,750000,1000000,2000000,3000000], -5:0.1:15);


%
fn = fieldnames(Polars);
output = [];
for n = 1:length(fn)
    Polar = Polars.(fn{n});
    output_temp = [[Polar.alpha]',[Polar.CL]',[Polar.CD]',[Polar.Re]',[Polar.CM]'];
    [~,idx] = sort([Polar.CL]');
    output_temp = output_temp(idx,:);
    output = [output;output_temp];
end

%
fileID = fopen('freewake_airfoil.dat','w');

fprintf(fileID, '%s,        free transition row= %i\n', AirfoilName, length(output(:,1)));
for n = 1:length(output(:,1));
    fprintf(fileID,'%.4f %.4f %.4f %i %.4f \n',...
        output(n,1),output(n,2),output(n,3),output(n,4),output(n,5));
end
fclose(fileID);














