clc
clear

load('Results/SC.mat');

Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*valWEIGHT/(valAREA*valDENSITY);

k = 1;

for wmaxth = 2:0.1:8
    
    j = 1;
    
    for i = LDindex:size(CL)
        wclimb(j,1) = fcnMAXCLIMB(CL(i), CD(i), Rrecip, wmaxth, WSroh);
        j = j + 1;
    end
    
    [wclimbMAX, indexWC] = max(wclimb);
    
    for i = 1:size(CL)
        V(i,1) = (Vcruise(i)*wclimbMAX)/(wglide(i)+wclimbMAX);
    end
    
    [VxcMAX, cruiseIndex] = max(V);
    VxcSC(k,:) = [wmaxth VxcMAX];
    k = k + 1;
    
end

clearvars -except VxcSC

%%

load('Results/SD1.mat');

Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*valWEIGHT/(valAREA*valDENSITY);

k = 1;

for wmaxth = 2:0.1:8
    
    j = 1;
    
    for i = LDindex:size(CL)
        wclimb(j,1) = fcnMAXCLIMB(CL(i), CD(i), Rrecip, wmaxth, WSroh);
        j = j + 1;
    end
    
    [wclimbMAX, indexWC] = max(wclimb);
    
    for i = 1:size(CL)
        V(i,1) = (Vcruise(i)*wclimbMAX)/(wglide(i)+wclimbMAX);
    end
    
    [VxcMAX, cruiseIndex] = max(V);
    VxcSD1(k,:) = [wmaxth VxcMAX];
    k = k + 1;
    
end

hFig14 = figure(14);
clf(14)

plot(VxcSC(:,1), 100.*(VxcSD1(:,2) - VxcSC(:,2))./VxcSC(:,2),'-k')
box on
grid on
axis tight

xlabel('Thermal Core Strength (m/s)','FontSize',15);
ylabel('% Increase in Cross-Country Speed','FontSize',15);
print(hFig14,'split1','-deps');

%% 

clc
clear

A = dlmread('optihistory2.txt');
B = A(A(:,17) > 7.45,:);

hFig15 = figure(15);
clf(15);

scatter(B(:,6),1./B(:,1))

box on
grid on
axis tight

xlabel('Projected winglet span (m)','FontSize',15);
ylabel('V_x_c in weak thermals (m/s)','FontSize',15);

print(hFig15,'span','-deps');






















