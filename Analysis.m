clc
clear

load('Results/14AOA/SC.mat');

Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*valWEIGHT/(valAREA*valDENSITY);

k = 1;

therm = [2:0.25:8];
alphas = [2:0.5:12];

for wmaxth = therm 
    
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

hFig16 = figure(16);
clf(16)
plot(CDfit(alphas), CLfit(alphas),'-k')
box on
grid on
axis tight

xlabel('Drag Coefficient','FontSize',15);
ylabel('Lift Coefficient','FontSize',15);

LDfitSC = LDfit;
CDfitSC = CDfit;
CdifitSC = Cdifit;
CLfitSC = CLfit;
VinffitSC = Vinffit;


%%

load('Results/14AOA/D1.mat');

Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*valWEIGHT/(valAREA*valDENSITY);

k = 1;

for wmaxth = therm 
    
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

plot(VxcSC(:,1), 100.*(VxcSD1(:,2) - VxcSC(:,2))./VxcSC(:,2),'-.r')
box on
grid on
axis tight

xlabel('Thermal Core Strength (m/s)','FontSize',15);
ylabel('% Increase in Cross-Country Speed','FontSize',15);
% print(hFig14,'split1','-deps');
%--------------------------------------------------------------------------------------------------
hFig15 = figure(15);
clf(15)
plot(Vinffit(alphas), 100.*((Vinffit(alphas)./LDfit(alphas)) - (VinffitSC(alphas)./LDfitSC(alphas)))./(VinffitSC(alphas)./LDfitSC(alphas)),'-.r')
box on
grid on
axis tight

xlabel('Airspeed (m/s)','FontSize',15);
ylabel('% Change in Sink Rate','FontSize',15);
%--------------------------------------------------------------------------------------------------

figure(hFig16);
hold on
plot(CDfit(alphas), CLfit(alphas),'-.r')
hold off

%--------------------------------------------------------------------------------------------------
hFig17 = figure(17);
clf(17)
plot(Vinffit(alphas), 100.*((Cdifit(alphas) - CdifitSC(alphas))./CdifitSC(alphas)),'-.r')
box on
grid on
axis tight

xlabel('Airspeed (m/s)','FontSize',15);
ylabel('% Change in Induced Drag','FontSize',15);

%--------------------------------------------------------------------------------------------------
hFig18 = figure(18);
clf(18)
plot(Vinffit(alphas), 100.*(((CDfit(alphas) - Cdifit(alphas)) - (CDfitSC(alphas) - CdifitSC(alphas)))./(CDfitSC(alphas) - CdifitSC(alphas))),'-.r')
box on
grid on
axis tight

xlabel('Airspeed (m/s)','FontSize',15);
ylabel('% Change in Profile Drag','FontSize',15);

%--------------------------------------------------------------------------------------------------
hFig19 = figure(19);
clf(19)
plot(Vinffit(alphas), 100.*((CDfit(alphas) - CDfitSC(alphas))./CDfitSC(alphas)),'-.r')
box on
grid on
axis tight

xlabel('Airspeed (m/s)','FontSize',15);
ylabel('% Change in Total Drag','FontSize',15);

%--------------------------------------------------------------------------------------------------
hFig20 = figure(20);
clf(20)
plot(alphas, 100.*((CDfit(alphas).*Vinffit(alphas) - CDfitSC(alphas).*VinffitSC(alphas))./(CDfitSC(alphas).*VinffitSC(alphas))),'-.r')
box on
grid on
axis tight

xlabel('Angle of Attack (Degrees)','FontSize',15);
ylabel('% Change in Power Required','FontSize',15);

%%
load('Results/14AOA/B1.mat');

Rthermal = 150;
Rrecip = 1/Rthermal;
WSroh = 2*valWEIGHT/(valAREA*valDENSITY);

k = 1;

for wmaxth = therm 
    
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
    VxcB1(k,:) = [wmaxth VxcMAX];
    k = k + 1;
    
end


wh = [5 5];  

figure(hFig14);
hold on
plot(VxcSC(:,1), 100.*(VxcB1(:,2) - VxcSC(:,2))./VxcSC(:,2),'--b')
legend('Split Winglet','Blended Winglet','Location','NorthEast');
saveFig2Latex(hFig14, 'res-vxc.pdf', wh);
hold off

figure(hFig15);
hold on
plot(Vinffit(alphas), 100.*((Vinffit(alphas)./LDfit(alphas)) - (VinffitSC(alphas)./LDfitSC(alphas)))./(VinffitSC(alphas)./LDfitSC(alphas)),'--b')
legend('Split Winglet','Blended Winglet','Location','NorthEast');
saveFig2Latex(hFig15, 'res-wsink.pdf', wh);
hold off

figure(hFig16);
hold on
plot(CDfit(alphas), CLfit(alphas),'--b')
legend('Base Aircraft','With Split Winglet','With Blended Winglet','Location','SouthEast');
saveFig2Latex(hFig16, 'res-clcd.pdf', wh);
hold off

figure(hFig17);
hold on
plot(Vinffit(alphas), 100.*((Cdifit(alphas) - CdifitSC(alphas))./CdifitSC(alphas)),'--b')
legend('Split Winglet','Blended Winglet','Location','SouthWest');
saveFig2Latex(hFig17, 'res-cdi.pdf', wh);
hold off

figure(hFig18);
hold on
plot(Vinffit(alphas), 100.*(((CDfit(alphas) - Cdifit(alphas)) - (CDfitSC(alphas) - CdifitSC(alphas)))./(CDfitSC(alphas) - CdifitSC(alphas))),'--b')
legend('Split Winglet','Blended Winglet','Location','NorthEast');
saveFig2Latex(hFig18, 'res-cdp.pdf', wh);
hold off

figure(hFig19);
hold on
plot(Vinffit(alphas), 100.*((CDfit(alphas) - CDfitSC(alphas))./CDfitSC(alphas)),'--b')
legend('Split Winglet','Blended Winglet','Location','NorthEast');
saveFig2Latex(hFig19, 'res-cd.pdf', wh);
hold off

figure(hFig20);
hold on
plot(alphas, 100.*((CDfit(alphas).*Vinffit(alphas) - CDfitSC(alphas).*VinffitSC(alphas))./(CDfitSC(alphas).*VinffitSC(alphas))),'--b')
legend('Split Winglet','Blended Winglet','Location','NorthWest');
saveFig2Latex(hFig20, 'res-preq.pdf', wh);
hold off

%% 
% 
% clc
% clear
% 
% A = dlmread('optihistory2.txt');
% B = A(A(:,17) > 7.45,:);
% 
% hFig15 = figure(15);
% clf(15);
% 
% scatter(B(:,6),1./B(:,1))
% 
% box on
% grid on
% axis tight
% 
% xlabel('Projected winglet span (m)','FontSize',15);
% ylabel('V_x_c in weak thermals (m/s)','FontSize',15);
% 
% print(hFig15,'span','-deps');






















