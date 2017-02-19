clc
clear


hFig20 = figure(20);
clf(20)
hFig20.Renderer = 'painters';

%%

fp1 = sprintf('maughmer.dat');
fp2 = sprintf('bikle.dat');
fp3 = sprintf('johnson.dat');

A = dlmread(fp1);
B = dlmread(fp2);
C = dlmread(fp3);

% % set(hFig1, 'Position', X)
% hold on
% scatter(A(:,3),A(:,2),50,'markeredgecolor','r');
% axis tight
% scatter(B(:,3),B(:,2),50,'^','markeredgecolor','b');
% scatter(C(:,3),C(:,2),50,'x','markeredgecolor',[0 0 0]);
% 
% grid on
% hold off

% hFig2 = figure(2);
% clf(2);
% % set(hFig2, 'Position', X)
% scatter(A(:,7),-A(:,8),50,'markeredgecolor','r');
% hold on
% scatter(B(:,7),-B(:,8),50,'^','markeredgecolor','b');
% 
% scatter(C(:,7),-C(:,8),50,'x','markeredgecolor',[0 0 0]);
% plot(D(:,5),-D(:,6),'LineWidth',1);
% plot(E(:,2)*3.6,-E(:,7),'--k');
% xlabel('V_i_n_f (km/h)','fontsize',15);
% ylabel('w_s_i_n_k (m/s)','fontsize',15);
% % title('Speed vs. Sink Rate for Standard Cirrus','fontsize',15);
% axis tight
% legend('Maughmer (Corrected to 334 N/m^2) [11]',...
%     'Bikle (Corrected to 334 N/m^2) [2]',...
%     'Johnson (Corrected to 334 N/m^2) [6]','FreeWake Prediction (334 N/m^2)'...
%     ,'New FW','Location', 'SouthWest');
% grid on
% hold off
% fclose all;

%%

hold on
yyaxis right

plot(A(:,1).*3.6, A(:,2)./A(:,3),'-ro'); 
plot(B(:,1).*3.6, B(:,2)./B(:,3),'--^b')
% scatter(C(:,1).*3.6, C(:,2)./C(:,3),50,'x','markeredgecolor','m'); 

load('ld30030.mat');
plot(ld30030(:,1), ld30030(:,2),'-k');

load('ld33033.mat');
plot(ld33033(:,1), ld33033(:,2),'-.k');

load('ld39039.mat');
plot(ld39039(:,1), ld39039(:,2),'--k');

load('VAP2.mat');
plot(vecVINF.*3.6, (vecCLv./vecCD),':^k');

ylim([15 45])
set(gca,'Ydir','reverse');
set(gca,'ycolor','k');

% plot([73 132.5],[38 38],'k','linewidth',1);

yyaxis left 



load('sink30030.mat');
plot(sink30030(:,1), sink30030(:,2),'-k');

load('sink33033.mat');
plot(sink33033(:,1), sink33033(:,2),'-.k');

load('sink39039.mat');
plot(sink39039(:,1), sink39039(:,2),'--k');

load('VAP2.mat');
plot(vecVINF.*3.6, (vecCD./vecCLv).*vecVINF,':^k');

hold off

grid on
grid minor
box on
axis tight

ylim([0 3]);
set(gca,'Ydir','reverse');
set(gca,'ycolor','k');

xlim([60 190])
