clc
clear


hFig20 = figure(20);
clf(20)

hFig20.Renderer = 'painters';

hold on
yyaxis right

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

