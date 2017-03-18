clc
clear

% A = dlmread('mem.txt')
% 
% plot(A./1000000)
% 
% grid on
% box on
% xlabel('Function Call #','FontSize',15);
% ylabel('FP_0 Size (Mb)','FontSize',15);
% axis tight

A = dlmread('gpumem.txt');

A(:,3) = floor(A(:,3))*3;

scatter(A(:,3), A(:,1), 'ob');
hold on
scatter(A(:,3), A(:,2), 'xk');

grid on
box on
xlabel('Array Elements','FontSize',15);
ylabel('Time','FontSize',15);
axis tight
