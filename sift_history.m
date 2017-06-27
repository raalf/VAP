clc
clear

% A = dlmread('optihistory2.txt');
% [~,idx] = sort(A(:,1),'ascend');
% A = A(idx,:);

% B = A(:,6:end);

% out = fcnOBJFUNC2(B(1,:))

%%
C = dlmread('optihistory1.txt');
[~,idx] = sort(C(:,1),'ascend');
C = C(idx,:);

D = C(:,6:end);

out = fcnOBJFUNC(D(1,:))

%% 
% out = fcnOBJFUNCSC()



