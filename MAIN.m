
clc
clear

lb = [];
ub = [];
A = [];
b = [];







Aeq = [];
beq = [];

nvars = 21;
TolCon_Data = 1e-6; 
TolFun_Data = 1e-08;

options = gaoptimset;
options = gaoptimset(options,'TolFun', TolFun_Data);
% options = gaoptimset(options,'CrossoverFraction', 0.5);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'InitialPopulation', Seed);
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv @gaplotexpectation @gaplotscorediversity @gaplotstopping });
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'UseParallel', 1 );
options = gaoptimset(options,'Generations',1000,'StallGenLimit', 50);
[x,fval,exitflag,output,population,score] = gamultiobj(@objectiveFunction,nvars,A,b,Aeq,beq,lb,ub,[],options);