% This script modifies the modifies geometry to the format accepted by the
% VAPTOR input file.

% The input data must should be columns in the following format:
% [r_R, c_R, Beta(deg), Midchordline]
clear,clc

%% Load input data
rotor = importdata('T_motor_w_zerolift.mat');
rotor_raduis = 0.4572/2; % Rotor raduis (m)

%% Calcualte each input required
% Calculate rotor raduis length and chord length (m)
r = rotor_raduis*rotor.r_R;
c = rotor_raduis*rotor.c_R;
len = size(r,1);

% Calculate x points
x = -rotor.MidChordLine-c/2;
x = reshape(repmat(x',[2,1]),[(2*len),1]);
x(1) = [];
x(end) = [];
x = reshape(x,2,len-1);

% Calculate y points
y = reshape(repmat(r',[2,1]),[(2*len),1]);
y(1) = [];
y(end) = [];
y = reshape(y,2,len-1);

% Aline leading edge
z = zeros(2,len-1);

% Reshape chords
c = reshape(repmat(c',[2,1]),[(2*len),1]);
c(1) = [];
c(end) = [];
c = reshape(c',2,len-1);

% Reshape beta
epsilon = reshape(repmat(rotor.Beta',[2,1]),[(2*len),1]);
epsilon(1) = [];
epsilon(end) = [];
epsilon = reshape(epsilon',2,len-1);

%% Output the 
for i = 1:(len-1)
test = fprintf('Panel #:%d\nNumber of spanwise elembers:\tvecN\t\t= 1\nNumber of chordwise elements:\tvecM\t\t= 1\nAirfoil number:\t\t\t\t\tvecAIRFOIL\t= 1;\nSymmetry edge (0,1 or 2):\t\tvecSYM\t\t= 0;\nx\t\t\ty\t\t\tz\t\t\tchord\t\tepsilon\n%f\t%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\t%f\n\n',i,x(1,i),y(1,i),z(1,i),c(1,i),epsilon(1,i),x(2,i),y(2,i),z(2,i),c(2,i),epsilon(2,i));
end

figure(1)
clf(1)
hold on
plot(x(1,:),y(2,:))
axis equal
grid on
hold off