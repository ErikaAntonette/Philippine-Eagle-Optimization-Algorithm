clc; clear; format longG; 

D = 2; % dimension number
f = @(x) x(1)^2 + x(2)^2; % function to be optimized
Space_x_max = 500*ones(1,D); % maximum bounds, must be a row vector
Space_x_min = -500*ones(1,D); % minimum bounds, must be a row vector

%% EGGHOLDER FUNCTION
% D = 2;
% f = @(x) -(x(2)+47) * sin(sqrt(abs(x(2)+x(1)/2+47))) + -x(1) * sin(sqrt(abs(x(1)-(x(2)+47))));
% Space_x_max = 512*ones(1,D); % maximum bounds, must be a row vector
% Space_x_min = -512*ones(1,D); % minimum bounds, must be a row vector

%% EASOM FUNCTION
% D = 2;
% f = @(x) -cos(x(1))*cos(x(2))*exp(-(x(1)-pi).^2 - (x(2) - pi).^2);
% Space_x_max = 100*ones(1,D); % maximum bounds, must be a row vector
% Space_x_min = -100*ones(1,D); % minimum bounds, must be a row vector

%% BEALE FUNCTION
% D = 2;
% f = @(x) (1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1) + x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;
% Space_x_max = 4.5*ones(1,D); % maximum bounds, must be a row vector
% Space_x_min = -4.5*ones(1,D); % minimum bounds, must be a row vector

%%
[fbest_pheagle, xbest_pheagle] = pheaglealgorithm(D, f, Space_x_max, Space_x_min)
