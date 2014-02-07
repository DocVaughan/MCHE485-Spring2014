#!/usr/local/bin/octave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mass_spring_step.m
%
% Script to analyze a spring-mass-damper system response to step input in position
% 
% Created: 1/30/13 
%   - Joshua Vaughan 
%   - joshua.vaughan@louisiana.edu
%   - http://www.ucs.louisiana.edu/~jev9637
%
% Modified:
%   *
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The baseline system, a simple mass-spring-damper system
%
%    +---> Y       +---> X
%    |             |
%    |    k     +-----+
%    +---/\/\---|     |
%    |          |  M  |
%    +----]-----|     |
%         c     +-----+
%
% Step input in Y, response of the mass is X

m = 1;          % kg
k = (2*pi)^2;   % N/m (Selected to give an undamped wn of 1Hz)
wn = sqrt(k/m)  % Natural Frequency (rad/s)


z = 0.25;        % Define a desired damping ratio
c = 2*z*wn*m;   % calculate the damping coeff. to create it (N/(m/s))

% Define the system to use in simulation - in transfer function form here
num = [2*z*wn wn^2];
den = [1 2*z*wn wn^2];

sys = tf(num,den);

% Define the time vector and step input
t = 0:0.001:5;                   % 5 seconds
U = [zeros(1,501) ones(1,4500)];  % start the step input at t=0.5s for prettier plotting

% run the simulation
[y,t] = lsim(sys,U,t);         

figure
plot(t,U,'--',t,y)
axis([0 5 -0.1 1.6])
% title('Damped Mass-Spring Response')
xlabel('Time (s)','fontsize',20)
ylabel('Displacement (m)','fontsize',20)
legend('Step Input','Response')
% pause()

% Figure setup for printing, probably will look terrible on screen
% h = figure
% set (h,'papertype', '<custom>')
% set (h,'paperunits','inches');
% set (h,'papersize',[9 6])
% set (h,'paperposition', [0,0,[9 6]])
% set (gca,'position', [0.16, 0.19, 0.8, 0.75])
% set (gca, "fontsize", 24)
% 
% FN = findall(h,'-property','FontName');
% set(FN,'FontName','CMU Serif');
% 
% plot(t,U,'--','LineWidth',10,t,y,'LineWidth',10)
% axis([0 5 -0.1 1.6])
% title('Damped Mass-Spring Response')
% xlabel('Time (s)','fontsize',28)
% ylabel('Displacement (m)','fontsize',28)
% grid on
% set (gca, 'LineWidth',4)
% 
% leg = legend('Step Input','Response');
% set (leg, "FontSize", 20,'FontName','CMU Serif') 



