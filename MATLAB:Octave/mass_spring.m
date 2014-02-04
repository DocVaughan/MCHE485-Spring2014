% #!/usr/local/bin/octave 
% uncomment the line above to allow running this directly from the terminal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to analyze a mass-spring system
% 
% Created: 1/23/13 
%	- Joshua Vaughan 
%	- joshua.vaughan@louisiana.edu
%	- http://www.ucs.louisiana.edu/~jev9637
%
% Modified:
%	* 2/3/14 - Joshua Vaughan - joshua.vaughan@louisiana.edu
%       - better commenting
%       - more explicit (less clever) coding style in some places
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with a simple mass-spring system
%				   +---> X
%				   |
%	/|			+-----+
%	/|	   k	|	  |
%	/|---/\/\---|  M  |
%	/|			|	  |
%	/|			+-----+

% Define system parameters
k = (2*pi)^2;		% Spring constant (N/m)
m = 1;				% Mass (kg) 

wn = sqrt(k/m);		% natural frequency (rad/s)

% Define the system to use in simulation - in state space form here
A = [0 1;
	-k/m 0];

B = [0;1];

C = [1 0;0 1];

sys = ss(A,B,C,[]);


% Set up simulation parameters
t = 0:0.01:5;			% time for simulation
U = zeros(1,length(t)); % no input, just initial condition simulation
x0 = [1;0];				% initial condition, x=1, x_dot=0

% run the simulation
[y,t] = lsim(sys,U,t,x0);		 

figure
plot(t,y(:,1))
title('Undamped Response')
xlabel('Time (s)')
ylabel('Displacement (m)')


% Add a second spring to the simple mass-spring system
%				   +---> X
%				   |
%	/|			+-----+			 |\
%	/|	  k1	|	  |	   k2	 |\
%	/|---/\/\---|  M  |---/\/\---|\
%	/|			|	  |			 |\
%	/|			+-----+			 |\

% let the springs be equal
k1 = k;
k2 = k;

% Define the system to use in simulation - in state space form here
A = [0 1;
	-(k1+k2)/m 0];

B = [0;1];

C = [1 0;0 1];

sys2 = ss(A,B,C,[]);

% run the simulation
[y2,t] = lsim(sys2,U,t,x0);		   

figure
plot(t,y(:,1),t,y2(:,1))
axis([0 5 0 1.25])
title('Mass-Spring Response')
xlabel('Time (s)')
ylabel('Displacement (m)')
legend('One Spring','Two Springs')


% Add a damper to the  simple mass-spring system
%				   +---> X
%				   |
%	/|	  k		+-----+
%	/|---/\/\---|	  |
%	/|			|  M  |
%	/|----]-----|	  |
%	/|	  c		+-----+

% Define system parameters
k = (2*pi)^2;		% Spring constant (N/m)
m = 1;				% Mass (kg) 

wn = sqrt(k/m);		% natural frequency (rad/s)

z = 0.2;            % Define the damping ratio
c = 2*z*wn*m;       % calculate the damping coeff from the damping ratio and nat. freq

% Define the system to use in simulation - in state space form here
A = [0 1;
	-k/m -c/m];

B = [0;1];

C = [1 0;0 1];

sys_damp = ss(A,B,C,[]);

% run the simulation
[y_damp,t] = lsim(sys_damp,U,t,x0);		   

figure
plot(t,y(:,1),t,y_damp(:,1))
axis([0 5 -1.25 1.25])
title('Damped Mass-Spring Response')
xlabel('Time (s)')
ylabel('Displacement (m)')
legend('Undamped','Damped')



