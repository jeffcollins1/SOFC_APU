
function [P,T] =intake(options.height,options.velocity)
Pamb = 107*exp(-0.0001*height)- 10; %Ambient Pressure as a function of altitude
SpeedSound = -0.004*height + 340.62; % Speed of sound at altitude
D = 1.3159*exp(-0.0001*height); %Air density in kg/m^3
M =options.velocity/SpeedSound; %Mach Number
P = Pamb*(1 + 0.5*(0.35)*M^2)^(1.35/0.35) %Isentropic Stagnation of air at a given velocity
end
