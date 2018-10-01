Ai = 0.25*pi*2.2^2; %Inlet area
Ap = 0.25*pi*2^2; %Prop area
Ae = 0.25*pi*1.8^2; %Exit area
gamma = 1.4; %Ratio of specific heats for air
%Ambient temperature, pressure, speed of sound and density through climb and cruise
p0_climb = 107*exp(-0.0001*OD.height)- 10;
p0_cruise = (107.*exp(-0.0001.*OD.height_cruise)- 10).*ones(10,10);
T0_climb =  -0.0065.*OD.height + 14.987 + 273.1;
T0_cruise =  -0.0065.*OD.height_cruise + 14.987 + 273.1;
d1 =  refproparray('D','T',T0_climb(1,1:10),'P',P0_climb(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]); 
d2 =  refproparray('D','T',T0_cruise(1,1:10),'P',P0_cruise(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
d1_climb = (ones(10,1)*d1(1,1:10));
d1_cruise = (ones(10,1)*d2(1,1:10));
a1_climb = refproparray('A','T',T0_climb(1,1:10),'P',P0_climb(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
a1_cruise = refproparray('A','T',T0_cruise(1,1:10),'P',P0_cruise(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
ss_climb = (ones(10,1)*a1_climb(1,1:10)); 
ss_cruise = (ones(10,1)*a1_cruise(1,1:10)); 
%Flight velocity 
V1_climb = OD.velocity_climb;%Magnitude of velocity during climb
V1_cruise = 0.85*ss_cruise; %Velocity during cruise

Mach_climb = V1_climb./ss_climb; %Mach number during climb phase

%Fluid properties at propeller (State 2)
mdot_climb = Ai.*V1_climb.*d1_climb
d2_climb = d1_climb.*(1 + 0.5.*((gamma-1)./gamma)).*Mach_climb.^2).^(1/(gamma-1)); %density of air at propeller
V2_climb = mdot_climb./(Ap.*d2_climb); %Velocity at propeller
p2_climb = p0_climb.*(1 + 0.5.*((gamma-1)./gamma)).*Mach_climb.^2).^(gamma/(gamma-1)); %Stagnation pressure at propeller
a2_climb = refproparray('A','D',d2_climb(1,1:10),'P',p2_climb(1,1:10),'OXYGEN','NITROGEN',[0.21,0.79]);
ss2_climb = (ones(10,1)*a2_climb(1,1:10)); 
roc = 7.5; %Max rate of climb, m/s