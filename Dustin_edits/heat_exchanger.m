function HX = heat_exchanger(hotin,hotout,coldin,coldout)
%% Use LMTD to estimate required surface area with an assumed overall heat transfer coefficient U
Q = property(hotin,'h','kJ') - property(hotout,'h','kJ');
dT1 = max(10,hotin.T - coldout.T); %LMTD Values for counterflow HX
dT2 = hotout.T - coldin.T; 
LMTD = (dT1 - dT2)./log(dT1./dT2); 
U = 40; %Upper heat transfer performance of a gas-to-gas counterflow HX based on Heat and Mass transfer, Cengel, 4e
A = Q./(U.*LMTD); %Surface area based on range of assumed overall heat transfer coefficients, U
t = 0.0018; %Total thickness of housing in m, based on NASA estimates, 2005
density = 2700; %Density of sintered silicon carbide, kg/m^3, chosen to replace SS 304 in NASA estimates with same plate and housing thickness
HX.mass = A*t*density; 

%% Calculate heat exchanger effectiveness
Ch = net_flow(hotin).*property(hotin,'C','kJ/(kmol K)'); % Net flow * Constant pressure specific heat for hot inlet
Cc = net_flow(coldin).*property(coldin,'C','kJ/(kmol K)'); %Constant pressure specific heat for cold inlet
Qmax = min(Ch,Cc).*(hotin.T - coldin.T);  
HX.effectiveness = Q./Qmax; % effectiveness, should not be greater than 1
end