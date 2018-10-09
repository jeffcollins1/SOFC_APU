function [HX] = heatexchanger(hotin,hotout,coldin,coldout)
Thin = hotin.T;
Thout = hotout.T;
Tcin = coldin.T;
Tcout = coldout.T;
%Use LMTD to estimate required surface area with an assumed overall heat transfer
%coefficient U

mflowh = net_flow(hotin); %mass flow of into hot side
mflowc = net_flow(coldin); %mass flow into cold side
Cph = property(hotin,'C','kJ/(kg K)'); %Constant pressure specific heat for hot inlet
Cpc = property(coldin,'C','kJ/(kg K)'); %Constant pressure specific heat for cold inlet
Ch = mflowh.*Cph;
Cc = mflowc.*Cpc; 
Cmin = min(min(Ch,Cc));
Qmax = Cmin*(Thin - Tcin); 
Q = property(hotin,'h','kJ') - property(hotout,'h','kJ'); 
e = Q./Qmax; % effectiveness, should not be greater than 1
if e >1
    HX = nan
else
dT1 = Thin - Tcout; %LMTD Values for counterflow HX
dT2 = Thout - Tcin; 
LMTD = (dT1 - dT2)./log(dT1./dT2); 
U = 40; %Upper heat transfer performance of a gas-to-gas counterflow HX based on Heat and Mass transfer, Cengel, 4e
As = Q./(U.*LMTD); %Surface area based on range of assumed overall heat transfer coefficients, U
t = 0.018; %Total thickness of housing in m, based on NASA estimates, 2005
density = 2300; %Density of silicon carbide, kg/m^3, chosen to replace SS 304 in NASA estimates
HX = As.*(t*density); 
end
end