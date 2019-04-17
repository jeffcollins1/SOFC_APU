function HX = heat_exchanger(hotin,hotout,coldin,coldout,options)
%% Use LMTD to estimate required surface area with an assumed overall heat transfer coefficient U
Q = enthalpy(hotin) - enthalpy(hotout); %property(hotin,'h','kJ') - property(hotout,'h','kJ');
dT1 = max(10,hotin.T - coldout.T); %LMTD Values for counterflow HX
dT2 = hotout.T - coldin.T; 
LMTD = (dT1 - dT2)./log(dT1./dT2); 
A = Q./(options.hx_U.*LMTD); %Surface area based on range of assumed overall heat transfer coefficients, U
HX.mass = options.hx_t.*options.hx_mat_density.*A; 

%% Calculate heat exchanger effectiveness
Ch = net_flow(hotin).*SpecHeat(hotin);%property(hotin,'C','kJ/(kmol K)'); % Net flow * Constant pressure specific heat for hot inlet
Cc = net_flow(coldin).*SpecHeat(coldin);%property(coldin,'C','kJ/(kmol K)'); %Constant pressure specific heat for cold inlet
Qmax = min(Ch,Cc).*(hotin.T - coldin.T);  
HX.effectiveness = Q./Qmax; % effectiveness, should not be greater than 1
end