function [HX, AC] = condenser(hotin,hotout,coldin,options)
%% Use LMTD to estimate required surface area with an assumed overall heat transfer coefficient U
[~,H_in] = enthalpy(hotin);
[H_out,~] = enthalpy(hotout);
cpcoldin = SpecHeat(coldin);
coldout = coldin;
dT2 = hotout.T - coldin.T; 
dT1 = dT2 + 10; 
coldout.T = hotin.T -10;
cpcoldout = SpecHeat(coldout);
cp = 0.5.*(cpcoldin + cpcoldout);
Q = H_in - H_out + hotin.H2O*2257/18.015; %property(hotin,'h','kJ') - property(hotout,'h','kJ');
molflow = Q./(cp.*dT1); 
fraction = molflow./net_flow(coldin);
coldin.O2 = fraction.*coldin.O2;
coldin.N2 = fraction.*coldin.N2;
AC = coldin;
%dT1 = max(10,hotin.T - coldout.T); %LMTD Values for counterflow HX
%dT2 = hotout.T - coldin.T; 
LMTD = (dT1 - dT2)./log(dT1./dT2); 
A = Q./(options.hx_U.*LMTD); %Surface area based on range of assumed overall heat transfer coefficients, U
HX.mass = options.hx_t.*options.hx_mat_density.*A; 

%% Calculate heat exchanger effectiveness
Ch = net_flow(hotin).*SpecHeat(hotin);%property(hotin,'C','kJ/(kmol K)'); % Net flow * Constant pressure specific heat for hot inlet
Cc = net_flow(coldin).*SpecHeat(coldin);%property(coldin,'C','kJ/(kmol K)'); %Constant pressure specific heat for cold inlet
Qmax = min(Ch,Cc).*(hotin.T - coldin.T);  
HX.effectiveness = Q./Qmax; % effectiveness, should not be greater than 1