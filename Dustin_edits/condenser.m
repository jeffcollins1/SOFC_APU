function [HX, cold_in] = condenser(hot_in,hot_out,cold_in,options)
%% Use LMTD to estimate required surface area with an assumed overall heat transfer coefficient U
[~,H_in] = enthalpy(hot_in);
[~,H_out] = enthalpy(hot_out);
Q = H_in - H_out + hot_in.H2O*2257/18.015; %sensible enthalpy difference + latent heat of vaporization
cpcoldin = SpecHeat(cold_in);
coldout = cold_in;
dT2 = hot_out.T - cold_in.T; 
dT1 = dT2 + 10; 
coldout.T = hot_in.T -10;
cpcoldout = SpecHeat(coldout);
cp = 0.5.*(cpcoldin + cpcoldout);

molflow = Q./(cp.*dT1); 
fraction = molflow./net_flow(cold_in);
cold_in.O2 = fraction.*cold_in.O2;
cold_in.N2 = fraction.*cold_in.N2;

%dT1 = max(10,hotin.T - coldout.T); %LMTD Values for counterflow HX
%dT2 = hotout.T - coldin.T; 
LMTD = (dT1 - dT2)./log(dT1./dT2); 
A = Q./(options.hx_U.*LMTD); %Surface area based on range of assumed overall heat transfer coefficients, U
HX.mass = options.hx_t.*options.hx_mat_density.*A; 

%% Calculate heat exchanger effectiveness
Ch = net_flow(hot_in).*SpecHeat(hot_in);%property(hotin,'C','kJ/(kmol K)'); % Net flow * Constant pressure specific heat for hot inlet
Cc = net_flow(cold_in).*SpecHeat(cold_in);%property(coldin,'C','kJ/(kmol K)'); %Constant pressure specific heat for cold inlet
Qmax = min(Ch,Cc).*(hot_in.T - cold_in.T);  
HX.effectiveness = Q./Qmax; % effectiveness, should not be greater than 1