function [air_in_cruise,air_in_climb,H2_used_cruise,H2_used_climb,P_airin_cruise,Pairin_climb,P_net_cruise,P_net_climb] = sizer_od(OD,param_od)
% P_den = param_od.P_den;
% maximum = max(max(P_den));
% [x,y] = find(P_den==maximum); 
Perr_cruise = param_od.Ppropmax_cruise - param_od.NetPower; 
Perr_climb = param_od.Pprop_climb - param_od.NetPower;
minimum1 = min(min(abs(Perr_cruise)))
minimum2 = min(min(abs(Perr_climb))); 
[w,x] = find(Perr_cruise==minimum1);
[y,z] = find(Perr_climb==minimum2); 
P_net_cruise = param_od.NetPower(w,x); 
P_airin_cruise = OD.P_non_perm(w,x);
air_in_cruise = param_od.OTMflux(w,x); 
H2_used_cruise = param_od.H2_used(w,x)*2; %Fuel consumption, kg/s 
P_net_climb = param_od.NetPower(y,z); 
P_airin_climb = OD.P_non_perm(y,z);
air_in_climb = param_od.OTMflux(y,z); 
H2_used_climb = param_od.H2_used(y,z)*2; %Fuel consumption, kg/s 

end