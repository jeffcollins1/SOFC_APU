function [OTM,A1,A2,A3,A4,A5,O1,O2,O3] = OxygenModule(options,A1)

[A2,C1_work] = compressor(A1,options.P_non_perm,options.C1_eff);
A3 = A2;
A3.T = options.T_otm;
OTM.heat_added = enthalpy(A3) - enthalpy(A2);
X_O2 = A3.O2/net_flow(A3);
OTM.Rt = 1-(1-X_O2)*(options.P_perm)/(X_O2*(options.P_non_perm-options.P_perm)); %Theoretical Percentage Recovery of O2 through OTM
OTM.Ra = options.OTM_perc_theoretical*OTM.Rt;
O1.T = A3.T;
O3.P = A3.P;
O3.O2 = OTM.Ra*A3.O2;
A4 = A3;
A4.O2 = A4.O2-O3.O2;
[A5,T1_work] = expander(A4,A1.P,options.T1_eff);
%%need iterative loop to find heat rejection + compression that arrives at
%%correct sofc inlet temperature

[O3,C2_work] = compressor(O2,options.P_sofc,options.C2_eff);

OTM.net_work = C1_work + C2_work + T1_work;