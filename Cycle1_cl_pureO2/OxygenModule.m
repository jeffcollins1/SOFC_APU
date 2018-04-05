function [OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1)

[A2,C1_work] = compressor(A1,options.P_non_perm,options.C1_eff);
A3 = A2;
A3.T = options.T_otm;
OTM.heat_added =  property(A3,'h','kJ') - property(A2,'h','kJ');
X_O2 = A3.O2/net_flow(A3);
OTM.Rt = 1-(1-X_O2)*(options.P_perm)/(X_O2*(options.P_non_perm-options.P_perm)); %Theoretical Percentage Recovery of O2 through OTM
OTM.Ra = options.OTM_perc_theoretical*OTM.Rt;
O1.T = A3.T;
O1.P = options.P_perm;
O1.O2 = OTM.Ra*A3.O2;
OTM.flux = O1.O2*1000/options.OTM_area;

O3 = O1;
O3.T = options.T_oxygen_pump;
[O4,C2_work] = compressor(O3,options.P_fc,options.C2_eff);
O5 = O4;
O5.T = options.T_fc;

Q_oxygen_HX =  property(O5,'h','kJ') - property(O4,'h','kJ');

O2 = O1;
H_O2 = property(O1,'h','kJ') - Q_oxygen_HX;
O2.T = find_T(O2, H_O2);

A4 = A3;
A4.O2 = A4.O2-O3.O2;
[A5,T1_work] = expander(A4,A1.P,options.T1_eff);

OTM.net_work = C1_work + C2_work + T1_work;