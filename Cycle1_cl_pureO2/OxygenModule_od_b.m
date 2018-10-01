function [OTM_od,A2,A3,A3b,A4,A5,O1,O2,O3,O4,O5] = OxygenModule_od_b(OD,A1,bypass,z)
[A2,C1_work] = compressor(A1,OD.P_non_perm,OD.C1_eff);
A3.O2 = (1-bypass).*A2.O2;
A3.N2 = (1-bypass).*A2.N2;
A3b.O2 = (bypass).*A2.O2;
A3b.N2 = (bypass).*A2.N2; 
A3b.T = OD.T_otm;
A3b.P = A2.P; 
A3.T = OD.T_otm;
A3.P = A2.P;
AX.N2 = A2.N2;
AX.O2 = A2.O2;
AX.T = OD.T_otm;
AX.P = A2.P;
OTM_od.heat_added =  property(AX,'h','kJ') - property(A2,'h','kJ');
X_O2 = A3.O2./net_flow(A3);
OTM_od.Rt = 1-(1-X_O2).*(OD.P_perm)./(X_O2.*(OD.P_non_perm-OD.P_perm)); %Theoretical Percentage Recovery of O2 through OTM

OTM_od.Ra = OD.OTM_perc_theoretical.*OTM_od.Rt;
O1.O2 = OTM_od.Ra.*A3.O2; 
O1.T = A3.T;
O1.P = OD.P_perm*ones(10,10);

%OTM_od.flux = O1.O2*1000./OD.OTM_area;
%OTM_od.cells = OD.OTM_area.*10000./81; %Number of cells as a function of OTM area;
O3 = O1;
O3.T = OD.T_oxygen_pump;
[O4,C2_work] = compressor(O3,OD.P_fc,OD.C2_eff);
O5 = O4;
O5.T = OD.T_fc;

Q_oxygen_HX =  property(O5,'h','kJ') - property(O4,'h','kJ');
OTM_od.Q_oxygen_HX = Q_oxygen_HX; 
O2 = O1;
H_O2 = property(O1,'h','kJ') - Q_oxygen_HX;
O2.T = find_T(O2, H_O2);
OTM_od.A2T = A2.T;
A4.N2 = A3.N2 + A3b.N2;
A4.O2 = A3.O2 + A3b.O2 -O3.O2;
A4.T = A3.T;
A4.P = A3.P; 
A4.Y_O2 = A4.O2./(A4.O2 + A4.N2); %Molar fraction of O2 in non permeate stream
A4.Y_N2 = A4.N2./(A4.O2 + A4.N2);
[A5,T1_work] = expander(A4,OD.P0,OD.T1_eff);
OTM_od.C1_work = C1_work;
OTM_od.T1_work = T1_work;
OTM_od.C2_work = C2_work; 
OTM_od.work_in = C1_work + C2_work;
OTM_od.net_work = C1_work + C2_work + T1_work;
%OTM.O5 = exergy(O5,options.T0,options.P0);
OTM_od.Q_out = property(O2,'h','kJ') - property(O3,'h','kJ');

end