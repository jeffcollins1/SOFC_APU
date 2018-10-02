function [OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1)
[A2,OTM.C1_work] = compressor(A1,options.PR_comp.*A1.P,options.C1_eff);
A3 = A2;
A3.T = max(options.T_otm,A2.T);
OTM.heat_added =  property(A3,'h','kJ') - property(A2,'h','kJ');
X_O2 = A3.O2./net_flow(A3);
OTM.Rt = 1-(1-X_O2).*(options.P_perm)./(X_O2.*(options.P_non_perm-options.P_perm)); %Theoretical Percentage Recovery of O2 through OTM

O1.T = A3.T;
O1.P = options.P_perm;
O1.O2 = 0;
A4 = A3;
n = 25;
for i = 1:1:n
    X_O2 = A4.O2./net_flow(A4);
    flux = options.j0_otm./log(options.P0_otm).*log(X_O2.*A4.P./O1.P).*options.OTM_area/n*1e4/60/1000/22.4/1000;%in kmol/s:  NmL/cm^2*min * 1e4 cm^2/m / 60 s/min / 1000 mL/L / 22.4 L/mol / 1000mol/kmol
    O1.O2 = O1.O2 + flux;
    A4.O2 = A4.O2 - flux;
end
OTM.Ra = O1.O2./A1.O2;
OTM.mean_flux = O1.O2./options.OTM_area*1000*22.4*60*1000/1e4;%kmol/s*m^2 * 1000 mol/kmol  * 22.4 L/mol * 60s/min * 1000mL/L / 1e4cm^2/m^2
O3 = O1;
O3.T = options.T_oxygen_pump;
[O4,OTM.C2_work] = compressor(O3,options.P_fc,options.C2_eff);
O5 = O4;
O5.T = options.T_fc;

OTM.Q_oxygen_HX = property(O5,'h','kJ') - property(O4,'h','kJ');
O2 = O1;
H_O2 = property(O1,'h','kJ') - OTM.Q_oxygen_HX;
O2.T = find_T(O2, H_O2);
OTM.A2T = A2.T;
[A5,OTM.T1_work] = expander(A4,A1.P,options.T1_eff);

OTM.work_in = OTM.C1_work + OTM.C2_work;
OTM.net_work = OTM.C1_work + OTM.C2_work + OTM.T1_work;
OTM.Q_out = property(O2,'h','kJ') - property(O3,'h','kJ');
end