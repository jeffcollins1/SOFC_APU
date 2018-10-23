function [OTM,C2,A3,A4,O1,O2,O3,O4,O5] = OxygenModule(options,A2)
A3 = A2;
A3.T = max(options.T_otm,A2.T);
OTM.heat_added = enthalpy(A3) - enthalpy(A2);  %property(A3,'h','kJ') - property(A2,'h','kJ');
X_O2 = A3.O2./net_flow(A3);
OTM.Rt = 1-(1-X_O2).*(options.P_perm)./(X_O2.*(A3.P-options.P_perm)); %Theoretical Percentage Recovery of O2 through OTM

O1.T = A3.T;
O1.P = options.P_perm;
O1.O2 = 0;
A4 = A3;
A4.T = A3.T + 0.5*options.dT_fc; 
n = 25;
for i = 1:1:n
    X_O2 = A4.O2./net_flow(A4);
    flux = options.j0_otm./log(options.P0_otm).*log(X_O2.*A4.P./O1.P).*options.OTM_area/n*1e4/60/1000/22.4/1000;%in kmol/s:  NmL/cm^2*min * 1e4 cm^2/m / 60 s/min / 1000 mL/L / 22.4 L/mol / 1000mol/kmol
    O1.O2 = O1.O2 + flux;
    A4.O2 = A4.O2 - flux;
end
OTM.Ra = O1.O2./A2.O2;
OTM.toal_flux = O1.O2;
OTM.mean_flux = OTM.toal_flux./options.OTM_area*1000*22.4*60*1000/1e4;%kmol/s*m^2 * 1000 mol/kmol  * 22.4 L/mol * 60s/min * 1000mL/L / 1e4cm^2/m^2

O3 = O1;
O3.T = options.T_oxygen_pump;
[O4,C2] = compressor(O3,options.P_fc,options.C2_eff);
O5 = O4;
O5.T = options.T_fc - .5*options.dT_fc;

OTM.Q_oxygen_HX = enthalpy(O5) - enthalpy(O4); %property(O5,'h','kJ') - property(O4,'h','kJ');
O2 = O1;
H_O2 = enthalpy(O1) - OTM.Q_oxygen_HX; %property(O1,'h','kJ') - OTM.Q_oxygen_HX;
O2.T = find_T(O2, H_O2);
OTM.Q_out = enthalpy(O2) - enthalpy(O3);  %property(O2,'h','kJ') - property(O3,'h','kJ');
end