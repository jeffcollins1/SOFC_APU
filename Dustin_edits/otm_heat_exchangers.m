function HX = otm_heat_exchangers(options,FC,OTM,HX,A1,O1,O2,O3,O4,O5)
HX.oxygen = heat_exchanger(O1,O2,O4,O5,options);
AC.O2 = 0.33*A1.O2;
AC.N2 = 0.33*A1.N2;
AC.T = A1.T; 
AC.P = A1.P;
ACout.O2 = AC.O2;
ACout.N2 = AC.N2;
ACout.P = A1.P; 
ACout.T = find_T(AC,property(AC,'h','kJ') + OTM.Q_out); 
HX.oxycompressor = heat_exchanger(O2,O3,AC,ACout,options); 
end