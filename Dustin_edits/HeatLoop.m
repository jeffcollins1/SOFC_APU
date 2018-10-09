function [HeatLoop,F2,F3,F4,F5,E2,E3,E4,HX] = HeatLoop(options,FC,OTM,E1,A1,O2,O3)
F2.T =  options.T_motor;
F2.P = options.P_fc  - options.Blower_dP;
F2.H2 = FC.H2_used;
Tdb_K = linspace(275,423);
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
F3.T = interp1(satP,Tdb_K,options.steamratio.*F2.P)+.02;%temperature of condensation at supply water ratio
F3.P = F2.P;
F3.H2 = FC.H2_supply;
F3.H2O = FC.H2O_supply;

AC.O2 = 0.5*0.21*ones(10,10); 
AC.N2 = 0.5*0.79*ones(10,10);
AC.T = A1.T;
AC.P = A1.P;
[F4,HeatLoop.blower_work] = compressor(F3,options.P_fc,options.Blower_eff);

F5 = F4;
F5.T = options.T_fc;
HeatLoop.Q_preheat = property(F5,'h','kJ') - property(F4,'h','kJ');
E2 = E1;

E2.T = F4.T;
Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
H_E2 = property(E1,'h','kJ') - HeatLoop.Q_preheat;
E2.T = E1.T - HeatLoop.Q_preheat./Q_removed.*(E1.T - F4.T);
E2.T = find_T(E2, H_E2);

HeatLoop.Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');

E3 = E2;
E3.H2 = E2.H2 + F2.H2;
H_E3 = H_E2 + property(F2,'h','kJ');
P_h2O_E3_sat = E3.H2O./net_flow(E3).*E3.P;
E3.T = interp1(satP,Tdb_K,P_h2O_E3_sat);%temperature of condensation with exhaust water concentration
H_E3_sat = property(E3,'h','kJ');
H2O_condense = max(0,(H_E3 - H_E3_sat)./(2260*18));%latent heat of water = 2260kJ/kg
E3.H2O = E3.H2O - H2O_condense;
water.H2O = H2O_condense;
water.T = E2.T;
water.P = E2.P;
E3_Tsat = E3.T;
E3.T = find_T(E3, H_E3-  property(water,'h','kJ'));
E3.T(H2O_condense>0) = E3_Tsat(H2O_condense>0);
HeatLoop.FCQbalance = FC.Qremove - OTM.heat_added; 
E4.T = F3.T;
E4.P = F5.P;
E4.H2O = E3.H2O - F3.H2O;
HeatLoop.Qremove_fuel = H_E3 - property(F3,'h','kJ') - property(E4,'h','kJ');
HACout = property(AC,'h','kJ') + HeatLoop.Qremove_fuel; 
ACout = AC; 
ACout.T = find_T(AC, HACout);

HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel;% - Q_addtl_fuel_heat;
[HX.fuel] = heatexchanger(E1,E2,F4,F5);
HX.condenser = heatexchanger(E3,F3,AC,ACout); 
HACout2 = property(AC,'h','kJ') + OTM.Q_out;
ACout2 = AC;
ACout2.T = find_T(ACout2,HACout2); 
HX.oxycompressor = heatexchanger(O2,O3,AC,ACout2); 
HX.HP = FC.Qremove./(0.628); %Weight of heat pipes
end
