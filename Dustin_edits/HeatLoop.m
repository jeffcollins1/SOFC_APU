function [HeatLoop,B1,F2,F3,F4,F5,E2,E3,E4] = HeatLoop(options,FC,OTM,E1)
F2.T =  options.T_motor;
F2.P = FC.pressure  - options.Blower_dP;
F2.H2 = FC.H2_used;

Tdb_K = linspace(275,423);
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
F3.T = interp1(satP,Tdb_K,options.steamratio.*F2.P)+.02;%temperature of condensation at supply water ratio
F3.P = F2.P;
F3.H2 = FC.H2_supply;
F3.H2O = FC.H2O_supply;

[F4,B1] = compressor(F3,FC.pressure,options.Blower_eff);

F5 = F4;
F5.T = options.T_fc - .5*options.dT_fc;
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
E3_sat = E3;
P_h2O_E3_sat = E3.H2O./net_flow(E3).*E3.P;
E3_sat.T = interp1(satP,Tdb_K,P_h2O_E3_sat)+.01;%temperature of condensation with exhaust water concentration
H_E3_sat = property(E3_sat,'h','kJ');
% E3.T = 0.5*E3.T + 0.5*E3_sat.T;
H2O_condense = min(E3.H2O,max(0,(H_E3_sat - H_E3)./(2260*18)));%latent heat of water = 2260kJ/kg
E3.H2O = E3.H2O - H2O_condense;
water.H2O = H2O_condense;
water.T = E2.T;
water.P = E2.P;
E3.T = E3_sat.T+.01;
E3.T = find_T(E3, H_E3 -  property(water,'h','kJ'));
E3.T(H2O_condense>0) = E3_sat.T(H2O_condense>0);
E4.T = F3.T;
E4.P = F5.P;
E4.H2O = E3.H2O - F3.H2O;
HeatLoop.Qremove_fuel = H_E3 - property(F3,'h','kJ') - property(E4,'h','kJ');
if ~isempty(OTM)
    HeatLoop.FCQbalance = FC.Qremove - OTM.heat_added; 
    HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel;% - Q_addtl_fuel_heat;
end
end%Ends function HeatLoop