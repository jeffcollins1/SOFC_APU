function [FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1)
F2.T =  options.T_motor;
F2.P = F5.P  - options.Blower_dP;
F2.H2 = F5.H2.*options.spu;


F3 = F5;
F3.T = 25.765*log(F3.P/1000)+354.84+5;% saturation temperature of 5% H2O in H2 @ 1MPa
F3.P = F2.P;

[F4,FL.blower_work] = compressor(F3,F5.P,options.Blower_eff);
B1 = FL.blower_work; 
FL.Q_preheat = property(F5,'h','kJ') - property(F4,'h','kJ');

E4.T = F3.T;
E4.P = F5.P;
E4.H2O = E1.H2O - F3.H2O;

FL.Qremove_fuel =  property(E1,'h','kJ') +  property(F2,'h','kJ') -  property(E4,'h','kJ') -  property(F3,'h','kJ') - FL.Q_preheat;

E2 = E1;
E2.T = F4.T;
H_E2 = property(E1,'h','kJ') - FL.Q_preheat;
E2.T = find_T(E2,H_E2);

E3 = E2;
E3.H2 = E2.H2 + F2.H2;
E3.T = find_T(E3,property(F2,'h','kJ') + H_E2);

% Tdb_K = linspace(275,500);
% satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
% E2.T = F4.T;
% Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
% H_E2 = property(E1,'h','kJ') - FL.Q_preheat;
% E2_sat = E2;
% P_h2O_E2_sat = E2.H2O./net_flow(E2).*E2.P;
% E2_sat.T = interp1(satP,Tdb_K,P_h2O_E2_sat)+.1;%temperature of condensation with exhaust water concentration
% E2_cp = property(E2_sat,'cp','kJ/K');
% H_E2_sat = property(E2_sat,'h','kJ');
% H2O_condense = min(E2.H2O,max(0,(H_E2_sat - H_E2)./(2260*18)));%latent heat of water = 2260kJ/kg
% E2.H2O = E2.H2O - H2O_condense;
% water.H2O = H2O_condense;
% water.T = E2.T;
% water.P = E2.P;
% E2.T = E2_sat.T+.01;
% % E2.T = E2_sat.T + (H_E2 - property(E2_sat,'h','kJ'))./E2_cp;
% % E2.T = E1.T - FL.Q_preheat./Q_removed.*(E1.T - F4.T);
% E2.T = find_T(E2, H_E2);
% 
% FL.Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
% 
% E3 = E2;
% E3.H2 = E2.H2 + F2.H2;
% H_E3 = H_E2 + property(F2,'h','kJ');
% E3_sat = E3;
% P_h2O_E3_sat = E3.H2O./net_flow(E3).*E3.P;
% E3_sat.T = interp1(satP,Tdb_K,P_h2O_E3_sat)+.1;%temperature of condensation with exhaust water concentration
% H_E3_sat = property(E3_sat,'h','kJ');
% H2O_condense = min(E3.H2O,max(0,(H_E3_sat - H_E3)./(2260*18)));%latent heat of water = 2260kJ/kg
% E3.H2O = E3.H2O - H2O_condense;
% water.H2O = H2O_condense;
% water.T = E3.T;
% water.P = E3.P;
% E3.T = E3_sat.T+.1;
% E3.T = find_T(E3, H_E3 -  property(water,'h','kJ'));
% E3.T(H2O_condense>0) = E3_sat.T(H2O_condense>0);
% E4.T = F3.T;
% E4.P = F5.P;
% E4.H2O = E3.H2O - F3.H2O;
% FL.Qremove_fuel = H_E3 - property(F3,'h','kJ') - property(E4,'h','kJ');

%% find AC mass flow based on condensor heat transfer
AC.O2 = 3*A1.O2;
AC.N2 = 3*A1.N2;
AC.T = A1.T;
AC.P = A1.P; 
ACout = AC; 
ACout.T = E3.T + 25;%find_T(AC, HACout);
HX.fuel = heat_exchanger(E1,E2,F4,F5,options);
HX.condenser = heat_exchanger(E3,F3,AC,ACout,options);
end%Ends function fuel_loop