function [FL,B1,F2,F3,F4,E2,E3,E4,HX] = fuel_loop(options,E1,F5,A1)
F2.T =  options.T_motor - 10;
F2.P = F5.P  - options.Blower_dP;
F2.H2 = F5.H2.*options.spu;
F1 = F2;
F1.T = 20; %LH2 storage temp;
FL.motor_cooling = enthalpy(F2) - enthalpy(F1);
F3 = F5;
F3.T = 25.765*log(F3.P/1000)+354.84+5;% saturation temperature of 5% H2O in H2 @ 1MPa
F3.P = F2.P;

[F4,B1] = compressor(F3,F5.P,options.Blower_eff);
%B1.work = FL.blower_work; 
FL.Q_preheat = enthalpy(F5) - enthalpy(F4);%property(F5,'h','kJ') - property(F4,'h','kJ');

E4.T = F3.T;
E4.P = F5.P;
E4.H2O = E1.H2O - F3.H2O;

FL.Qremove_fuel = enthalpy(E1) + enthalpy(F2) - enthalpy(E4) - enthalpy(F3) - FL.Q_preheat; %property(E1,'h','kJ') +  property(F2,'h','kJ') -  property(E4,'h','kJ') -  property(F3,'h','kJ') - FL.Q_preheat;

E2 = E1;
E2.T = F4.T;
H_E2 = enthalpy(E1) - FL.Q_preheat; %property(E1,'h','kJ') - FL.Q_preheat;
E2.T = find_T(E2,H_E2);

E3 = E2;
E3.H2 = E2.H2 + F2.H2;
H_E3 = enthalpy(F2) + H_E2;
E3.T = find_T(E3,H_E3);%find_T(E3,property(F2,'h','kJ') + H_E2);

%% find AC mass flow based on condenser heat transfer
AC.O2 = A1.O2; %initial guess for molar airflow through condenser
AC.N2 = A1.N2;
AC.T = A1.T;
AC.P = A1.P; 
ACout = AC; 
ACout.T = E3.T + 25;%find_T(AC, HACout);
HX.fuel = heat_exchanger(E1,E2,F4,F5,options);
[HX.condenser,AC] = condenser(E3,F3,AC,options);
end%Ends function fuel_loop