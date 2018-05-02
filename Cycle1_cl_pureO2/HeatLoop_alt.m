function [HeatLoop,F1,F2,F3,F4,F5,E2,E3,E4] = HeatLoop_alt(options,FC,OTM,E1)
F1.T = options.T_motor;
F1.P = options.P_fc  - options.Blower_dP;
F1.H2 = FC.H2_used;

F2.T = 354;%temperature of condensation with water partial pressure of 50kPa
F2.P = options.P_fc  - options.Blower_dP;
F2.H2 = FC.H2_supply;
F2.H2O = FC.H2O_supply;

[F3,HeatLoop.blower_work] = compressor(F2,options.P_fc,options.Blower_eff);

F5 = F3;
F5.T = options.T_fc;
%Q_preheat = property(F5,'h','kJ') - property(F3,'h','kJ');
E2 = E1;
E1.T = options.T_fc + options.dT_fc;
cp1 = property(E1,'c','kJ/(mol K)');
E2.T = E1.T - OTM.heat_added/(E1*cp1);
F4.T = E2.T - 10; %F4 temperature based on 10 degree HX pinch point;
E2.T = F3.T;
Q_preheat = property(F5,'h','kJ') - property(F3,'h','kJ');
Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
Q_addtl_fuel_heat = 0;
if Q_removed > Q_preheat
    H_E2 = property(E1,'h','kJ') - Q_preheat;
    E2.T = E1.T - Q_preheat/Q_removed*(E1.T - F3.T);
    E2.T = find_T(E2, H_E2);
else
    Q_addtl_fuel_heat = Q_preheat - Q_removed;
end
HeatLoop.Q_preheat = Q_preheat; 
E3 = E2;
E3.H2 = E3.H2 + F1.H2;
H_E3 = H_E2 + property(F1,'h','kJ');
E3.T = find_T(E3, H_E3);
E3.Y_H2O = E1.H2O/(E3.H2 + E1.H2O);
E3.Y_H2 = E2.H2/(E3.H2 + E1.H2O);

E4.T = F2.T;
E4.P = F5.P;
%E4 = E3.H2O - F2.H2O;
HeatLoop.Qremove_fuel = H_E3 - property(F2,'h','kJ') - property(E4,'h','kJ');
HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel - Q_addtl_fuel_heat;
HeatLoop.F4 = exergy(F5,options.T0,options.P0);
HeatLoop.E2 = exergy(E2,options.T0,options.P0);
HeatLoop.F3 = exergy(F3,options.T0,options.P0); 
HeatLoop.E1 = exergy(E1,options.T0,options.P0);
end