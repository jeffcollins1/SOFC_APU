function [HeatLoop,F1,F2,F3,F4,F5,E2,E3,E4,E5] = HeatLoop(options,FC,OTM,E1)
F1.T = options.T_motor;
F1.P = options.P_fc  - options.Blower_dP;
F1.H2 = FC.H2_used;

F2.T = 354;%temperature of condensation with water partial pressure of 50kPa
F2.P = options.P_fc  - options.Blower_dP;
F2.H2 = FC.H2_supply;
F2.H2O = FC.H2O_supply;

[F3,HeatLoop.blower_work] = compressor(F2,options.P_fc,options.Blower_eff);
F4 = F3;
F5 = F3;
F5.T = options.T_fc;
Q_preheat = property(F4,'h','kJ') - property(F3,'h','kJ');
cp1 = property(E1,'c','kJ/(kmol K)');
E2 = E1;
E2.T = E1.T - OTM.heat_added/(E1*cp1);
E3 = E2;
E3.T = F3.T + 10;
Q_preheat2 = property(E2,'h','kJ') - property(E3,'h','kJ');
cp2 = property(F3,'c','kJ/(kmol K)');
F4.T = F3.T + Q_preheat2/(F3*cp2); 
E4.T = F2.T;
Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
Q_addtl_fuel_heat = 0;
if Q_removed > Q_preheat
    H_E2 = property(E1,'h','kJ') - Q_preheat;
    E3.T = E1.T - Q_preheat/Q_removed*(E1.T - F3.T);
    E3.T = find_T(E3, H_E3);
else
    Q_addtl_fuel_heat = Q_preheat - Q_removed;
end
HeatLoop.Q_preheat = Q_preheat; 

E4.H2 = E3.H2 + F1.H2;
H_E4 = H_E3 + property(F1,'h','kJ');
E4.T = find_T(E4, H_E4);
E4.Y_H2O = E1.H2O/(E4.H2 + E1.H2O);
E4.Y_H2 = E4.H2/(E4.H2 + E1.H2O);
E5 = E3.H2O - F2.H2O;
E5.T = F2.T;
E5.P = F4.P;

HeatLoop.Qremove_fuel = H_E3 - property(F2,'h','kJ') - property(E5,'h','kJ');
HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel - Q_addtl_fuel_heat;
HeatLoop.F4 = exergy(F4,options.T0,options.P0);
HeatLoop.E2 = exergy(E2,options.T0,options.P0);
HeatLoop.F3 = exergy(F3,options.T0,options.P0); 
HeatLoop.E1 = exergy(E1,options.T0,options.P0);
HeatLoop.E3 = exergy(E3,options.T0,options.P0);
HeatLoop.E4 = exergy(E4,options.T0,options.P0);
HeatLoop.E5 = exergy(E5,options.T0,options.P0);
end