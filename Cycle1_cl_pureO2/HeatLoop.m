function [HeatLoop,F1,F2,F3,F4,E2,E3,E4,HX] = HeatLoop(options,FC,OTM,E1)
F1.T =  216*ones(10,10); %options.T_motor + 100;
F1.P = options.P_fc  - options.Blower_dP;
F1.H2 = FC.H2_used;

F2.T = 354*ones(10,10);%temperature of condensation with water partial pressure of 50kPa
F2.P = options.P_fc  - options.Blower_dP;
F2.H2 = FC.H2_supply;
F2.H2O = FC.H2O_supply;

[F3,HeatLoop.blower_work] = compressor(F2,options.P_fc,options.Blower_eff);

F4 = F3;
F4.T = options.T_fc;
Q_preheat = property(F4,'h','kJ') - property(F3,'h','kJ');
E2 = E1;
E2.T = 400;

E2.T = F3.T;
Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
HeatLoop.Q_removed = Q_removed; 
% if Q_removed > Q_preheat
%     H_E2 = property(E1,'h','kJ') -Q_preheat;
%     E2.T = E1.T - Q_preheat/Q_removed*(E1.T-F3.T);
%     E2.T = find_T(E2,H_E2);
% else Q_addtl_fuel_heat = Q_preheat - Q_removed;
% end
need_heat = Q_removed < Q_preheat;
H_E2 = property(E1,'h','kJ') - Q_preheat;
E2.T = E1.T - Q_preheat./Q_removed.*(E1.T - F3.T);
E2.T = find_T(E2, H_E2);
E2.T(need_heat) = F3.T(need_heat);
Q_addtl_fuel_heat = max(0,Q_preheat - Q_removed);

HeatLoop.Q_preheat = Q_preheat; 
E3 = E2;
E3.H2 = E2.H2 + F1.H2;
H_E3 = H_E2 + property(F1,'h','kJ');
E3.T = find_T(E3, H_E3);
E3.Y_H2O = E1.H2O./(E3.H2 + E1.H2O);
E3.Y_H2 = E2.H2./(E3.H2 + E1.H2O);
HeatLoop.FCQbalance = FC.Qremove - OTM.heat_added; 
E4.T = F2.T;
E4.P = F4.P;
E4.H2O = E3.H2O - F2.H2O;
HeatLoop.Qremove_fuel = H_E3 - property(F2,'h','kJ') - property(E4,'h','kJ');
HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel - Q_addtl_fuel_heat;
[HX.fuel] = heatexchanger(E1,E2,F3,F4)
end
