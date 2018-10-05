function [HeatLoop,OTM,F1,F2,F3,F4,E2,E3,E4,A6,A7] = HeatLoopQrecover(options,FC,OTM,E1,A2,A5,A6)
F1.T = options.T_motor;
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
E3.H2 = E3.H2 + F1.H2;
H_E3 = H_E2 + property(F1,'h','kJ');
E3.T = find_T(E3, H_E3);
E3.Y_H2O = E1.H2O./(E3.H2 + E1.H2O);
E3.Y_H2 = E2.H2./(E3.H2 + E1.H2O);
HeatLoop.FCQbalance = FC.Qremove - (OTM.heat_added);  
HeatLoop.QavailFC = FC.Qremove; 
A6c = A6;
A6c.T = A2.T+10*ones(10,10); 
HeatLoop.QavailTurb = property(A6c,'h','kJ') - property(A5,'h','kJ');
HeatLoop.Qturbrecover = zeros(10,10); 
HeatLoop.H2heat = zeros(10,10); 
HA6 = zeros(10,10);
for w = 1:10
    for x = 1:10
        if HeatLoop.QavailFC(w,x) < OTM.heat_added(w,x)| OTM.heat_added(w,x) < (HeatLoop.QavailFC(w,x) + HeatLoop.QavailTurb);
            HeatLoop.Qturbrecover(w,x) = OTM.heat_added(w,x) - HeatLoop.QavailFC(w,x);
            A6b.O2 = A6.O2(w,x);
            A6b.N2 = A6.N2(w,x);
            A6b.T = A6.T(w,x);
            A6b.P = 980; 
            
            HA6(w,x) = property(A6b,'h','kJ') - HeatLoop.Qturbrecover(w,x);
            A6.T(w,x) = find_T(A6b,HA6(w,x));
        else if HeatLoop.QavailTurb(w,x) < (OTM.heat_added(w,x) - HeatLoop.QavailFC(w,x))
                HeatLoop.H2heat(w,x) = (OTM.heat_added(w,x) - HeatLoop.QavailFC(w,x) - HeatLoop.QavailTurb(w,x))/(FC.hrxnmol(1,1)); 
                A6.T(w,x) = A2.T(w,x) + 10; 
        end
    end
    end
end
OTM.Qrecover = HeatLoop.Qturbrecover;
[A7,T1_work] = expander(A6,options.Pamb,options.T1_eff);          
OTM.T1_workb = T1_work;
OTM.net_work = OTM.T1_workb  + OTM.C1_work + OTM.C2_work; 
E4.T = F2.T;
E4.P = F4.P;
E4.H2O = E3.H2O - F2.H2O;
HeatLoop.Qremove_fuel = H_E3 - property(F2,'h','kJ') - property(E4,'h','kJ');
HeatLoop.Qexcess = FC.Qremove - OTM.heat_added + OTM.Q_out + HeatLoop.Qremove_fuel - Q_addtl_fuel_heat;
end
