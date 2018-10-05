function [HeatLoop,F1,F2,F3,F4,E1,E2,E3,E4,A2a,A2b,A3,A4,A5] = HeatLoop_3b(intake,options,FCArray,A2)
F1.H2 = FCArray.H2used;
 F1.T = 20*ones(10,10);
 F1.P = 1000*ones(10,10); 
F2.T = options.T_motor;
F2.P = options.P_fc;
F2.H2 = FCArray.H2used;
F3.T = options.Tamb - 10; 
F3.P = options.P_fc;
%HeatLoop.QinFuel1 = property(F3,'h','kJ') - property(F2,'h','kJ');
F3.H2 = FCArray.H2in;

F4.H2 = F3.H2;
F4.H2O = 0.05*F3.H2;
F4.T = options.T_fc; 
F4.P = options.P_fc;
E1.H2 = FCArray.H2out;
E1.H2O = FCArray.H2Oout;
A3.N2 = A2.N2;
A3.O2 = A2.O2;
A3.T = A2.T;
A3.P = A2.P;
A4.O2 = A2.O2;
A4.N2 = A2.N2;
A4.T = options.T_fc;
A4.P = options.P_fc;
A5.O2 = FCArray.O2out;
A5.N2 = FCArray.N2out;
A5.P = 1000*ones(10,10); 
A5.T = 1073*ones(10,10);
E1.P = 980*ones(10,10);
E1.T = 1073*ones(10,10);
E1.H2O = FCArray.H2used;
E1.H2 = FCArray.H2out; 
HeatLoop.Qinfuel1 = property(F3,'h','kJ') - property(F2,'h','kJ'); 
HeatLoop.Qinfuel2 = property(F4,'h','kJ') - property(F3,'h','kJ');
% E2.H2 = E1.H2;
% E2.H2O = E1.H2O;
% E2.P = E1.P; 
% %need_heat = HeatLoop.Qinfuel2;
% H_E2 = property(E1,'h','kJ') - HeatLoop.Qinfuel2;
% %E2.T = E1.T - Q_preheat./Q_removed.*(E1.T - F3.T);
% E2.T = F3.T;
%E2.T = find_T(E2, H_E2);
%E2.T(need_heat) = F3.T(need_heat);
%Q_addtl_fuel_heat = max(0,Q_preheat - Q_removed);
%HeatLoop.Qinfuel2 = property(F3,'h','kJ') - property(F2,'h','kJ');
HeatLoop.Qairin1 = property(A4,'h','kJ') - property(A3,'h','kJ'); %Total preheating required between compressor outlet and FC inlet
Cathode_h_in = property(A3,'h','kJ');
HeatLoop.Qoutanode = property(E1,'h','kJ') - property(F4,'h','kJ'); %Heat Leaving FC by fuel exhaust
HeatLoop.Qoutcathode = property(A5,'h','kJ') - property(A4,'h','kJ'); %Heat leaving fc by cathode exhaust
HeatLoop.Qexcess = FCArray.Qgen - HeatLoop.Qoutanode - HeatLoop.Qoutcathode - HeatLoop.Qairin1; %Excess heat in fuel cell
HeatLoop.Qavailable = FCArray.Qgen - HeatLoop.Qoutanode - HeatLoop.Qoutcathode; %Heat available to incoming air stream
HeatLoop.Qairin2 = zeros(10,10); %Consider using exhaust for preheating 
A2b.N2 = zeros(10,10);
A2b.O2 = zeros(10,10);
A2a.N2 = zeros(10,10);
A2a.O2 = zeros(10,10);
HeatLoop.QinA2b = zeros(10,10);
A2a.T = A2.T(1,1);
A2a.P = A2.P(1,1)
A2b.T = options.T_fc(1,1);
A2b.P = A2.P(1,1); 
for i = 1:10
    for j = 1:10
        if HeatLoop.Qexcess(i,j) < 0 && HeatLoop.Qavailable(i,j) > 0
            HeatLoop.Qairin2(i,j) = -HeatLoop.Qexcess(i,j); 
            A2bl.N2 = A2.N2(i,j)*(1 - HeatLoop.Qavailable(i,j)/(HeatLoop.Qairin1(i,j)));
            A2bl.O2 = A2.O2(i,j)*(1 - HeatLoop.Qavailable(i,j)/(HeatLoop.Qairin1(i,j)));
            A2al.N2 = A2.N2(i,j) - A2bl.N2; 
            A2al.O2 = A2.O2(i,j) - A2bl.O2;
            A2al.T = A2.T(1,1);
            A2al.P = A2.P(1,1);
            HeatLoop.QinA2b(i,j) = Cathode_h_in(i,j) - property(A2al,'h','kJ'); 
            A2a.N2(i,j) = A2al.N2;
        A2a.O2(i,j) = A2al.O2;
        A2b.N2(i,j)= A2bl.N2;
        A2b.O2(i,j) = A2bl.O2;
        end
        
    end
end

HeatLoop.CoolingLoad = property(F2,'h','kJ') - property(F1,'h','kJ'); 
HeatLoop.Qexchangetotal = HeatLoop.Qairin1 + HeatLoop.Qinfuel1 + HeatLoop.Qinfuel2 +HeatLoop.QinA2b;
P_out = options.Pamb; 
%Qexcess = FCArray.Qgen - Total; 
Ti = 1023; 
error1 = 1000;
% Temp = zeros(10,10)
% %Solve for exhaust temperature based on excess heat from fuel cell
% for y = 1:10
%     for z = 1:10
%        if Qexcess(y,z) >0 
%     while  error1 >10
%     A3b.N2 = A3.N2(y,z);
%     A3b.O2 = A3.O2(y,z);
%     A3b.T = 1023;
%     A3b.P = 1000; 
%     F3b.H2 = F3.H2(y,z);
%     F3b.H2O = F3.H2O(y,z);
%     F3b.T = 1023;
%     F3b.P = 1000; 
%     Ti = Ti + 1;
%     E1b.H2 = E1.H2(y,z);
%     E1b.H2O = E1.H2(y,z);
%     E1b.P = 980;
%     A4b.N2 = A4.N2(y,z);
%     A4b.O2 = A4.O2(y,z);
%     A4b.P = 980; 
%     E1b.T = Ti;
%     A4b.T = Ti;
%     Hg = property(E1b,'h','kJ') - property(F3b,'h','kJ') + property(A4b,'h','kJ') - property(A3b,'h','kJ');
%     error1 = Qexcess(y,z) - Hg; 
%     end
%        else Ti = 1023;
%        end
%     E1.T(y,z) = Ti;
%     A4.T(y,z) = Ti;
% end
% end
% 
% 
% E1.T = 1073*ones(10,10);

ET.N2 =  A5.N2;
ET.O2 = A5.O2;
ET.H2 = E1.H2;
ET.H2O = E1.H2O; 
ET.P = E1.P;

dH = zeros(10,10);
%Solve for composition and total enthalpy after combustion
for i = 1:10
    for j = 1:10
if A5.O2(i,j) > 2*E1.H2(i,j)
    dH(i,j) = FCArray.hrxnmol(1,1)*E1.H2(i,j);
    E2.H2(i,j) = 0;
    E2.H2O(i,j) = ET.H2O(i,j) + ET.H2(i,j);
    E2.O2(i,j) = A5.O2(i,j) - 2*E1.H2(i,j);
    E2.N2(i,j) = A5.N2(i,j); 
else 
    dH(i,j) = 2*FCArray.hrxnmol(1,1)*A5.O2(i,j);
    E2.O2(i,j) = 0;
    E2.H2O(i,j) = ET.H2O(i,j) + 2*ET.O2(i,j);
    E2.H2(i,j) = E1.H2(i,j) - 0.5*A5.O2(i,j); 
    E2.N2(i,j) = A5.N2(i,j); 
end
    end
end
 
H = property(E1,'h','kJ') + property(A5,'h','kJ') + dH - HeatLoop.QinA2b; 

T = 1023; 
E2b.P = ET.P(1,1);
Treal = zeros(10,10); 
 for w = 1:10
     error = 1000;
while error>10
    T = T + 1;
    E2b.H2O = E2.H2O(w,1);
    E2b.H2 = E2.H2(w,1);
    E2b.O2 = E2.O2(w,1);
    E2b.N2 = E2.N2(w,1);
    E2b.T = T;
    HET = property(E2b,'h','kJ');
    error = H(w,1) - HET;
   
end
Treal(w,1:10) = T; 
 end
E2.T = Treal; 
E2.P = E1.P;

H_E3 = property(E2,'h','kJ') - HeatLoop.Qinfuel2 - HeatLoop.QinA2b;
%E3.T = E1.T - Q_preheat./Q_removed.*(E1.T - F3.T);
E3 = E2;
E3.P = E2.P;
E3.T = find_T(E3, H_E3);
P_out = options.Pamb; 
[E4,T1_work] = expander(E3,P_out,options.T1_eff);
HeatLoop.T1_work = T1_work; 



%E3 = E2;
%E3.T = F3.T + 10;
% Q_preheat2 = property(E2,'h','kJ') - property(E3,'h','kJ');
% cp2 = property(F3,'c','kJ/(kmol K)');
% F4.T = F3.T + Q_preheat2/(F3*cp2); 
% E4.T = F2.T;
% Q_removed = property(E1,'h','kJ') - property(E2,'h','kJ');
% Q_addtl_fuel_heat = 0;
% if Q_removed > Q_preheat
%     H_E2 = property(E1,'h','kJ') - Q_preheat;
%     E3.T = E1.T - Q_preheat/Q_removed*(E1.T - F3.T);
%     E3.T = find_T(E3, H_E3);
% else
%     Q_addtl_fuel_heat = Q_preheat - Q_removed;
% end
% HeatLoop.Q_preheat = Q_preheat; 


% HeatLoop.Qoutanode = property(E1,'h','kJ') - property(F4,'h','kJ'); %Heat Leaving FC by fuel exhaust
% HeatLoop.Qoutcathode = property(A4,'h','kJ') - property(A3,'h','kJ'); %Heat leaving fc by cathode exhaust
% HeatLoop.Qexcess = FCArray.Qgen - HeatLoop.Qoutanode - HeatLoop.Qoutcathode - HeatLoop.Qairin; 
% HeatLoop.H2used = zeros(10,10); 
% for i = 1:10
%     for j = 1:10
%         if HeatLoop.Qexcess(i,j) < 0
%             HeatLoop.H2used(i,j) = -HeatLoop.Qexcess(i,j)/(FCArray.hrxnmol(1,1)); 
%             HeatLoop.Qexcess(i,j) = 0;
%         end
%     end
% end
% HeatLoop.CoolingLoad = property(F2,'h','kJ') - property(F1,'h','kJ'); 
% HeatLoop.Qexchangetotal = HeatLoop.Qairin1 + HeatLoop.Qinfuel1 + HeatLoop.Qinfuel2;

% HeatLoop.F4 = exergy(F4,options.T0,options.P0);
% HeatLoop.E2 = exergy(E2,options.T0,options.P0);
% HeatLoop.F3 = exergy(F3,options.T0,options.P0); 
% HeatLoop.E1 = exergy(E1,options.T0,options.P0);
% HeatLoop.E3 = exergy(E3,options.T0,options.P0);
% HeatLoop.E4 = exergy(E4,options.T0,options.P0);
end