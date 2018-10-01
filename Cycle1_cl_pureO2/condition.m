function [alt1,OD] = condition(OD,param_od,z,A2,A6)
P_req = OD.altdata(z,6);
P_gen = param_od.Eout;
error = abs(P_req - P_gen);
min_err = min(min(error));
Qbal = param_od.Qbalance; 
H2_used = param_od.H2_used;
%H2_heat=param_od.H2heat;
hrxnmol = param_od.HRXNmol;
for a = 1:10
    for b = 1:10
        if Qbal(a,b) < 0
            H2_heat(a,b) = -Qbal(a,b)/hrxnmol(a,b); %
            Qbal(a,b) = 0; 
        end
    end
end
OD.Qbal = Qbal; 
H2_total = H2_heat*2 + 2.*H2_used;
OD.H2_total = H2_total;
FTEtrue = param_od.Eout./(OD.H2_total.*hrxnmol/2);
OD.FTEtrue = FTEtrue; 
OD.Voltage = param_od.FCVoltage;
OD.iden = param_od.iden; 
OD.T1_work = param_od.T1_work;
OD.C1_work = param_od.C1_work;
OD.BlowerWork = param_od.BlowerWork; 
OD.O2used = param_od.O2_used; 
OD.FCeff = param_od.FC_eff;
OD.FCPower = param_od.FCPower; 
OD.Qhxtotal = param_od.HXtotal;  
OD.Qbalmotors = param_od.Qbalmotor;
selection = zeros(18,17);
for c = 1:9
    for d = 1:10
        if P_req(1,1) <= param_od.Eout(c,d) && P_req(1,1) >= param_od.Eout((c+1),d)
            fraction = (P_req(1,1) - param_od.Eout((c+1),d))/(param_od.Eout(c,d) - param_od.Eout((c+1),d));
            selection(d,1) = P_req(1,1);
            selection(d,2) = OD.P_perm((c+1),d) + fraction*(OD.P_perm(c,d) - OD.P_perm((c+1),d)); 
            selection(d,3) = OD.FTEtrue((c+1),d) + fraction*(OD.FTEtrue(c,d) - OD.FTEtrue((c+1),d));
            selection(d,4) = OD.air_in((c+1),d) + fraction*(OD.air_in(c,d) - OD.air_in((c+1),d));
            selection(d,5) = OD.H2_total((c+1),d) + fraction*(OD.H2_total(c,d) - OD.H2_total((c+1),d));
            selection(d,6) = OD.Qbal((c+1),d) + fraction*(OD.Qbal(c,d) - OD.Qbal((c+1),d));
            selection(d,7) = H2_heat((c+1),d) + fraction*(H2_heat(c,d) - H2_heat((c+1),d));
            selection(d,8) = OD.FCeff((c+1),d) + fraction*(OD.FCeff(c,d) - OD.FCeff((c+1),d));
            selection(d,9) = OD.Voltage((c+1),d) + fraction*(OD.Voltage(c,d) - OD.Voltage((c+1),d));
            selection(d,10) = OD.iden((c+1),d) + fraction*(OD.iden(c,d) - OD.iden((c+1),d));
            selection(d,11) = OD.T1_work((c+1),d) + fraction*(OD.T1_work(c,d) - OD.T1_work((c+1),d));
            selection(d,12) = OD.C1_work((c+1),d) + fraction*(OD.C1_work(c,d) - OD.C1_work((c+1),d));
             selection(d,13) = OD.BlowerWork((c+1),d) + fraction*(OD.BlowerWork(c,d) - OD.BlowerWork((c+1),d));
             selection(d,14) = OD.O2used((c+1),d) + fraction*(OD.O2used(c,d) - OD.O2used((c+1),d));
              selection(d,15) = OD.FCPower((c+1),d) + fraction*(OD.FCPower(c,d) - OD.FCPower((c+1),d));
              selection(d,16) = OD.Qhxtotal((c+1),d) + fraction*(OD.Qhxtotal(c,d) - OD.Qhxtotal((c+1),d));
              selection(d,17) = 0.79*((OD.air_in((c+1),d) + fraction*(OD.air_in(c,d) - OD.air_in((c+1),d))))*28 + 0.21*(OD.air_in((c+1),d) + fraction*(OD.air_in(c,d) - OD.air_in((c+1),d)) - OD.O2used((c+1),d) + fraction*(OD.O2used(c,d) - OD.O2used((c+1),d)))*32; %Mass flow into turbine
                  
              selection(d,18) = OD.Qbalmotors((c+1),d) + fraction*(OD.Qbalmotors(c,d) - OD.Qbalmotors((c+1),d)); 
                  
           end
        end
     end
              
maxEff = max(selection(:,3)); 
[x] = find(selection(:,3)==maxEff);
if size(x) == 1
alt1(z,1:18) = selection(x,1:18); 
OD.alt1(z,1:18) = alt1(z,1:18); 

else    
alt1(z,1:18) = NaN; 
OD.alt1(z,1:18) = NaN; 
end
      
%FC_od,OTM_od,HL_od
%  P = P_gen;
%  X = OD.P_perm;
%  Y = OD.air_in; 
%  Z = OD.FTEtrue; 
%  p = P_req(1,1);
% [x,y,z] = constrained_max_val(X,Y,Z,P,p);

% % err_tol = min_err + 1000;
% alt1.P_req = P_req; 
% % [x,y] = find(error<=err_tol); 
% % Effrange = OD.FTEtrue(x,y);
% % Effmax = max(max(abs(Effrange)));
% [a,b] = find(OD.P_perm==x);
% [c,d] = find(OD.air_in ==y);
% g = a(1,1);
% h = d(1,1);
% %alt1.P_err = error(x,y);
% %alt1.P_err_percent =(param_od.Eout(g,h)- alt1.P_err)/param_od.Eout(g,h);
% alt1.P_act = param_od.Eout(g,h); 
% alt1.air_in = y;
% alt1.H2used = OD.H2_total(g,h);  
% alt1.Qbalance = Qbal(g,h); 
% alt1.Eff = OD.FTEtrue(g,h); 
% alt1.P_perm = x; 
% alt1.H2_heat = H2_heat(g,h);
end