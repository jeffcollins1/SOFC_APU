function [OffDesign,Performance,Comp,Turb] = condition_3_exp(P,paramod,FCArray,y,OffDesign,Weight)
P_req = P;
P_gen = paramod.NetPower;
error = abs(P_req - P_gen);
min_err = min(min(error));
Qbal = paramod.Qbalance; 
H2_used = paramod.H2_used*2;
OffDesign.Qbal = Qbal; 
OffDesign.H2_total = H2_used;
FTEtrue = paramod.FTE;
OffDesign.FTEtrue = FTEtrue; 
OffDesign.air_in = FCArray.airin;
OffDesign.iden = FCArray.iDenArray; 
OffDesign.Voltage = FCArray.FCVoltage;
OffDesign.FCeff = FCArray.Efficiency;
OffDesign.O2utilization = FCArray.O2util;
OffDesign.Qairpreheat = paramod.preheat_air; 
OffDesign.Comp1Work = paramod.C1_work;
OffDesign.Turb1Work = paramod.T1_work;
OffDesign.FCPower = FCArray.Power;
selection = zeros(18,16);
for i = 1:10
    for j = 1:10
        if FCArray.FCVoltage(i,j) < 0.65
            FCArray.FCVoltage(i,j) = 0;
            FCArray.Power(i,j) = 0;
            paramod.NetPower(i,j) = 0; 
        end
    end
end
% maxP = max(max(paramod.NetPower));
% [y,z] = find(paramod.NetPower ==maxP);
% Pshort = P - paramod.NetPower(y,z);
% Pshortturbine = Pshort - OffDesign.Turb1Work(y,z);
% PtotalTurbine = Pshort + OffDesign.Turb1Work(y,z); 

for c = 1:9
    for d = 1:10
        if P_req <= P_gen(c,d) && P_req >= P_gen((c+1),d)
            fraction = -(P_req - paramod.NetPower((c+1),d))/(-(paramod.NetPower(c,d) - paramod.NetPower((c+1),d)));
            selection(d,1) = P_req(1,1);
            selection(d,2) = OffDesign.FTEtrue((c+1),d) + fraction*(OffDesign.FTEtrue(c,d) - OffDesign.FTEtrue((c+1),d));
            selection(d,3) = OffDesign.air_in((c+1),d) + fraction*(OffDesign.air_in(c,d) - OffDesign.air_in((c+1),d));
            selection(d,4) = OffDesign.H2_total((c+1),d) + fraction*(OffDesign.H2_total(c,d) - OffDesign.H2_total((c+1),d));
            selection(d,5) = OffDesign.Qbal((c+1),d) + fraction*(OffDesign.Qbal(c,d) - OffDesign.Qbal((c+1),d));
            %selection(d,6) = OffDesign.H2_heat((c+1),d) + fraction*(OffDesign.H2_heat(c,d) - OffDesign.H2_heat((c+1),d));
            selection(d,7) = OffDesign.iden((c+1),d) + fraction*(OffDesign.iden(c,d) - OffDesign.iden((c+1),d)); 
             selection(d,8) = OffDesign.Qairpreheat((c+1),d) + fraction*(OffDesign.Qairpreheat(c,d) - OffDesign.Qairpreheat((c+1),d)); 
              selection(d,9) = OffDesign.Comp1Work((c+1),d) + fraction*(OffDesign.Comp1Work(c,d) - OffDesign.Comp1Work((c+1),d)); 
               selection(d,10) = OffDesign.Turb1Work((c+1),d) + fraction*(OffDesign.Turb1Work(c,d) - OffDesign.Turb1Work((c+1),d)); 
                selection(d,11) = OffDesign.FCeff((c+1),d) + fraction*(OffDesign.FCeff(c,d) - OffDesign.FCeff((c+1),d)); 
                 selection(d,12) = OffDesign.O2utilization((c+1),d) + fraction*(OffDesign.O2utilization(c,d) - OffDesign.O2utilization((c+1),d)); 
                  selection(d,13) = OffDesign.Voltage((c+1),d) + fraction*(OffDesign.Voltage(c,d) - OffDesign.Voltage((c+1),d)); 
                   selection(d,14) = OffDesign.FCPower((c+1),d) + fraction*(OffDesign.FCPower(c,d) - OffDesign.FCPower((c+1),d));
                     selection(d,15) = Weight.comp((c+1),d) + fraction*(Weight.comp(c,d) - Weight.comp((c+1),d)); 
                     selection(d,16) = Weight.turb((c+1),d) + fraction*(Weight.turb(c,d) - Weight.turb((c+1),d)); 
        end
    end
end
              
maxEff = max(selection(:,2)); 
[x] = find(selection(:,2)==maxEff);
Performance(y,1:15) = selection(x,1:15);
OffDesign.H2usedactual = selection(x,4); 
Comp = selection(y,16);
Turb = selection(y,17); 
OffDesign.specs = zeros(15,2); 
OffDesign.specs = {'Net Power',selection(x,1);'System Efficiency',selection(x,2);'FC Efficiency',selection(x,11);'Current Density',selection(x,7); 'Voltage',selection(x,13);'FC Power',selection(x,14);'Turbine Work',selection(x,10);'Compressor Work',selection(x,9);'Recirculation Work',selection(x,15); 'O2 utilization',selection(x,12);'H2 used',selection(x,4);'Air in',selection(x,3);'Heat Balance',selection(x,5);};
end 