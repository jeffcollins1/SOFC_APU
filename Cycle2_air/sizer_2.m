function [Design] = sizer_2(options,param,Weight,HL)
% for i = 1:10
%     for j = 1:10
%         if param.FCVoltage(i,j) < 0.7
%             param.FCVoltage(i,j) = 0;
%             HL.Qexcess(i,j) = 1e7; 
%             param.FTE(i,j) = 0;
%             param.P_den(i,j) = 0;
%         end
%     end
% end
% point = max(max(param.FTE));
% [x,y] = find(param.FTE == point);
% Qbalance = abs(HL.Qexcess);
% point = min(min(Qbalance));
% [x,y] = find(Qbalance==point); 
Qbal = HL.Qexcess;
H2_heat_to = zeros(10,10);
for i = 1:10
    for j =1:10
if Qbal(i,j) < 0
    H2_heat_to(i,j) = -Qbal(i,j)/param.hrxnmol(1,1); 
    Qbal(i,j) = 0; 
end
    end
end
% Qbalance = abs(Qbal);
% point = min(min(Qbalance));
% [x,y] = find(Qbalance==point); 

% point = max(max(param.P_den))
% [x,y] = find(param.P_den ==point);

P_net = param.NetPower; 
Design.H2_usednom = param.H2_used*2; %Fuel consumption, kg/s
Design.H2usedtotal = Design.H2_usednom + 2*H2_heat_to; %Actual fuel burn at max takeoff power, kg/s
Design.FTE = P_net./(Design.H2usedtotal.*param.hrxnmol/2); 
point = max(max(Design.FTE));
[x,y] = find(Design.FTE == point);

batteryboost = 0; 
Pactual = P_net(x,y);
scale = (1-batteryboost)*options.Pprop_to./Pactual; 
Design.H2usedactual = scale*Design.H2usedtotal(x,y); 
Design.Qbal = scale*Qbal(x,y);
Design.point = point;
Design.H2used = Design.H2usedactual;
Design.air_in = scale*options.airflow(x,y); 
Design.T1_work = scale*param.T1_work(x,y);
Design.C1_work = scale*param.C1_work(x,y); 

Design.mass = scale*param.weight(x,y); 
Design.O2used = scale*param.O2used(x,y);
Design.FCeff = param.FC_eff(x,y); 
Design.Voltage = param.FCVoltage(x,y);
Design.O2util = param.O2util(x,y);
Design.Comp_mass = scale*Weight.comp(x,y);
Design.Turb_mass = scale*Weight.turb(x,y);
%Design.SOFC_size = scale*param.sofc_area(x,y); 
Design.SOFC_mass = scale*Weight.sofc(x,y); 
Design.SOFC_power = scale*param.FCPower(x,y);
Design.HTSM_mass = scale*Weight.motor(x,y); 
Design.HX_mass = scale*Weight.hx(x,y); 
%Design.O2max = scale*param.O2_used(x,y); 
Design.MaxPower = scale*Pactual; 
Design.Penalty = Design.mass - 0.7*4*(9670/2.2); %Total Weight Penalty of option 2
Design.P_den = scale*Pactual/Weight.Total(x,y); 
Design.FTE_to = scale*Pactual/(Design.H2usedactual*param.hrxnmol(1,1)/2);

Design.preheat_air = scale*param.preheat_air(x,y);
Design.FCQgen = scale*param.Qgen(x,y);
Design.Qoutanode = scale*param.Qoutanode(x,y);
Design.Qoutcathode = scale*param.Qoutcathode(x,y); 
%Design.FCQremove = scale*param.FCQout(x,y);
Design.Qfuelout =  scale*param.Qremovefuel(x,y);
Design.HXmass = scale*Weight.hx(x,y); 
Design.Cells = scale*param.Cells(x,y);
Design.scale = scale; 
Design.iden = param.iden(x,y);
Design.BlowerWork = scale*param.BlowerWork(x,y);
H2_heat_to = 0


batteryweight = batteryboost*options.Pprop_to*options.time_to/1260; %battery weight required to assist with takeoff
Design.fuelstorageactual = (101100 +  - Design.Penalty)*(1660/(1660+270)); %Total LH2 storage including weight of insulated container
%Design.Ebal = Design.H2usedactual* param.hrxnmol(1,1)/2 - Design.MaxPower - Design.OTMQout - Design.Qfuelout; 
%Design.FCQbalance = Design.FCQremove - Design.OTMQin; 
Design.specs = {'Net Power',options.Pprop_to;'System Efficiency',Design.FTE_to;'FC Efficiency',Design.FCeff;'Current Density',Design.iden; 'Voltage',Design.Voltage;'FC Power',Design.SOFC_power;'Turbine Work',Design.T1_work;'Compressor Work',Design.C1_work;'Recirculation Work',Design.BlowerWork; 'O2 utilization',Design.O2util;'H2 used',Design.H2usedactual;'Air in',Design.air_in;'Heat Balance',Design.Qbal;};
end
