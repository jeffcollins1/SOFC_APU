function [Design] = sizer_3_b(options,param,Weight,HL)

Qbal = HL.Qexcess; 
% % point = min(min(Qbal));
% [x,y] = find(Qbal ==point);
FTE = param.FTE;
point = max(max(FTE));
 [x,y] = find(FTE==point); 
% point = max(max(param.P_den)); 
% [x,y] = find(param.P_den ==point); 
batteryboost = 0; 
P_net = param.NetPower(x,y); 
scale = (1-batteryboost)*options.Pprop_to/P_net; 
% maximumEff = max(max(param.FTE));
% [x,y] = find(param.FTE == maximumEff);
Design.Qbal = scale*Qbal(x,y);
Design.H2_usednom = scale*param.H2_used(x,y)*2; %Fuel consumption by fuel cell, kg/s
Design.H2usedactual = scale*param.H2_used(x,y)*2; %Actual fuel burn at max cruise efficiency, kg/s

Design.point = point;

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
Design.MaxPower = scale*P_net; 
Design.Penalty = scale*Weight.Total(x,y) - 0.7*4*(9670/2.2); %Total Weight Penalty of option 2
Design.P_den = scale*P_net/Weight.Total(x,y); 


Design.preheat_air = scale*param.preheat_air(x,y);
Design.FCQgen = scale*param.Qgen(x,y);
Design.Qoutanode = scale*param.Qoutanode(x,y);
Design.Qoutcathode = scale*param.Qoutcathode(x,y); 
%Design.FCQremove = scale*param.FCQout(x,y);
Design.HXmass = scale*Weight.hx(x,y); 
Design.Cells = scale*param.Cells(x,y);
Design.scale = scale; 
Design.iden = param.iden(x,y);
Design.QoutMotor = Design.MaxPower/0.986 - Design.MaxPower; %Cooling load required by motor
Design.QbalanceMotor = scale*HL.CoolingLoad - Design.QoutMotor; 
batteryweight = batteryboost*options.Pprop_to*options.time_to/1260; %battery weight required to assist with takeoff
Design.fuelstorageactual = (101100 +  - Design.Penalty)*(1660/(1660+270)); %Total LH2 storage including weight of insulated container
%Design.Ebal = Design.H2usedactual* param.hrxnmol(1,1)/2 - Design.MaxPower - Design.OTMQout - Design.Qfuelout; 
%Design.FCQbalance = Design.FCQremove - Design.OTMQin; 
Design.specs = {'Net Power',options.PrCruise;'System Efficiency',Design.point;'FC Efficiency',Design.FCeff;'Current Density',Design.iden; 'Voltage',Design.Voltage;'FC Power',Design.SOFC_power;'Turbine Work',Design.T1_work;'Compressor Work',Design.C1_work; 'O2 utilization',Design.O2util;'H2 used',Design.H2usedactual;'Air in',Design.air_in;'Heat Balance',Design.Qbal;};
end