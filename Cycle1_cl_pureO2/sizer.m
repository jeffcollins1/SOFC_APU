function [Design] = sizer(options,param,Weight)
P_den = param.P_den;
maximumPden = max(max(P_den));
[x,y] = find(P_den==maximumPden); 
% maximumEff = max(max(param.FTE));
% [x,y] = find(param.FTE == maximumEff);
P_net = param.NetPower(x,y); 
Design.P_airin = options.P_non_perm(x,y);
batteryboost = 0.0; 
scale = (1-batteryboost)*param.Pprop_to./P_net; 
Design.air_in = scale*options.airflow(x,y); 
Design.H2_usednom = scale*param.H2_used(x,y)*2; %Fuel consumption, kg/s
Design.mass = scale*param.weight(x,y); 
Design.Comp_mass = scale*Weight.comp(x,y);
Design.Turb_mass = scale*Weight.turb(x,y);
Design.SOFC_size = scale*param.sofc_area(x,y); 
Design.SOFC_mass = scale*Weight.sofc(x,y); 
Design.OTM_mass = scale*Weight.otm(x,y); 
Design.HTSM_mass = scale*Weight.motor(x,y); 
Design.HX_mass = scale*Weight.hx(x,y); 
Design.O2max = scale*param.O2_used(x,y); 
Design.MaxPower = scale*P_net; 
Design.Penalty = scale*Weight.Total(x,y) - 0.7*4*(9670/2.2); %Total Weight Penalty of option 1
Design.P_den = param.P_den(x,y); 
Design.FTE_to = param.FTE(x,y);
Design.Qbal = scale*param.Qbalance(x,y);
Design.FCQgen = scale*param.FCQgen(x,y);
Design.FCQremove = scale*param.FCQout(x,y);
Design.OTMQin =  scale*param.OTMheat_in(x,y);
Design.OTMQout = scale*param.OTMheat_out(x,y);
Design.Qfuelout =  scale*param.Qremovefuel(x,y);
Design.HXmass = scale*Weight.hx(x,y); 
Design.FCeff = param.FC_eff(x,y); 
Design.iden = param.iden(x,y);
Design.Voltage = param.FCVoltage(x,y); 
Design.SOFC_power = scale*param.FCPower(x,y);
Design.T1_work = scale*param.T1_work(x,y);
Design.C1_work = scale*param.C1_work(x,y); 
Design.BlowerWork = scale*param.BlowerWork(x,y); 
Design.scale = scale; 
H2_heat_to = 0
if Design.Qbal < 0
    H2_heat_to = -Design.Qbal/param.hrxnmol; 
end
Design.H2usedactual = Design.H2_usednom + H2_heat_to*2; %Actual fuel burn at max takeoff power, kg/s
batteryweight = batteryboost*param.Pprop_to*param.time_to/1260; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
Design.FuelMassBalance = 101000 - Design.Penalty; 
Design.fuelstorageactual = Design.FuelMassBalance/1.15; 
Design.payload = 112760 + Design.FuelMassBalance - Design.fuelstorageactual; 
Design.Ebal = Design.H2usedactual* param.HRXNmol(1,1)/2 - Design.MaxPower - Design.OTMQout - Design.Qfuelout; 
Design.FCQbalance = Design.FCQremove - Design.OTMQin; 
Design.specs = {'Net Power',param.Pprop_to;'System Efficiency',Design.FTE_to;'FC Efficiency',Design.FCeff;'Current Density',Design.iden; 'Voltage',Design.Voltage;'FC Power',Design.SOFC_power;'Turbine Work',Design.T1_work;'Compressor Work',Design.C1_work;'Recirculation Work',Design.BlowerWork; 'O2 utilization',Design.O2max ;'H2 used',Design.H2usedactual;'Air in',Design.air_in;'Heat Balance',Design.Qbal;};
end
