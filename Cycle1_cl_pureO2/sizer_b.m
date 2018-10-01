function [Design] = sizer_b(options,param,Demand,P_non_perm,A2,A6)
%P_den = param.P_den;
% maximumPden = max(max(P_den));
% [x,y] = find(P_den==maximumPden); 
H2_heat_to = zeros(10,10);
for i = 1:10
    for j = 1:10
if param.Qbalance(i,j) < 0
    H2_heat_to(i,j) = -param.Qbalance(i,j)/param.hrxnmol(1,1);  %Check to see if additional preheating is needed
    param.Qbalance(i,j) = 0; 
end
    end
end

H2heat = H2_heat_to; 
FTEreal = param.NetPower./((param.H2used + H2heat).*param.HRXNmol); 

% for r = 1:10
%     for s=1:10
%         if param.iden(r,s) > 1.5 || param.iden(r,s) < 1.30
%             FTEreal(r,s) = 0;
%         
%   
%         end
%         end
% end

maximumEff = max(max(FTEreal)); %Size for max eff
[x,y] = find(FTEreal == maximumEff);
% minSOFC = min(min(param.sofc_area)); %size for min sofc area within a range of current densities
% [x,y] = find(param.sofc_area==minSOFC);
%  l = w(1,1);
% FTE = FTEreal(l,z); 
% maxFTE = max(FTE);
% [x,y] = find(FTEreal == maxFTE);
P_net = param.NetPower(x,y); 
Design.P_airin = P_non_perm;
batteryboost = 0.0; 
scale = (1-batteryboost)*Demand.PrCruise./P_net; 
Design.P_perm = options.P_perm;
Design.P_non_perm = P_non_perm; 
Design.air_in = scale*options.airflow(x,y); 
Design.H2_usednom = scale*param.H2_used(x,y)*2; %Fuel consumption, kg/s
Design.H2_usedheat = scale.*H2heat*2; 
%Design.mass = scale*param.weight(x,y); 
% Design.Comp_mass = scale*Weight.comp(x,y);
% Design.Turb_mass = scale*Weight.turb(x,y);
Design.SOFC_size = scale*param.sofc_area(x,y); 
% Design.SOFC_mass = scale*Weight.sofc(x,y); 
% Design.OTM_mass = scale*Weight.otm(x,y); 
% Design.HTSM_mass = scale*Weight.motor(x,y); 
% Design.HX_mass = scale*Weight.hx(x,y); 
Design.O2max = scale*param.O2_used(x,y); 
Design.MaxPower = scale*P_net; 
%Design.Penalty = scale*Weight.Total(x,y) - 0.7*4*(9670/2.2); %Total Weight Penalty of option 1
%Design.P_den = param.P_den(x,y); 
Design.FTE = FTEreal(x,y);
Design.Qbal = scale*param.Qbalance(x,y);
Design.FCQgen = scale*param.FCQgen(x,y);
Design.FCQremove = scale*param.FCQout(x,y);
Design.OTMQin =  scale*param.OTMheat_in(x,y);
Design.OTMQout = scale*param.OTMheat_out(x,y);
Design.Qfuelout =  scale*param.Qremovefuel(x,y);
%Design.HXmass = scale*Weight.hx(x,y); 
Design.FCeff = param.FC_eff(x,y); 
Design.iden = param.iden(x,y);
Design.Voltage = param.FCVoltage(x,y); 
Design.SOFC_power = scale*param.FCPower(x,y);
Design.T1_work = scale*param.T1_work(x,y);
Design.C1_work = scale*param.C1_work(x,y); 
Design.BlowerWork = scale*param.BlowerWork(x,y); 
Design.scale = scale; 



 

Design.H2usedactual = Design.H2_usednom + Design.H2_usedheat; %Actual fuel burn at max takeoff power, kg/s
%batteryweight = batteryboost*Demand.Pr_to*param.time_to/1260; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;
%Design.FuelMassBalance = 101000 - Design.Penalty; 
%Design.fuelstorageactual = Design.FuelMassBalance/1.15; 
%Design.payload = 112760 + Design.FuelMassBalance - Design.fuelstorageactual; 
Design.Ebal = Design.H2usedactual* param.HRXNmol(1,1)/2 - Design.MaxPower - Design.OTMQout - Design.Qfuelout; 
Design.FCQbalance = Design.FCQremove - Design.OTMQin; 
Design.specs = {'Net Power',Design.MaxPower;'System Efficiency',Design.FTE;'FC Efficiency',Design.FCeff;'Current Density',Design.iden; 'Voltage',Design.Voltage;'FC Power',Design.SOFC_power;'Turbine Work',Design.T1_work;'Compressor Work',Design.C1_work;'Recirculation Work',Design.BlowerWork; 'O2 utilization',Design.O2max ;'H2 used',Design.H2usedactual;'Air in',Design.air_in;'Heat Balance',Design.Qbal;};
end
