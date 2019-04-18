function param = NetParam(options,FC,turbo,OTM,FL)
F = 96485.33; %Faraday's Constant in Coulombs/mol
Ein = FC.i_total./(2000.*F).*FC.hrxnmol + FC.Q_pre_combustor; %Lower heating value of fuel intake, kJ/kmol
param.sofc_area = options.SOFC_area;
param.Eout = FC.Power;
for i = 1:1:length(turbo)
    param.Eout = param.Eout + turbo{i}.work;
end
param.FTE = param.Eout./Ein;
if ~isempty(OTM)
    param.OTM.net_work = turbo{1}.work+turbo{2}.work+turbo{4}.work; 
    param.OTM.heat_in = OTM.heat_added;
    param.OTM.heat_out = OTM.Q_out;
    param.OTM.Ra = OTM.Ra;
    param.OTM.mean_flux = OTM.mean_flux;
    param.fuel_for_OTM_preheat = -min(0,FC.Qremove - OTM.heat_added)./FC.hrxnmol;
    param.P_perm = options.P_perm; 
    param.FTE = param.Eout./(Ein  + param.fuel_for_OTM_preheat.*FC.hrxnmol);
    param.Qbalance = FC.Qremove - OTM.heat_added + OTM.Q_out + FL.Qremove_fuel;
    param.FCQout = FC.Qremove; 
end
param.FC_eff = FC.Power./Ein;

param.NetPower = param.Eout;
param.iden = FC.i_den; 
param.O2_used = FC.O2; 
param.T1_work = turbo{2}.work;
param.C1_work = turbo{1}.work; 
param.H2_used = FC.i_total./(2000.*F) + FC.Q_pre_combustor./FC.hrxnmol; 
param.FCVoltage = FC.V; 
param.FCPower = FC.Power; 
param.FCQgen = FC.Qgen;
param.FCV = FC.V;
param.hrxnmol = FC.hrxnmol; 
param.i_den = FC.i_den; 
param.Qremovefuel = FL.Qremove_fuel;
param.Q_preheat = FL.Q_preheat; 
end