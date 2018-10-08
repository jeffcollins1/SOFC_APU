function param = NetParam(options,FC,turbo,OTM,HL)
Ein = FC.H2_used.*FC.hrxnmol; %Lower heating value of fuel intake, kJ/kmol
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
    param.fuel_for_OTM_preheat = -min(0,HL.FCQbalance)./FC.hrxnmol;
    param.P_perm = options.P_perm; 
    param.FTE = param.Eout./(Ein - min(0,HL.FCQbalance));
    param.Qbalance = HL.Qexcess;
    param.FCQout = FC.Qremove; 
end

param.FC_eff = FC.Power./Ein;
param.NetPower = param.Eout;
param.iden = FC.i_den; 
param.O2_used = FC.O2; 
param.T1_work = turbo{2}.work;
param.C1_work = turbo{1}.work; 
param.H2_used = FC.H2_used; 
param.FCVoltage = FC.V; 
param.FCPower = FC.Power; 
param.Qremovefuel = HL.Qremove_fuel;
param.Q_preheat = HL.Q_preheat; 
param.FCQgen = FC.Qgen;
param.FCV = FC.V;
param.hrxnmol = FC.hrxnmol; 
param.sofc_area = options.SOFC_area;
param.i_den = FC.i_den; 
end