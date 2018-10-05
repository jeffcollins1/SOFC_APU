function param = NetParam(options,FC,OTM,HL)
Ein = FC.H2_used.*FC.hrxnmol; %Lower heating value of fuel intake, kJ/kmol
param.Eout = FC.Power + OTM.net_work + HL.blower_work;
param.FTE = param.Eout./(Ein - min(0,HL.FCQbalance));
param.FC_eff = FC.Power./Ein;
param.NetPower = param.Eout;
param.iden = FC.i_den; 
param.O2_used = FC.O2; 
param.T1_work = OTM.T1_work;
param.C1_work = OTM.C1_work; 
param.BlowerWork = HL.blower_work; 
param.H2_used = FC.H2_used; 
param.FCVoltage = FC.V; 
param.FCPower = FC.Power; 
param.Qbalance = HL.Qexcess;
param.Qremovefuel = HL.Qremove_fuel;
param.Q_preheat = HL.Q_preheat; 
param.FCQgen = FC.Qgen;
param.i_total = FC.i_total;
param.OTM.net_work = OTM.net_work; 
param.OTM.heat_in = OTM.heat_added;
param.OTM.heat_out = OTM.Q_out;
param.OTM.Ra = OTM.Ra;
param.OTM.mean_flux = OTM.mean_flux;
param.fuel_for_OTM_preheat = -min(0,HL.FCQbalance)./FC.hrxnmol;
param.P_non_perm = options.P_non_perm; 
param.P_perm = options.P_perm; 
param.i_cell = 81*FC.i_den;
param.FCQout = FC.Qremove; 
param.FCV = FC.V;
param.hrxnmol = FC.hrxnmol; 
param.sofc_area = options.SOFC_area;
param.i_den = FC.i_den; 
end