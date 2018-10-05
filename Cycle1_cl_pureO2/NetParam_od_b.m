function [param] = NetParam_od_b(OD,FC,OTM,HL,A4)
Ein = FC.H2_used.*FC.hrxnmol; %Lower heating value of fuel intake, kJ/kmol
Eout = FC.Power + OTM.net_work + HL.blower_work;
param.FTE = Eout./Ein;
param.FC_eff = FC.Power./Ein;
param.NetPower = Eout;
param.iden = FC.i_den; 
param.O2_used = FC.O2; 
param.T1_work = OTM.T1_work;
param.C1_work = OTM.C1_work; 
param.BlowerWork = HL.blower_work; 
param.HXtotal = OTM.Q_out + OTM.heat_added + OTM.Q_oxygen_HX + HL.Q_preheat + HL.Qremove_fuel;% + OTM.Qrecover;
%[Weight] = weight_b(OD,param,FC,OTM,HL,A4);
%[Weight] = weight(options,param,FC,OTM,HL,A4);
% param.weight = Weight.Total;
% param.weightOTM = Weight.otm;
% param.weightComp = Weight.comp;
% param.weightTurb = Weight.turb;
% param.weightFC = Weight.sofc; 
% param.weightHX = Weight.hx; 
% param.fuelstoragenominal = Weight.Fuel; 
param.FCVoltage = FC.V; 
param.FCPower = FC.Power; 
end