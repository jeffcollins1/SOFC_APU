function param = NetParam(options,FC,OTM,HL)
Ein = FC.H2_used*FC.hrxnmol; %Lower heating value of fuel intake, kJ/kmol
Eout = FC.Power + OTM.net_work + HL.blower_work;
param.FTE = Eout/Ein;
param.FC_eff = FC.Power/Ein;
%param.subsystem1 = [OTM.heat_added,OTM.Q_out,OTM.OTM.work_in,OTM.T1_work,1.321161660208921e+03,2.284421581676992e+03]';

end