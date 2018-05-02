function [intake,A2,A3,A4,A5] = intake_cycle(options,A1)

[A2,C1_work] = compressor(A1,options.P_non_perm,options.C1_eff);
A3 = A2;
A3.T = options.T_otm;
intake.heat_added =  property(A3,'h','kJ') - property(A2,'h','kJ');
X_O2 = A3.O2/net_flow(A3);



intake.A2T = A2.T;
A4 = A3;
[A5,T1_work] = expander(A4,A1.P,options.T1_eff);
intake.C1_work = C1_work;
intake.T1_work = T1_work;
intake.work_in = C1_work;
intake.net_work = C1_work + T1_work;
intake.A5 = exergy(A5,options.T0,options.P0);

