function [intake,A2,A3,A4,A5] = intake_cycle_odb(options,FCArray)
Intake.O2 = 0.21*FCArray.airin;
Intake.N2 = 0.79*FCArray.airin; 
Intake.T = FCArray.Tamb;
Intake.P = FCArray.Pamb; 

[A2,C1_work] = compressor(Intake,options.P_fc,options.C1_eff);
A3 = A2;
A3.T = options.T_fc;
A3.P = options.P_fc; 
intake.heat_added =  property(A3,'h','kJ') - property(A2,'h','kJ');
intake.Y_O2 = A3.O2./net_flow(A3);
intake.Y_N2 = A3.N2./net_flow(A3); 
intake.Pratio = 1000./Intake.P; 
intake.C1_work = C1_work;
A4.O2 = FCArray.O2out;
A4.N2 = FCArray.N2out;
A4.T = A3.T+ 50*ones(10,10); 
A4.P = 980*ones(10,10); 
P_out = FCArray.P0; 
[A5,intake.T1_work] = expander(A4,P_out,options.T1_eff);
end