function [intake,A2,A3] = intake_cycle_3(options,FCArray,A1)
Intake.O2 = 0.21*ones(10,10);
Intake.N2 = 0.79*ones(10,10); 
Intake.T = A1.T*ones(10,10);
Intake.P = A1.P*ones(10,10); 

[A2,C1_work] = compressor(Intake,options.P_fc,options.C1_eff);
A3 = A2;
A3.T = options.T_fc;
A3.P = options.P_fc; 
intake.heat_added =  property(A3,'h','kJ') - property(A2,'h','kJ');
intake.Y_O2 = A3.O2./net_flow(A3);
intake.Y_N2 = A3.N2./net_flow(A3); 

intake.C1_work = C1_work;


end
