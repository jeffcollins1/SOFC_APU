function param = run_cycle_2(options,mission)
A1.T = -0.0065*options.height + 14.987 + 273.1; %Ambient Temperature Kelvin;
A1.P = 107*exp(-0.0001*options.height)- 10; %Ambient Pressure as a function of altitude
A1.O2 = .21*options.airflow;
A1.N2 = .79*options.airflow;
[FCair,E1,A1] = fuelcell_air(options,Cathode);
[intake,A2,A3,A4,A5] = intake_cycle(options,FCair);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FCair,intake,E1);
param = NetParam(options,FCair,intake,HL);

param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;};
[m,n] = size(param.states);
for i = 1:1:m
%     eval(strcat(param.states{i,1}, ' = exergy(',param.states{i,1},',options.T0,options.P0)'));
   param.output.(param.states{i,1}) = exergy(eval(param.states{i,1}), options.T0, options.P0);
   Tin = getfield(param.states{6,2},'T') %Temperature of heat added by FC Exhaust;
   Tout = getfield(param.states{16,2},'T') %Temperature of heat rejected before 02 recompression
   Xout1 = 1.321161660208921e+03;
   %Xout1 = struct(param(4,1).field('A5',9);('.getfield(param.output{'A5','X'}); %Exergy of mass leaving turbine
   Xout2 = 2.284421581676992e+03;
   %Xout2 = getfield(param.output{18,2},'X'); %Exergy of oxygen leaving intake system to FC system
 param.X_dest_intake = -intake.net_work + intake.heat_added*(1- options.T0/Tin) - intake.Q_out*(1-options.T0/Tout)  - Xout1 - Xout2;
 param.X_dest_FC = intake.O5.X + HL.F3.X -FCair.Power + FCair.H2_used*235000 - intake.heat_added*(1- options.T0/intake.A2T)-HL.Qexcess*(1-options.T0/1073) - HL.E2.X;
 param.X_dest_recirc = -HL.blower_work -15850 + FCair.Qremove*(1-options.T0/350) + 19109 - 133.11;
 param.Qbalance = HL.Qexcess;
 param.Q_preheat = HL.Q_preheat; 
 param.FCPower = FCair.Power; 
 param.FCQgen = FCair.Qgen;
 param.Xch_in = FCair.H2_used*(FCair.G - FCair.G0);
 param.H2_used = FCair.H2_used; %H2 flow in mol/s
 param.HRXNmol = FCair.hrxnmol; 
 param.O2_used = FCair.O2;
 param.OTMnetwork = intake.net_work;
 param.blowerwork = HL.blower_work; 
end


