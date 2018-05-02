function param = run_cycle(options,mission)
A1.T = -0.0065*options.height + 14.987 + 273.1; %Ambient Temperature Kelvin;
A1.P = 107*exp(-0.0001*options.height)- 10; %Ambient Pressure as a function of altitude
A1.O2 = .21*options.airflow;
A1.N2 = .79*options.airflow;
[OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1);
[FC,E1] = fuelcell(options,O5);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,OTM,E1);
param = NetParam(options,FC,OTM,HL);

param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
[m,n] = size(param.states);
for i = 1:1:m
    eval(strcat(param.states{i,1}, ' = exergy(',param.states{i,1},',options.T0,options.P0);'));
%    param.states{i,2} = exergy(param.states{i,2}, options.T0, options.P0);
%    Tin = getfield(param.states{6,2},'T'); %Temperature of heat added by FC Exhaust;
%    Tout = getfield(param.states{16,2},'T'); %Temperature of heat rejected before 02 recompression
%    Xout1 = 1.321161660208921e+03;
%    %Xout1 = struct(param(4,1).field('A5',9);('.getfield(param.output{'A5','X'}); %Exergy of mass leaving turbine
%    Xout2 = 2.284421581676992e+03;
%    %Xout2 = getfield(param.output{18,2},'X'); %Exergy of oxygen leaving intake system to FC system
end
param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};

 param.X_dest_intake = -OTM.net_work + OTM.heat_added*(1- options.T0/E1.T) - OTM.Q_out*(1-options.T0/O3.T)  - A5.X - O5.X;
 param.X_dest_FC = OTM.O5.X + F3.X -FC.Power + FC.H2_used*235000 - OTM.heat_added*(1- options.T0/OTM.A2T)-HL.Qexcess*(1-options.T0/1073) - E2.X;
 param.X_dest_recirc = -HL.blower_work -15850 + FC.Qremove*(1-options.T0/350) + 19109 - 133.11;
 param.Qbalance = HL.Qexcess;
 param.Q_preheat = HL.Q_preheat; 
 param.FCPower = FC.Power; 
 param.FCQgen = FC.Qgen;
 param.Xch_in = FC.H2_used*(FC.G - FC.G0);
 param.H2_used = FC.H2_used; %H2 flow in mol/s
 param.HRXNmol = FC.hrxnmol; 
 param.O2_used = FC.O2;
 param.OTMnetwork = OTM.net_work;
 param.blowerwork = HL.blower_work; 
 param.OTMheat_in = OTM.heat_added - OTM.Q_out;
end


