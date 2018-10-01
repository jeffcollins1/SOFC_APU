function [Design,Demand,Weight,param] = run_cycle(options,mission)
A1.T = 298*ones(10,10); %Ambient Temperature Kelvin;
A1.P = 100*ones(10,10); %Ambient Pressure as a function of altitude
A1.O2 = .21*options.airflow;
A1.N2 = .79*options.airflow;
[OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1);
[FC,E1] = fuelcell(options,O5);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,OTM,E1);
[Weight,param] = NetParam(options,FC,OTM,HL,A4);

param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
% [m,n] = size(param.states);
% for i = 1:1:m
%     eval(strcat(param.states{i,1}, ' = exergy(',param.states{i,1},',options.T0,options.P0);'));
% 
% end
% param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
% 
%  param.X_dest_intake = -OTM.net_work + OTM.heat_added.*(1- options.T0./E1.T) - OTM.Q_out.*(1-options.T0./O3.T)  - A5.X - O5.X;
%  param.X_dest_FC = OTM.O5.X + F3.X -FC.Power + FC.H2_used*235000 - (OTM.heat_added-OTM.Q_out).*(1- options.T0./OTM.A2T)-HL.Qexcess.*(1-options.T0./1073) - E2.X;
%  param.X_dest_recirc = -HL.blower_work -15850 + FC.Qremove.*(1-options.T0/350) + 19109 - 133.11;
 param.Qbalance = HL.Qexcess;
 param.Qremovefuel = HL.Qremove_fuel;
 param.Q_preheat = HL.Q_preheat; 
 param.FCPower = FC.Power; 
 param.FCQgen = FC.Qgen;
 param.i_total = FC.i_total;
 param.Xch_in = FC.H2_used.*(FC.G - FC.G0);
 param.H2_used = FC.H2_used; %H2 flow in kmol/s
 param.HRXNmol = FC.hrxnmol; 
 param.O2_used = FC.O2;
 param.OTMnetwork = OTM.net_work;
 param.blowerwork = HL.blower_work; 
 param.OTMheat_in = OTM.heat_added;
 param.OTMheat_out = OTM.Q_out;
 %param.OTMflux = OTM.flux; 
 param.P_non_perm = options.P_non_perm; 
 param.P_perm = options.P_perm; 
 param.i_cell = 81*FC.i_den;
 param.FCQout = FC.Qremove; 
 param.FCV = FC.V;
 param.hrxnmol = FC.hrxnmol; 
 param.sofc_area =options.SOFC_area;
 param.Eout = FC.Power + HL.blower_work + OTM.net_work;
 param.P_den = param.NetPower./param.weight;
 param.i_den = FC.i_den; 
[Demand,param.Ppropmax,param.Pprop_to,param.TrTO,param.time_to,param.velocityclimb,param.time_climb,param.Palt,param.Talt,param.dalt,param.Malt,param.PrC,param.altitude,param.dcruise,param.PCruise,param.Tcruise,param.Vcruise,param.t1,param.t2,param.t3,param.t4,param.t5,param.t6,param.t7,param.t8,param.Thr_to,param.Thr_cruise,clcruise,cdcruise] = craft(options);
% param.velocity = options.velocity_to;
[Design] = sizer(options,param,Weight);
% param.fb_to = 4*param.H2_used*param.time_to; %number of engines, fuel consumption per engine at max power, time that max power is consumed
end

%Demand,Ppropmax,P_to,TrTO,time_to,velocityclimb,time_climb,P,T,dclimb,M,PrC,altitude,dcruise,PCruise,Tcruise,Vcruise,t1,t2,t3,t4,t5,t6,t7,t8,Fnet_to,TrCruise
%param.mass,param.air_in,param.H2_used_to,param.FC_area_fixed,param.P_int_max,param.O2max,param.scale,param.MaxPower,param.FTE_to,param.batteryweight,param.FuelStorageActual,param.Pden,param.Qbal,param.H2usedactual