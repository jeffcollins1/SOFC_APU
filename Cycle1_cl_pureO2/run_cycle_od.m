function param_od = run_cycle_od(OD,options)
A1.T = OD.Talt; %Ambient Temperature Kelvin;
A1.P = OD.Palt; %Ambient Pressure as a function of altitude
A1.O2 = .21*OD.air_in;
A1.N2 = .79*OD.air_in;
[OTM_od,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule_od(OD,A1);
for i = 1:10
    for j = 1:10
if OTM_od.Rt(i,j) < 0
   O5.O2(i,j) = NaN;
   O5.T(i,j) = NaN;
   O5.P(i,j) = NaN;
end
    end
end
[FC_od,E1] = fuelcell_od(OD,O5);
[HL_od,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(OD,FC_od,OTM_od,E1);
[param_od] = NetParam_od(OD,FC_od,OTM_od,HL_od);

param_od.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
[m,n] = size(param_od.states);
for i = 1:1:m
    eval(strcat(param_od.states{i,1}, ' = exergy(',param_od.states{i,1},',options.T0,options.P0);'));

end
param_od.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
% 
%  param_od.X_dest_intake = -OTM_od.net_work + OTM_od.heat_added.*(1- OD.T0./E1.T) - OTM_od.Q_out.*(1-OD.T0./O3.T)  - A5.X - O5.X;
%  param_od.X_dest_FC = OTM_od.O5.X + F3.X -FC_od.Power + FC_od.H2_used*235000 - (OTM_od.heat_added-OTM_od.Q_out).*(1- OD.T0./OTM_od.A2T)-HL_od.Qexcess.*(1-OD.T0./1073) - E2.X;
%  param_od.X_dest_recirc = -HL_od.blower_work -15850 + FC_od.Qremove.*(1-OD.T0/350) + 19109 - 133.11;
 param_od.Qbalance = HL_od.Qexcess;
 param_od.Qremovefuel = HL_od.Qremove_fuel;
 param_od.Q_preheat = HL_od.Q_preheat; 
 param_od.FCPower = FC_od.Power; 
 param_od.FCQgen = FC_od.Qgen;
 param_od.i_total = FC_od.i_total;
 param_od.Xch_in = FC_od.H2_used.*(FC_od.G - FC_od.G0);
 param_od.H2_used = FC_od.H2_used; %H2 flow in kmol/s
 param_od.HRXNmol = FC_od.hrxnmol; 
 param_od.O2_used = FC_od.O2;
 param_od.OTMnetwork = OTM_od.net_work;
 param_od.blowerwork = HL_od.blower_work; 
 param_od.OTMheat_in = OTM_od.heat_added;
 param_od.C1_work = OTM_od.C1_work;
 param_od.C2_work = OTM_od.C2_work;
 
 param_od.OTMheat_out = OTM_od.Q_out;
 %param_od.OTMflux = OTM_od.flux; 
 param_od.P_non_perm = OD.P_non_perm; 
 param_od.P_perm = OD.P_perm; 
 param_od.i_cell = 100*FC_od.i_den;
 param_od.FCQout = FC_od.Qremove; 
 param_od.FCV = FC_od.V;
 param_od.hrxnmol = FC_od.hrxnmol; 
 param_od.sofc_area =OD.SOFC_area;
 param_od.Eout = FC_od.Power + HL_od.blower_work + OTM_od.net_work;
 param_od.i_den = FC_od.i_den; 
 param_od.air_in = OD.air_in; 
 param_od.QexcessFC = FC_od.Qremove;
 param_od.Ebalcheck = FC_od.H2_used.*FC_od.hrxnmol - param_od.Qbalance - param_od.Eout;

% 
% param_od.fb_climb = mean(param_od.H2_used)*param_od.time_climb; 
end