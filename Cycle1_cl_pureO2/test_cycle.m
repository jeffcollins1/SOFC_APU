%%  Test Cycle
%Create options to test
%Initial sizing to meet max power demand at sea level
n1 = 10; % number of points in test dimension 1
n2 = 10; % number of points in test dimension 1
options.height = 10*ones(n1,n2); %Altitude, meters
% options.OTM_area = 2.5e3*ones(n1,n2); %membrane area in m^2
options.SOFC_area = (ones(n1,1)*linspace(1800,3000,n2))';%4e3*ones(n1,n2); 
options.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
options.asr = 0.15*ones(n1,n2); % Area specific resistance, ohm-cm^2
options.P_fc = 1000*ones(n1,n2); %Operating pressure for SOFC
options.T_fc = 1023*ones(n1,n2); %Inlet temperature for SOFC
options.T_otm = options.T_fc; %Operating temperature for OTM
options.T_oxygen_pump = 323*ones(n1,n2); %Inlet temperature of vacuum pump
options.T_motor = 77*ones(n1,n2); %temperture of H2 gas after cooling superconducting motors
options.spu = 0.2*ones(n1,n2); 
options.steamratio = 0.05*ones(n1,n2); %Percentage of humidification at fuel inlet
% options.velocity_to = (ones(n1,1)*linspace(0,70,n2)); %Velocity, m/s 
options.P_non_perm = ones(n1,1)*linspace(1350,3500,n2); %Range of intake pressures for OTM, kPa
options.C1_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 1
options.airflow = ones(n1,n2); %Initial airflow, kmol/s
P_permMaxNom = 200;
options.P_perm = P_permMaxNom*ones(n1,n2); %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
options.OTM_perc_theoretical = 0.8*ones(n1,n2); %Actual percentage of theoretical O2 recovered
options.T1_eff = 0.88*ones(n1,n2); %Mechanical efficiency of turbine
options.C2_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 2
options.Blower_eff = 0.5*ones(n1,n2); %efficiency of blower
options.Blower_dP = 20*ones(n1,n2); %Pressure rise in blower in kPa
options.T0= 298*ones(n1,n2); %Ambient Temperature Kelvin;
GibbsRP.H2 = -48786*2; %kJ/kg at 2kg/kmol;
GibbsRP.H2O = -4544.6*18; %kJ/kg at 18 kg/kmol
GibbsRP.O2 = -6194*32; % kJ/kg at 32 kg/kmol
GibbsRP.rxn = -GibbsRP.H2O + GibbsRP.H2 + 0.5*GibbsRP.O2; 
options.P0 = 100*ones(n1,n2); %Ambient Pressure as a function of altitude
% options.O2ref = 0.2062*ones(n1,n2); % molar fraction of oxygen in reference atmosphere
% options.O2Xref = 3.9*ones(n1,n2); %kJ/mol
% options.N2ref = 0.7651*ones(n1,n2);
% options.H2Oref = 1.3*ones(n1,n2);% kJ/mol
% options.H2Oref = 0.0190*ones(n1,n2);
% options.rel_hum_ref = 0.6*ones(n1,n2);
mission = [];
tic
[Design,Demand,Weight,param] = run_cycle(options,mission);
toc
% x = param.plot_P;
% y = param.plot_T;
% y2 = param.plot_V;
% yyaxis left
% plot(x,y)
% yyaxis right
% plot(x,y2)
%Find Power Plant Conditions to Meet demand for climb and cruise

OD.altdata(1,:) = [param.altitude(1,1),param.velocityclimb(1,1),param.Palt(1,1),param.Talt(1,1),param.dalt(1,1),param.PrC(1,1),param.Malt(1,1)];
OD.altdata(2,:) = [param.altitude(1,5),param.velocityclimb(1,5),param.Palt(1,5),param.Talt(1,5),param.dalt(1,5),param.PrC(1,5),param.Malt(1,5)];
OD.altdata(3,:) = [param.altitude(1,15),param.velocityclimb(1,15),param.Palt(1,15),param.Talt(1,15),param.dalt(1,15),param.PrC(1,15),param.Malt(1,15)];
OD.altdata(4,:) = [param.altitude(1,45),param.velocityclimb(1,45),param.Palt(1,45),param.Talt(1,45),param.dalt(1,45),param.PrC(1,45),param.Malt(1,45)];
OD.altdata(5,:) = [param.altitude(1,95),param.velocityclimb(1,95),param.Palt(1,95),param.Talt(1,95),param.dalt(1,95),param.PrC(1,95),param.Malt(1,95)];
OD.altdata(6,:) = [param.altitude(1,195),param.velocityclimb(1,195),param.Palt(1,15),param.Talt(1,195),param.dalt(1,195),param.PrC(1,15),param.Malt(1,195)];
OD.altdata(7,:) = [param.altitude(1,295),param.velocityclimb(1,295),param.Palt(1,295),param.Talt(1,295),param.dalt(1,295),param.PrC(1,295),param.Malt(1,295)];
OD.altdata(8,:) = [param.altitude(1,395),param.velocityclimb(1,395),param.Palt(1,395),param.Talt(1,395),param.dalt(1,395),param.PrC(1,395),param.Malt(1,395)];
OD.altdata(9,:) = [param.altitude(1,450),param.velocityclimb(1,450),param.Palt(1,450),param.Talt(1,450),param.dalt(1,450),param.PrC(1,450),param.Malt(1,450)];
OD.altdata(10,:) =[param.altitude(1,450),param.Vcruise,param.PCruise,param.Tcruise,param.dcruise,Demand.PrCruise,0.82];
OD.Performance = zeros(11,15);
for z = 10:10
m1 = 10;
m2 = 10; 
OD.SOFC_area = Design.SOFC_size*ones(m1,m2); 
Pair_int_max = Design.P_airin;
OD.air_inMax = Design.air_in(1,1); %Max molar flow into compressor

OD.vflowairMax = 1.5*OD.air_inMax*28.84/1.22; % max volume flow of compressor
OD.vflowairMin = 0.75*OD.air_inMax*28.84/1.22; %Actual volume flow at sea level assuming 50% power at sea level
OD.molflowairMax = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMax/28.84; %Mol flow at max compressor volume flow with stagnation effects
OD.molflowairMin = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMin/28.84; 
OD.air_in = ones(m2,1)*linspace(OD.molflowairMin,OD.molflowairMax,m1); 
OD.Pratio = Design.P_airin/100; %Max pressure ratio
OD.height = ones(m1,1)*param.altitude; 
OD.Palt = OD.altdata(z,3)*((1 + 0.5*(1.4-1)*OD.altdata(z,7)^2)^(1.4/0.4))*ones(m1,m1);%Stagnation pressure at compressor inlet
OD.Talt = OD.altdata(z,4)*(1 + 0.5*0.4*OD.altdata(z,7)^2)*ones(m1,m1); %stagnation temperature at compressor inlet
OD.P_non_perm = OD.Pratio.*OD.Palt; 
OD.velocity = OD.altdata(z,2); %Velocity at altitude
OD.dT_fc = 50*ones(m1,m1); %Maximum temperature differential, Kelvin
OD.asr = 0.15*ones(m1,m1); % Area specific resistance, ohm-cm^2
OD.P_fc = 1000*ones(m1,m2); %Operating pressure for SOFC
OD.T_fc = 1023*ones(m1,m2); %Inlet temperature for SOFC
OD.T_otm = 1043*ones(m1,m2); %Operating temperature for OTM
OD.T_oxygen_pump = 323*ones(m1,m2); %Inlet temperature of vacuum pump
OD.T_motor = 77*ones(m1,m2); %temperture of H2 gas after cooling superconducting motors
OD.spu = 0.2*ones(m1,m1); 
OD.steamratio = 0.05*ones(m1,m1); %Percentage of humidification at fuel inlet
OD.C1_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 1
OD.P_permMin = 35;
%if 0.18*OD.Pratio*OD.alt1(z,3) < 240
OD.P_permMax = 0.18*OD.Pratio*OD.altdata(z,3); 
%else OD.P_permMax = 240;
%end
OD.P_perm = (ones(m2,1)*linspace(OD.P_permMin,OD.P_permMax,10))'; %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
OD.OTM_perc_theoretical = 0.8*ones(m1,m2); %Actual percentage of theoretical O2 recovered
OD.T1_eff = 0.88*ones(m1,m2); %Mechanical efficiency of turbine
OD.C2_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 2
OD.Blower_eff = 0.5*ones(m1,m2); %efficiency of blower
OD.Blower_dP = 20*ones(m1,m2); %Pressure rise in blower in kPa
mission = [];
tic
param_od = run_cycle_od(OD,options);
toc

[alt1,OD] = condition(OD,param_od,z);
% OD.Performance(1,1:8) = ['Altitude','Power','Pnonperm','FTE','air in','H2 used','Qbalance','H2 heat'];
% OD.Performance(1,1) = ['Power'];
% OD.Performance(1,2) = 'Non Permeate Pressure';
% OD.Performance(1,3) = 'System Efficiency'; 
% OD.Performance(1,4) = 'Air in';
% OD.Performance(1,5) = 'H2 Used';
% OD.Performance(1,6) = 'Heat Rejected';
% OD.Performance(1,7) = 'Heat Added';
% OD.Performance(1,8) = 'FC Efficiency';
% OD.Performance(1,9) = 'FC Voltage';
% OD.Performance(1,10) = 'Current Density';
% OD.Performance(1,11) = 'Turbine Work';
% OD.Performance(1,12) = 'Compressor Work'; 
% OD.Performance(1,13) = 'Recirculation Work';
% OD.Performance(1,14) = 'O2 Used';
% OD.Performance(1,15) = 'FC Power'; 
OD.Performance(z,1:15) = alt1(z,1:15); %Altitude, Power, P perm, FTE, air in, H2 Used, Qbalance, H2 used for heating
end
OD.specs = {'Net Power',OD.alt1(10,1);'System Efficiency',OD.Performance(10,3);'FC Efficiency',OD.Performance(10,8);'Current Density',OD.Performance(10,10); 'Voltage',OD.Performance(10,9);'FC Power',OD.alt1(10,15);'Turbine Work',OD.Performance(10,11);'Compressor Work',OD.alt1(10,12);'Recirculation Work',OD.Performance(10,13); 'O2 Used',OD.Performance(10,15) ;'H2 used',OD.Performance(10,5);'Air in',OD.Performance(10,4);'Heat Rejected',OD.Performance(10,6);};
%Find total fuel burned from takeoff and climb
OD.FB_to = Design.H2usedactual*param.time_to; %Fuel burned during takeoff, kg
%Fuel burned in each altitude band
OD.FB_b1 = 0.5*(OD.Performance(1,5) + OD.Performance(2,5))*param.t1; 
OD.FB_b2 = 0.5*(OD.Performance(2,5) + OD.Performance(3,5))*param.t2; 
OD.FB_b3 = 0.5*(OD.Performance(3,5) + OD.Performance(4,5))*param.t3; 
OD.FB_b4 = 0.5*(OD.Performance(4,5) + OD.Performance(5,5))*param.t4; 
OD.FB_b5 = 0.5*(OD.Performance(5,5) + OD.Performance(6,5))*param.t5; 
OD.FB_b6 = 0.5*(OD.Performance(6,5) + OD.Performance(7,5))*param.t6; 
OD.FB_b7 = 0.5*(OD.Performance(7,5) + OD.Performance(8,5))*param.t7;
OD.FB_b8 = 0.5*(OD.Performance(8,5) + OD.Performance(9,5))*param.t8; 
OD.FB_climb = OD.FB_b1 + OD.FB_b2 +OD.FB_b3 + OD.FB_b4 +OD.FB_b5+OD.FB_b6+OD.FB_b7 +OD.FB_b8;
OD.FuelReserve = OD.FB_to + OD.FB_climb + OD.Performance(11,6)*3600; %Reserve fuel sufficient to takeoff, climb and cruise for one hour; 
%Maximum Cruise distance
OD.Cruise_fuel = Design.fuelstorageactual - OD.FB_to - OD.FB_climb - OD.FuelReserve; %Total fuel available after takeoff and climb
OD.Cruise_time = (OD.Cruise_fuel)/OD.Performance(10,5); %Cruise time with useable fuel mass
OD.Cruise_distance_m = OD.Cruise_time*OD.altdata(10,2); %Total Cruise distance, meters
OD.Cruise_distance_nm = OD.Cruise_distance_m/1852; %Total Cruise distance, nautical miles 
OD.Cruise_nm_matchstandard = 7260; %standard range of 747
OD.Cruise_matchtime = 7260*1852/OD.altdata(10,2); 
OD.Cruise_matchfuelmass = OD.Cruise_matchtime*OD.Performance(10,5); %Fuel required to match standard 747 cruise range
OD.PayloadatRange = (101000 + 112760) - OD.Cruise_matchfuelmass - Design.Penalty - OD.FuelReserve; %Available payload after mass penalties 
OD.tsfc_to =Design.H2usedactual*1000000/param.TrTO; %g/kN*s
OD.tsfc_cruise = OD.Performance(10,5)*1000000/param.Thr_cruise; %g/kN*s

W0 = 396000*9.81; %Gross weight of craft
W1 = W0 - 9.81*OD.Cruise_fuel; %Empty weight of craft
TSFCcruise = OD.Performance(10,5)*1000000/Cruise.TrCruise; %g/kNs
TSFCtakeoff = Design.H2usedactual*1000000/TrTO; %g/kNs
ct = 9.81*TSFCcruise/1000000; %N/N*s
Range = 2*sqrt(2/(Cruise.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))*(sqrt(W0) - sqrt(W1)); %Total Range according to eqn 6.77, introduction to flight, at standard payload
Rangenm = Range/1852;
W1b = (sqrt(W0) - (7260*1852)/(2*sqrt(2/(Cruise.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))))^2; %Total fuel consumption to match standard 747 range
matchfuelmass = (W0 - W1b)/9.81;
payloadboost = OD.Cruise_fuel - matchfuelmass;
payloadatrange = 112760 + payloadboost; 
% m1 = 10;
% m2 = 10; 
% Cruise.SOFC_area = 0.33*param.FC_area_fixed*ones(m1,m2); 
% Pair_int_max = param.P_int_max;
% Cruise.air_inMax = param.air_in(1,1); %Max molar flow into compressor
% 
% Cruise.vflowairMax = Cruise.air_inMax*28.84/1.22; % max volume flow at sea level
% Cruise.vflowairMin = 0.2*Cruise.vflowairMax; %Min volume flow
% Cruise.molflowairMax = param.dcruise*Cruise.vflowairMax/28.84;
% Cruise.molflowairMin = param.dcruise*Cruise.vflowairMin/28.84; 
% Cruise.air_in = ones(m2,1)*linspace(Cruise.molflowairMax,Cruise.molflowairMin,m1); 
% Cruise.Pratio = Pair_int_max./100; %Max pressure ratio
% Cruise.OTM_area = 2.5e3*ones(m1,m2); %membrane area in m^2
% Cruise.height = ones(m1,1)*param.altitude; 
% Cruise.Palt = param.Pcruise*ones(m1,m1);%ones(m1,1)*param.Palt; %Pressure as a function of altitude
% Cruise.Talt = param.Tcruise*ones(m1,m1); %ones(m1,1)*param.Talt; %Temperature as a function of altitude
% Cruise.P_non_perm = Cruise.Pratio.*Cruise.Palt; 
% Cruise.velocity = param.Vcruise; %Velocity at altitude
% Cruise.dT_fc = 50*ones(m1,m1); %Maximum temperature differential, Kelvin
% Cruise.asr = 0.2*ones(m1,m1); % Area specific resistance, ohm-cm^2
% Cruise.P_fc = 1000*ones(m1,m2); %Operating pressure for SOFC
% Cruise.T_fc = 1023*ones(m1,m2); %Inlet temperature for SOFC
% Cruise.T_otm = 1043*ones(m1,m2); %Operating temperature for OTM
% Cruise.T_oxygen_pump = 323*ones(m1,m2); %Inlet temperature of vacuum pump
% Cruise.T_motor = 77*ones(m1,m2); %temperture of H2 gas after cooling superconducting motors
% Cruise.spu = 0.2*ones(m1,m1); 
% Cruise.steamratio = 0.05*ones(m1,m1); %Percentage of humidification at fuel inlet
% Cruise.C1_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 1
% Cruise.P_permMin = 50;
% Cruise.P_permMax = 175; 
% Cruise.P_perm = (ones(m2,1)*linspace(75,175,10))'; %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
% Cruise.OTM_perc_theoretical = 0.7*ones(m1,m2); %Actual percentage of theoretical O2 recovered
% Cruise.T1_eff = 0.88*ones(m1,m2); %Mechanical efficiency of turbine
% Cruise.C2_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 2
% Cruise.Blower_eff = 0.5*ones(m1,m2); %efficiency of blower
% Cruise.Blower_dP = 20*ones(m1,m2); %Pressure rise in blower in kPa
% Cruise.Preq = param.PrCruise;
% mission = [];
% tic
% param_od = run_cycle_od(Cruise,options);
% toc
% err_cruise = param_od.NetPower - Cruise.Preq*ones(m1,m1); 

% H2_heat = zeros(m1,m1);
% param_od.H2_total = zeros(m1,m2);


%plot(param.plot_D,param.plot_T)
%plot(param.P_non_perm(1,1:10),param.sofc_area(1,1:10))
%plot(options.P_non_perm(1,1:10),param.Qbalance(1:10,1))
%plot_case(param,param,'thrust_err','velocity','FCPower')
%plot_case(param,param,'Qbalance','P_non_perm','i_den')
%plot(param.i_total, param.P_non_perm)
%plot_case(param,param,'Qbalance','P_non_perm','sofc_area')
%plot(param.FTE(1,1:10),param.P_den(1,1:10)); 
%plot_case(param,param,'P_den','FTE','i_den')
%plot_case(param_od,param_od,'Eout','P_perm','air_in'); 