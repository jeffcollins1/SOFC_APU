%%  Test Cycle b sized
%Create options to test
%This version Tests Range for a previously designed cycle with added
%battery assist

n1 = 10; % number of points in test dimension 1
n2 = 10; % number of points in test dimension 1
options = [];
[Demand,param.Ppropmax,param.Pprop_to,param.TrTO,param.time_to,param.velocityclimb,param.time_climb,param.Palt,param.Talt,param.dalt,param.Malt,param.PrC,param.altitude,param.dcruise,param.Pcruise,param.Tcruise,param.Vcruise,param.t1,param.t2,param.t3,param.t4,param.t5,param.t6,param.t7,param.t8,param.Thr_to,param.Thr_cruise,clcruise,cdcruise] = craft(options);
%options.height = 10*ones(n1,n2); %Altitude, meters
% options.OTM_area = 2.5e3*ones(n1,n2); %membrane area in m^2
options.SOFC_area = (ones(n1,1)*linspace(1000,4000,n2))';%4e3*ones(n1,n2); 
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
options.P_non_perm = (param.Pcruise*((1+ 0.5*0.4*0.82^2)^(1/0.4))/100)*ones(n1,1)*linspace(1500,3000,n2); %Range of intake pressures for OTM, kPa
options.C1_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 1
options.airflow = ones(n1,n2); %Initial airflow, kmol/s
%P_permmax = 0.18*OD.Pratio*OD.altdata(z,3)
P_permMaxNom = 200;
%options.P_perm = P_permMaxNom*ones(n1,n2); %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
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
options.Pamb = param.Pcruise*ones(10,10); %Ambient Pressure as a function of altitude
options.Tamb = param.Tcruise*ones(10,10); %Ambient Temperature as a function of altitude
options.P0 = options.Pamb.*((1+ 0.5*0.4*0.82^2)^(1/0.4));
options.T0 = options.Tamb.*(1 + 0.5*0.4*0.82^2);
% options.O2ref = 0.2062*ones(n1,n2); % molar fraction of oxygen in reference atmosphere
% options.O2Xref = 3.9*ones(n1,n2); %kJ/mol
% options.N2ref = 0.7651*ones(n1,n2);
% options.H2Oref = 1.3*ones(n1,n2);% kJ/mol
% options.H2Oref = 0.0190*ones(n1,n2);
% options.rel_hum_ref = 0.6*ones(n1,n2);
mission = [];

P_permmin = options.P0;
%options.P_perm = P_permMaxNom*ones(n1,n2);
DesignOptions = zeros(16,10);
DesignOptionsTotalPower = zeros(16,10);
DesignOptionsMass = zeros(7,10); 
DesignOptionsTotalMass = zeros(7,10); 
%Determine most efficient operating conditions for 10 different OTM intake
%pressures at altitude.  
% for j = 1:10
%   
%    P_permmax = 0.19*options.P_non_perm(1,j); 
%     P_non_perm = options.P_non_perm(1,j);
% for y = 1:10
%     options.P_perm =  35 + (y/10)*(P_permmax - 35); %ones(10,1)*linspace(35,P_permmax,10); 
%     [Design,param] = run_cycle_b(options,param,Demand,P_non_perm);
%     DesignOptions(1,y) = Design.P_perm;
%     DesignOptions(2,y) = P_non_perm; 
%     DesignOptions(3,y) = Design.MaxPower; 
%     DesignOptions(4,y) = Design.FTE;
%     DesignOptions(5,y) = Design.FCeff;
%     DesignOptions(6,y) = Design.iden;
%     DesignOptions(7,y) = Design.Voltage;
%     DesignOptions(8,y) = Design.SOFC_power; 
%     DesignOptions(9,y) = Design.T1_work;
%     DesignOptions(10,y) = Design.C1_work;
%     DesignOptions(11,y) = Design.BlowerWork;
%     DesignOptions(12,y) = Design.H2_usednom;
%     DesignOptions(13,y) = Design.Qbal; 
%     DesignOptions(14,y) = Design.SOFC_size;%*10000/81; 
%     DesignOptions(15,y) = P_non_perm/options.P0(1,1); 
%     DesignOptions(16,y) = Design.air_in;
% %Design.H2_usednom = scale*param.H2_used(x,y)*2; %Fuel consumption, kg/s
% % DesignOptionsMass(1,y) = Design.mass;
% % DesignOptionsMass(2,y) = Design.Comp_mass; 
% % DesignOptionsMass(3,y) = Design.Turb_mass; 
% % % Design.SOFC_size = scale*param.sofc_area(x,y); 
% % DesignOptionsMass(4,y)=  Design.SOFC_mass; 
% % DesignOptionsMass(5,y) = Design.OTM_mass; 
% % DesignOptionsMass(6,y)=  Design.HTSM_mass;
% % DesignOptionsMass(7,y) = Design.HX_mass; 
% 
% % Design.O2max = scale*param.O2_used(x,y); 
% % Design.MaxPower = scale*P_net; 
% % Design.Penalty = scale*Weight.Total(x,y) - 0.7*4*(9670/2.2); %Total Weight Penalty of option 1
% % Design.P_den = param.P_den(x,y); 
% % Design.FTE_to = param.FTE(x,y);
% % Design.Qbal = scale*param.Qbalance(x,y);
% % Design.FCQgen = scale*param.FCQgen(x,y);
% % Design.FCQremove = scale*param.FCQout(x,y);
% % Design.OTMQin =  scale*param.OTMheat_in(x,y);
% % Design.OTMQout = scale*param.OTMheat_out(x,y);
% % Design.Qfuelout =  scale*param.Qremovefuel(x,y);
% % Design.HXmass = scale*Weight.hx(x,y); 
% % Design.FCeff = param.FC_eff(x,y); 
% % Design.iden = param.iden(x,y);
% % Design.Voltage = param.FCVoltage(x,y); 
% % Design.SOFC_power = scale*param.FCPower(x,y);
% % Design.T1_work = scale*param.T1_work(x,y);
% % Design.C1_work = scale*param.C1_work(x,y); 
% % Design.BlowerWork = scale*param.BlowerWork(x,y); 
% % Design.scale = scale; 
% % tic
% % [Design,Weight,param] = run_cycle_b(options,param,Demand);
% % toc
% 
% end
% maxeff = max(DesignOptions(4,1:10));
% [x,v] = find(DesignOptions(4,1:10)==maxeff);
% DesignOptionsTotalPower(1:16,j) = DesignOptions(1:16,v); 
% DesignOptionsTotalMass(1:7,j) = DesignOptionsMass(1:7,v); 
% end
[Demand,param.Ppropmax,param.Pprop_to,param.TrTO,param.time_to,param.velocityclimb,param.time_climb,param.Palt,param.Talt,param.dalt,param.Malt,param.PrC,param.altitude,param.dcruise,param.Pcruise,param.Tcruise,param.Vcruise,param.t1,param.t2,param.t3,param.t4,param.t5,param.t6,param.t7,param.t8,param.Thr_to,param.Thr_cruise,clcruise,cdcruise] = craft(options);
[Battery] = battery(options);
% x = param.plot_P;
% y = param.plot_T;
% y2 = param.plot_V;
% yyaxis left
% plot(x,y)
% yyaxis right
% plot(x,y2)

 %%sized for maximum cruise efficiency 
  

% DesignOptionsTotalPower = [35	35	35	35	35	35	35	35	35	35; ...
% 315.2475203	350.2750226	385.3025249	420.3300271	455.3575294	490.3850317	525.4125339	560.4400362	595.4675384	630.4950407;...
% 62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849;...
% 0.675958646	0.681310144	0.684990877	0.687682998	0.689745589	0.691384188	0.692724183	0.693846019	0.69480369	0.695634495;...
% 0.703182597	0.695424197	0.689191655	0.684078447	0.67980908	0.676191279	0.673086817	0.670393697	0.668035664	0.665953789;...
% 0.859392658	0.944003153	1.011692877	1.067076257	1.113229678	1.152283001	1.185757591	1.214769136	1.240154416	1.26255333;...
% 0.904191003	0.894214825	0.886200679	0.879625834	0.874136045	0.869484077	0.865492188	0.862029227	0.858997138	0.856320148;...
% 56029.68314	56711.56782	57176.69533	57511.61887	57762.40293	57955.81371	58108.45532	58231.16389	58331.30774	58414.06173;...
% 25703.91605	23713.73356	22331.58417	21315.22163	20536.0411	19919.42379	19419.10115	19004.86127	18656.1228	18358.39902;...
% -18717.23899	-17439.40981	-16554.41049	-15905.1438	-15408.37669	-15015.9101	-14697.91409	-14434.94623	-14213.78271	-14025.13179;...
% -200.8922599	-205.6056324	-209.16653	-211.9643562	-214.2256368	-216.0929479	-217.6613898	-218.9972692	-220.148239	-221.1497548;...
% 0.642239035	0.65730737	0.66869132	0.677635782	0.684864944	0.690834612	0.695848815	0.700119531	0.703799104	0.707000883;...
% 0	0	0	0	0	0	0	0	0	0; ...
% 7210.516017	6718.252826	6377.321037	6127.201462	5935.829902	5784.638446	5662.135585	5560.831438	5475.631741	5402.956999;...
% 15	16.66666667	18.33333333	20	21.66666667	23.33333333	25	26.66666667	28.33333333	30; ...
% 1.802629004	1.679563207	1.594330259	1.531800366	1.483957476	1.446159611	1.415533896	1.390207859	1.368907935	1.35073925;]

CycleE = [37.48970289; 315.2475203; 62368.94849; 0.631368939; 0.63487119; 0.797907612; 0.816352425; 62001.8732; 38789.84971; -34672.85361; -246.2251115; 0.787165111; 0; 9518.631518; 15; 2.379657879]; %Cycle designation from previously generated designs

% maxpayload = max(Final.Payload); 
% [l,m] = find(Final.Payload==maxpayload);
OD.altdata(1,:) = [param.altitude(1,1),param.velocityclimb(1,1),param.Palt(1,1),param.Talt(1,1),param.dalt(1,1),Demand.PrCruise,param.Malt(1,1)];
OD.altdata(2,:) = [param.altitude(1,5),param.velocityclimb(1,5),param.Palt(1,5),param.Talt(1,5),param.dalt(1,5),Demand.PrCruise,param.Malt(1,5)];
OD.altdata(3,:) = [param.altitude(1,15),param.velocityclimb(1,15),param.Palt(1,15),param.Talt(1,15),param.dalt(1,15),Demand.PrCruise,param.Malt(1,15)];
OD.altdata(4,:) = [param.altitude(1,45),param.velocityclimb(1,45),param.Palt(1,45),param.Talt(1,45),param.dalt(1,45),Demand.PrCruise,param.Malt(1,45)];
OD.altdata(5,:) = [param.altitude(1,95),param.velocityclimb(1,95),param.Palt(1,95),param.Talt(1,95),param.dalt(1,95),Demand.PrCruise,param.Malt(1,95)];
OD.altdata(6,:) = [param.altitude(1,195),param.velocityclimb(1,195),param.Palt(1,15),param.Talt(1,195),param.dalt(1,195),Demand.PrCruise,param.Malt(1,195)];
OD.altdata(7,:) = [param.altitude(1,295),param.velocityclimb(1,295),param.Palt(1,295),param.Talt(1,295),param.dalt(1,295),Demand.PrCruise,param.Malt(1,295)];
OD.altdata(8,:) = [param.altitude(1,395),param.velocityclimb(1,395),param.Palt(1,395),param.Talt(1,395),param.dalt(1,395),Demand.PrCruise,param.Malt(1,395)];
OD.altdata(9,:) = [param.altitude(1,450),param.velocityclimb(1,450),param.Palt(1,450),param.Talt(1,450),param.dalt(1,450),Demand.PrCruise,param.Malt(1,450)];
OD.altdata(10,:) =[param.altitude(1,450),param.Vcruise,param.Pcruise,param.Tcruise,param.dcruise,Demand.PrCruise,0.82];
OD.altdata(11,:) =[param.altitude(1,450),param.Vcruise,param.Pcruise,param.Tcruise,param.dcruise,Demand.PrCruise2,0.82];
OD.Performance = zeros(11,18);
OD.Mass = zeros(10,6);
Final.TurbineMass = zeros(11,1);
Final.CompMass = zeros(11,1); 
Final.HXMass = zeros(11,1); 
Final.OTMMass = zeros(11,1);
Final.Mass = zeros(11,1);
Final.SOFCMass =zeros(11,1);
Final.MassPenalty = zeros(11,1);
Final.FuelPenalty = zeros(11,1);
Final.Payload = zeros(11,1);
Final.TSFC_to = zeros(11,1);
Final.TSFC_cruise = zeros(11,10); 
Performance = zeros(11,18,10); 
PerformanceTable = zeros(110,18); 
%For 10 unique SOFC sizes, OD performance can be examined for each.
%After allSOFC sizes have been examined across all power demands and
%ambient conditions, total range and payload are stored for each.  Pressure
%ratio into the OTM is constant from Design conditions, intake flow rate and OTM permeate
%pressure are varied to meet changing demands.  Turbomachinery and OTM masses are also calculated so that the final system mass, payload and range can be determined  
DemandRange = [Demand.Pr_to,Demand.PrCb1,Demand.PrCb2,Demand.PrCb3,Demand.PrCb4,Demand.PrCb5,Demand.PrCb6,Demand.PrCb7,Demand.PrCb8,Demand.PrCruise,Demand.PrCruise2];
for b = 1:10 %b index selects scalar of SOFC size for range of current densities 
    
for z = 1:11 %z index selects ambient conditions and flight velocity at a particular altitude band
m1 = 10;
m2 = 10; 
OD.SOFC_area = CycleE(14,1)*(0.35 + b/20)*ones(10,10); %Fuel cell active area based index b
    Pair_int_max = CycleE(15,1)*(OD.altdata(z,3)*(1 + 0.5*(1.4-1)*OD.altdata(z,7)^2)^(1.4/0.4)); %OTM intake pressure based pressure ratio from system b and stagnation pressure at condition z 
    OD.air_inmax = CycleE(16,1); %Flow rate of ambient air into system based on airflow of system b
    P_perm_max = 0.2*Pair_int_max; %Set upper limit for varying permeate pressure
    P_perm_min = 35; %Set lower limit for permeate pressure based on outlet temperature of oxygen compressor
    OD.P_perm = (ones(m2,1)*linspace(P_perm_min,P_perm_max,10))'; %Create range from min to max permeate pressure
    %OD.P_perm = DesignOptionsTotalPower(1,b); 
  
% Pair_int_max = Design.P_airin;
% OD.air_inMax = Design.air_in(1,1); %Max molar flow into compressor
OD.Pinmotor = DemandRange(1,z); 
OD.vflowairMax = OD.air_inmax*(28.84/OD.altdata(10,5)); % max volume flow of compressor
OD.vflowairMin = 0.5*OD.air_inmax*28.84/(OD.altdata(10,5)); %Actual volume flow at sea level assuming 50% power at off design condition
OD.molflowairMax = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMax/28.84; %Mol flow at max compressor volume flow with stagnation effects
OD.molflowairMin = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMin/28.84; %Mol flow at 
OD.air_in = ones(m2,1)*linspace(OD.molflowairMin,OD.molflowairMax,m1); 
OD.Pratio = CycleE(15,1); %Max pressure ratio
OD.Pratiot = Pair_int_max/(OD.altdata(z,3)); 
OD.height = ones(m1,1)*param.altitude; 
OD.Palt = OD.altdata(z,3)*((1 + 0.5*(1.4-1)*OD.altdata(z,7)^2)^(1.4/0.4))*ones(m1,m1);%Stagnation pressure at compressor inlet
OD.Talt = OD.altdata(z,4)*(1 + 0.5*0.4*OD.altdata(z,7)^2)*ones(m1,m1); %stagnation temperature at compressor inlet
OD.P0 = OD.altdata(z,3)*ones(10,10);
OD.T0 = OD.altdata(z,4)*ones(10,10);
OD.P_non_perm = OD.Pratio.*OD.Palt; 
OD.velocity = OD.altdata(z,2); %Velocity at altitude
OD.dT_fc = 50*ones(m1,m1); %Maximum temperature differential, Kelvin
OD.asr = 0.3*ones(m1,m1); % Area specific resistance, ohm-cm^2
OD.P_fc = 1000*ones(m1,m2); %Operating pressure for SOFC
OD.T_fc = 1023*ones(m1,m2); %Inlet temperature for SOFC
OD.T_otm = 1043*ones(m1,m2); %Operating temperature for OTM
OD.T_oxygen_pump = 323*ones(m1,m2); %Inlet temperature of vacuum pump
OD.T_motor = 77*ones(m1,m2); %temperture of H2 gas after cooling superconducting motors
OD.spu = 0.2*ones(m1,m1); 
OD.steamratio = 0.05*ones(m1,m1); %Percentage of humidification at fuel inlet
OD.C1_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 1
%OD.P_permMin = 35;
%if 0.18*OD.Pratio*OD.alt1(z,3) < 240
%OD.P_permMax = 0.18*OD.Pratio*OD.altdata(z,3); 
%else OD.P_permMax = 240;
%end
%OD.P_perm = (ones(m2,1)*linspace(OD.P_permMin,OD.P_permMax,10))'; %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
OD.OTM_perc_theoretical = 0.8*ones(m1,m2); %Actual percentage of theoretical O2 recovered
OD.T1_eff = 0.88*ones(m1,m2); %Mechanical efficiency of turbine
OD.C2_eff = 0.80*ones(m1,m2); %Mechanical efficiency of compressor 2
OD.Blower_eff = 0.5*ones(m1,m2); %efficiency of blower
OD.Blower_dP = 20*ones(m1,m2); %Pressure rise in blower in kPa
mission = [];
tic
[param_od] = run_cycle_od_b_battery(OD,options,z);
toc

[alt1,OD] = condition(OD,param_od,z);
% OD.Performance(1,1:8) = ['Altitude','Power','Pnonperm','FTE','air in','H2 used','Qbalance','H2 heat'];
% OD.Performance(1,1) = ['Power'];
% OD.Performance(1,2) = 'Permeate Pressure';
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
%OD.Performance(1,16) = Q hx total
%OD.Performance(1,17) = mass flow into turbine;
OD.Performance(z,1:18) = alt1(z,1:18); %Altitude, Power, P perm, FTE, air in, H2 Used, Qbalance, H2 used for heating
OD.CompIn = OD.Performance(z,4); %molar flow of air into comp
OD.HXLoad = OD.Performance(z,16); %Total heat exchanged
OD.OTMflux = OD.Performance(z,14); %Total OTM Flux
OD.TurbIn = OD.Performance(z,17); %Total molar flow into turbine
OD.H2used = OD.Performance(z,5); % Total fuel flow 
 
if OD.Performance(z,1) == NaN
    OD.Mass(z,1:6) = NaN;
else
     [Weight] = weight_b_sized(OD); 
OD.Mass(z,1:6) = [Weight.Total,Weight.turb,Weight.comp,Weight.hx,Weight.otm,Weight.sofc(1,1)]; %Total mass for a single power plant across a variety of conditions
end
l = 11*b;
m = l-10; 
PerformanceTable(m:l,1:18) = OD.Performance; 
Performance(:,:,b) = OD.Performance; 
end
Final.TurbineMass(b,1) = max(OD.Mass(1:10,2));
Final.CompMass(b,1) = max(OD.Mass(1:10,3)); 
Final.HXMass(b,1) = max(OD.Mass(1:10,4));
Final.OTMMass(b,1) = max(OD.Mass(1:10,5));
Final.SOFCMass(b,1) = max(OD.Mass(1:10,6)); 
Final.Mass(b,1) = (Final.TurbineMass(b,1) + Final.CompMass(b,1) + Final.HXMass(b,1) + Final.OTMMass(b,1) + Final.SOFCMass(b,1))*1.1; %Find total mass required for power plant index b to perform across all power demands and ambient conditions with 10% for hardware
%OD.specs = {'Net Power',OD.alt1(10,1);'System Efficiency',OD.Performance(10,3);'FC Efficiency',OD.Performance(10,8);'Current Density',OD.Performance(10,10); 'Voltage',OD.Performance(10,9);'FC Power',OD.alt1(10,15);'Turbine Work',OD.Performance(10,11);'Compressor Work',OD.alt1(10,12);'Recirculation Work',OD.Performance(10,13); 'O2 Used',OD.Performance(10,15) ;'H2 used',OD.Performance(10,5);'Air in',OD.Performance(10,4);'Heat Rejected',OD.Performance(10,6);};
%Find total fuel burned from takeoff and climb
OD.FB_to = OD.Performance(1,5)*param.time_to; %Fuel burned during takeoff, kg
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
OD.FuelReserve = OD.FB_to + OD.FB_climb + OD.Performance(10,5)*3600; %Reserve fuel sufficient to takeoff, climb and cruise for one hour; 
%Maximum Cruise distance
W0 = 396000*9.81; %Gross weight of craft
%W1 = W0 - 9.81*OD.Cruise_fuel; %Empty weight of craft
OD.Cruise_matchtime = 7260*1852/OD.altdata(10,2); %total time to travel desired cruise range at cruise velocity
TSFCcruise1 = OD.Performance(10,5)*1000000/param.Thr_cruise; %g/kNs
TSFCcruise2 = OD.Performance(11,5)*1000000/param.Thr_cruise;
TSFCcruise = ((Battery.Total/1000)*TSFCcruise2 + TSFCcruise1*(OD.Cruise_matchtime - Battery.Total/1000))/OD.Cruise_matchtime; %True TSFC at cruise based on weighted average of time spent recharging the battery at 1 MJ/s
TSFCtakeoff = OD.Performance(1,5)*1000000/param.Thr_to; %g/kNs
ct = 9.81*TSFCcruise/1000000; %N/N*s
W1b = (sqrt(W0) - (7260*1852)/(2*sqrt(2/(param.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))))^2; %Total craft weight after burning cruise fuel, with 1 nautical mile = 1852 meters
Range = 2*sqrt(2/(param.dcruise*511))*(sqrt(clcruise)/(ct*cdcruise))*(sqrt(W0) - sqrt(W1b)); %Total Range according to eqn 6.77, introduction to flight, at standard payload
Rangenm = Range/1852;

matchfuelmass = (W0 - W1b)/9.81; %Total fuel required for cruise distance
matchfuelstoragemass = 1.15*matchfuelmass; 
Final.MassPenalty(b,1) = Final.Mass(b,1) - 0.7*4*(9670/2.2); %Total Penalty of Power Plant
Final.FuelPenalty(b,1) = 101000 - Final.MassPenalty(b,1) - matchfuelstoragemass - OD.FuelReserve; %Available fuel storage without cutting into payload
Final.Payload(b,1) = 112760 + Final.FuelPenalty(b,1); %Total Payload Reduction at standard range
Final.TSFC_to(b,1) = TSFCtakeoff;
Final.TSFC_cruise(b,1) = TSFCcruise; 
Final.H2storedtotal(b,1) = matchfuelmass + OD.FuelReserve;
Final.H2reserve(b,1) = OD.FuelReserve;

% OD.Cruise_fuel = W1b - OD.FB_to - OD.FB_climb - OD.FuelReserve; %Total fuel available after takeoff and climb
% OD.Cruise_time = (OD.Cruise_fuel)/OD.Performance(10,5); %Cruise time with useable fuel mass
% OD.Cruise_distance_m = OD.Cruise_time*OD.altdata(10,2); %Total Cruise distance, meters
% OD.Cruise_distance_nm = OD.Cruise_distance_m/1852; %Total Cruise distance, nautical miles 
% OD.Cruise_nm_matchstandard = 7260; %standard range of 747
% OD.Cruise_matchfuelmass = OD.Cruise_matchtime*OD.Performance(10,5); %Fuel required to match standard 747 cruise range
% OD.PayloadatRange = (101000 + 112760) - OD.Cruise_matchfuelmass - Design.Penalty - OD.FuelReserve; %Available payload after mass penalties 
% OD.tsfc_to =Design.H2usedactual*1000000/param.TrTO; %g/kN*s
% OD.tsfc_cruise = OD.Performance(10,5)*1000000/param.Thr_cruise; %g/kN*s

end
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