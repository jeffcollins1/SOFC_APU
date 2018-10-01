%%  Test Cycle b with Battery Assist
%Create options to test
%This version sizes for max cruise efficiency, 'Design', and examines
%takeoff performance in 'OD'.  The design portion selects an OTM intake
%pressure at cruise, then compares system performance with 10 OTM oxygen
%pressures and SOFC stack sizes to give 100 unique operating points and
%system metrics.  The system with highest efficiency is selected, values
%are saved, then repeats the process for another intake pressure.  All
%system masses are saved to 'DesignOptionsTotalMass', and performance metrics
%are saved to 'DesignOptionsTotalPower'


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
options.T_otm = options.T_fc + 0.25*options.dT_fc; %Operating temperature for OTM
options.T_oxygen_pump = 323*ones(n1,n2); %Inlet temperature of vacuum pump
options.T_motor = 77*ones(n1,n2); %temperture of H2 gas after cooling superconducting motors
options.spu = 0.2*ones(n1,n2); 
options.steamratio = 0.05*ones(n1,n2); %Percentage of humidification at fuel inlet
% options.velocity_to = (ones(n1,1)*linspace(0,70,n2)); %Velocity, m/s 
options.P_non_perm = (param.Pcruise*((1+ 0.5*0.4*0.82^2)^(1/0.4))/100)*ones(n1,1)*linspace(1350,3000,n2); %Range of intake pressures for OTM, kPa
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
%and permeate pressures at altitude.  
% for j = 1:10
%   
%    P_permmax = 0.19*options.P_non_perm(1,j); 
%     P_non_perm = options.P_non_perm(1,j);
% for y = 1:10
%     options.P_perm = 35 + ((y/10)*(P_permmax-35)); 
%     [Design,param] = run_cycle_b(options,param,Demand,P_non_perm);
%     DesignOptions(1,y) = options.P_perm;
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


% x = param.plot_P;
% y = param.plot_T;
% y2 = param.plot_V;
% yyaxis left
% plot(x,y)
% yyaxis right
% plot(x,y2)

 %%sized for maximum cruise efficiency 
 
% DesignOptionsTotalPower = [37.10700983	40.60110875	48.97639568	51.75057386	59.40594004	62.87366277	66.34138549	76.77092986	73.27683094	85.0934644;...
% 283.7227683	322.2530208	360.7832733	399.3135258	437.8437783	476.3740307	514.9042832	553.4345357	591.9647882	630.4950407;...
% 62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849;...
% 0.788299789	0.789005255	0.791772356	0.7912093	0.792993704	0.793976107	0.793101055	0.793135005	0.793571462	0.795195326;...
% 0.717407584	0.71392035	0.721042362	0.716541658	0.720989627	0.718282323	0.715990695	0.72322592	0.712324731	0.720464341;...
% 0.703437853	0.741923833	0.663143261	0.713007464	0.663729311	0.693758235	0.719092655	0.638836673	0.759486082	0.669564422;...
% 0.921949663	0.917468174	0.926620762	0.92083685	0.926552991	0.923073801	0.920128801	0.92942688	0.915417623	0.925877941;...
% 56515.55501	56379.01745	56797.45397	56278.06082	56601.176	56381.137	56104.52848	56581.47182	55983.54602	56330.45345;...
% 31775.92531	29243.36	35352.79136	31221.3048	35199.74243	32896.21208	30863.92515	37184.60542	28373.11406	34728.35369;...
% -22882.78445	-20320.62522	-27085.68026	-22507.96208	-26970.21827	-24510.96117	-22265.7976	-29224.93464	-19755.64476	-26626.70383;...
% -198.7311778	-199.2194397	-198.715645	-198.1352071	-198.0434172	-198.0170674	-197.6762581	-197.3623152	-198.2651367	-197.23988;...
% 0.635330201	0.636891142	0.635280544	0.633424923	0.633131477	0.633047238	0.631957694	0.630954039	0.633840299	0.630562622;...
% 0	0	50.45883742	0	0	0	0	0	40.44049832	0;...
% 8714.351073	8282.609248	9243.138934	8571.609103	9203.73689	8804.186907	8479.38666	9529.479336	8052.325371	9086.510678;...
% 13.5	15.33333333	17.16666667	19	20.83333333	22.66666667	24.5	26.33333333	28.16666667	30;...
% 2.178587768	2.070652312	2.310784734	2.142902276	2.300934223	2.201046727	2.119846665	2.382369834	2.013081343	2.27162767;]

%%sized with corrrected OTM temperature

% DesignOptionsTotalPower = [36.8907326	37.6228074	38.35488219	39.08695699	39.81903179	40.55110658	47.56636276	49.03051236	50.49466195	51.95881155;...
% 283.7227683	322.2530208	360.7832733	399.3135258	437.8437783	476.3740307	514.9042832	553.4345357	591.9647882	630.4950407;...
% 62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849;...
% 0.647680417	0.654452984	0.658441828	0.661039332	0.662860484	0.66421322	0.66528544	0.666506253	0.667550385	0.668462637;...
% 0.711168603	0.70163197	0.694273584	0.688432371	0.683685828	0.679754198	0.68466238	0.682092651	0.679872215	0.677934508;...
% 0.709584822	0.814928445	0.895573269	0.959294148	1.010912597	1.053577487	1.000302791	1.028210789	1.052297927	1.073298599;...
% 0.914973443	0.902703827	0.893236694	0.88572152	0.879614725	0.87455638	0.880871136	0.87756498	0.874708217	0.872215207;...
% 60724.54751	61727.16769	62354.89984	62782.98898	63091.66292	63323.04761	62767.32437	62895.76625	62996.19886	63075.52026;...
% 39006.24171	34571.89672	31815.27464	29934.74631	28569.17064	27531.88642	28718.89303	28009.13925	27422.50478	26929.16274;...
% -34062.68317	-30558.9877	-28387.70493	-26910.41343	-25840.044	-25028.50379	-25942.81689	-25385.59494	-24925.26577	-24538.27325;...
% -215.1597435	-221.6849887	-226.3128689	-229.7999926	-232.5330578	-234.7357383	-231.0077	-232.3524991	-233.4835878	-234.4457727;...
% 0.68785122	0.70871199	0.723507011	0.734655111	0.743392537	0.750434358	0.73851607	0.742815303	0.746431317	0.749507357;...
% 0	0	0	0	0	0	0	0	0	0; ...
% 9353.011774	8390.958818	7794.762881	7389.124701	7095.220142	6872.38552	7123.439765	6970.43645	6844.038183	6737.776865;...
% 13.5	15.33333333	17.16666667	19	20.83333333	22.66666667	24.5	26.33333333	28.16666667	30; ...
% 2.338252944	2.097739705	1.94869072	1.847281175	1.773805036	1.71809638	1.780859941	1.742609112	1.711009546	1.684444216];


% DesignOptionsTotalPower = [36.60700983	37.30055437	37.99409892	38.68764346	39.38118801	40.07473255	46.5365542	47.92364329	49.31073238	58.5467322;...
% 283.7227683	322.2530208	360.7832733	399.3135258	437.8437783	476.3740307	514.9042832	553.4345357	591.9647882	630.4950407;...
% 62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849	62368.94849;...
% 0.651131605	0.657732745	0.661614535	0.664138088	0.665904563	0.667214611	0.668399547	0.669565812	0.670564689	0.671242299;...
% 0.716074556	0.706530866	0.699170597	0.693329306	0.688583796	0.684653349	0.688879746	0.68631552	0.684100048	0.69056894;...
% 0.717632392	0.82273607	0.903199657	0.966778837	1.018283473	1.060854949	1.01507493	1.042860593	1.066842553	0.996752854;...
% 0.920768197	0.908496394	0.899032153	0.891521098	0.885419059	0.880365074	0.885799608	0.882502385	0.879653608	0.887971667;...
% 61318.45864	62248.4737	62831.77268	63230.00588	63517.33679	63732.77696	63253.40941	63369.34989	63459.62918	62829.46741;...
% 38167.80604	33840.69472	31148.78707	29311.58633	27977.0664	26963.17361	27947.5154	27264.37557	26699.43674	28262.65144;...
% -33796.07455	-30329.9948	-28180.43538	-26717.25503	-25656.77661	-24852.57565	-25619.91755	-25076.35914	-24627.07797	-25852.59578;...
% -215.8967667	-222.1317914	-226.5736015	-229.9306232	-232.5672907	-234.695767	-231.5014283	-232.7922851	-233.8789098	-229.387363;...
% 0.690207434	0.710140388	0.724340556	0.735072728	0.743501977	0.750306572	0.740094487	0.744221269	0.747695135	0.733335962;...
% 0	0	0	0	0	0	0	0	0	0;...
% 9279.805752	8328.081409	7737.85032	7336.086814	7044.898151	6824.078757	7034.777305	6885.525755	6762.161073	7098.666642;...
% 13.5	15.33333333	17.16666667	19	20.83333333	22.66666667	24.5	26.33333333	28.16666667	30;...
% 2.319951438	2.082020352	1.93446258	1.834021703	1.761224538	1.706019689	1.758694326	1.721381439	1.690540268	1.774666666]

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
OD.Performance = zeros(10,18);
OD.Mass = zeros(10,6);
Final.TurbineMass = zeros(10,1);
Final.CompMass = zeros(10,1); 
Final.HXMass = zeros(10,1); 
Final.OTMMass = zeros(10,1);
Final.Mass = zeros(10,1);
Final.SOFCMass =zeros(10,1);
Final.MassPenalty = zeros(10,1);
Final.FuelPenalty = zeros(10,1);
Final.Payload = zeros(10,1);
Final.UseableFuel = zeros(10,1);
Final.TSFC_to = zeros(10,1);
Final.TSFC_cruise = zeros(10,1); 
Performance = zeros(10,18,10); 
PerformanceTable = zeros(100,18);
Final.HTSMMass= zeros(10,1); 
Final.Batterymass = zeros(10,1); 
%For 10 unique system sizes, OD performance can be examined for each.
%After all system sizes have been examined across all power demands and
%ambient conditions, total range and payload are stored for each.  Pressure
%ratio into the OTM is constant from Design conditions, intake flow rate and OTM permeate
%pressure are varied to meet changing demands.  Turbomachinery and OTM masses are also calculated so that the final system mass, payload and range can be determined  
DemandRange = [Demand.Pr_to,Demand.PrCb1,Demand.PrCb2,Demand.PrCb3,Demand.PrCb4,Demand.PrCb5,Demand.PrCb6,Demand.PrCb7,Demand.PrCb8,Demand.PrCruise]; %,Demand.PrCruise2];
for b = 1:10 %b index selects power cycle 
    
for z = 1:10 %z index selects ambient conditions and flight velocity at a particular altitude band
m1 = 10;
m2 = 10; 
OD.SOFC_area = DesignOptionsTotalPower(14,b); %Fuel cell active area based index b
    Pair_int_max = DesignOptionsTotalPower(15,b)*(OD.altdata(z,3)*(1 + 0.5*(1.4-1)*OD.altdata(z,7)^2)^(1.4/0.4)); %OTM intake pressure based pressure ratio from system b and stagnation pressure at condition z 
    OD.air_inmax = DesignOptionsTotalPower(16,b); %Flow rate of ambient air into system based on airflow of system b
    P_perm_max = 0.19*Pair_int_max; %Set upper limit for varying permeate pressure
    P_perm_min = 35 ; %Set lower limit for permeate pressure based on outlet temperature of oxygen compressor
    OD.P_perm = (ones(m2,1)*linspace(P_perm_min,P_perm_max,10))'; %Create range from min to max permeate pressure
    %OD.P_perm = DesignOptionsTotalPower(1,b); 
  OD.Pinmotor = DemandRange(1,z);
 
% Pair_int_max = Design.P_airin;
% OD.air_inMax = Design.air_in(1,1); %Max molar flow into compressor

OD.vflowairMax = OD.air_inmax*(28.84/OD.altdata(10,5)); % max volume flow of compressor
OD.vflowairMin = 0.5*OD.air_inmax*28.84/(OD.altdata(10,5)); %Actual volume flow at sea level assuming 50% power at off design condition
OD.molflowairMax = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMax/28.84; %Mol flow at max compressor volume flow with stagnation effects
OD.molflowairMin = OD.altdata(z,5)*((1+ 0.5*0.4*OD.altdata(z,7)^2)^(1/0.4))*OD.vflowairMin/28.84; %Mol flow at 
OD.air_in = ones(m2,1)*linspace(OD.molflowairMin,OD.molflowairMax,m1); 
OD.Pratio = DesignOptionsTotalPower(15,b); %Max pressure ratio
OD.Pratiot = Pair_int_max/(OD.altdata(z,3)); 
OD.height = ones(m1,1)*param.altitude; 
OD.Palt = OD.altdata(z,3)*((1 + 0.5*(1.4-1)*OD.altdata(z,7)^2)^(1.4/0.4))*ones(m1,m1);%Stagnation pressure at compressor inlet
OD.Talt = OD.altdata(z,4)*(1 + 0.5*0.4*OD.altdata(z,7)^2)*ones(m1,m1); %stagnation temperature at compressor inlet
OD.P0 = OD.altdata(z,3)*ones(10,10);
OD.T0 = OD.altdata(z,4)*ones(10,10);
OD.P_non_perm = OD.Pratio.*OD.Palt; 
OD.velocity = OD.altdata(z,2); %Velocity at altitude
OD.dT_fc = 50*ones(m1,m1); %Maximum temperature differential, Kelvin
OD.asr = 0.15*ones(m1,m1); % Area specific resistance, ohm-cm^2
OD.P_fc = 1000*ones(m1,m2); %Operating pressure for SOFC
OD.T_fc = 1023*ones(m1,m2); %Inlet temperature for SOFC
OD.T_otm = 1035.5*ones(m1,m2); %Operating temperature for OTM
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
[param_od] = run_cycle_od_b(OD,options,z);
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
     [Weight] = weight_b(OD); 
OD.Mass(z,1:6) = [Weight.Total,Weight.turb,Weight.comp,Weight.hx,Weight.otm,Weight.sofc(1,1)]; %Total mass for a single power plant across a variety of conditions
end
l = 10*b;
m = l-9; 
PerformanceTable(m:l,1:18) = OD.Performance; 
Performance(:,:,b) = OD.Performance; 
end
[Battery] = battery(options);
Weight.Battery = Battery.Total/1260; 
Final.BatteryMass(b,1) = Weight.Battery; 
Final.TurbineMass(b,1) = max(OD.Mass(1:10,2));
Final.CompMass(b,1) = max(OD.Mass(1:10,3)); 
Final.HXMass(b,1) = max(OD.Mass(1:10,4));
Final.OTMMass(b,1) = max(OD.Mass(1:10,5));
Final.SOFCMass(b,1) = max(OD.Mass(1:10,6)); 
Final.HTSMMass(b,1) = 106000/20; 
Final.Mass(b,1) = (Final.TurbineMass(b,1) + Final.CompMass(b,1) + Final.HXMass(b,1) + Final.OTMMass(b,1) + Final.SOFCMass(b,1) + Final.BatteryMass(b,1)); %Find total mass required for power plant index b to perform across all power demands and ambient conditions with 10% for hardware
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
TSFCcruise = OD.Performance(10,5)*1000000/param.Thr_cruise; %g/kNs
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
Final.UseableFuel(b,1) = matchfuelmass; 
% OD.Cruise_fuel = W1b - OD.FB_to - OD.FB_climb - OD.FuelReserve; %Total fuel available after takeoff and climb
% OD.Cruise_time = (OD.Cruise_fuel)/OD.Performance(10,5); %Cruise time with useable fuel mass
% OD.Cruise_distance_m = OD.Cruise_time*OD.altdata(10,2); %Total Cruise distance, meters
% OD.Cruise_distance_nm = OD.Cruise_distance_m/1852; %Total Cruise distance, nautical miles 
% OD.Cruise_nm_matchstandard = 7260; %standard range of 747
% OD.Cruise_matchtime = 7260*1852/OD.altdata(10,2); 
% OD.Cruise_matchfuelmass = OD.Cruise_matchtime*OD.Performance(10,5); %Fuel required to match standard 747 cruise range
% OD.PayloadatRange = (101000 + 112760) - OD.Cruise_matchfuelmass - Design.Penalty - OD.FuelReserve; %Available payload after mass penalties 
% OD.tsfc_to =Design.H2usedactual*1000000/param.TrTO; %g/kN*s
% OD.tsfc_cruise = OD.Performance(10,5)*1000000/param.Thr_cruise; %g/kN*s

end

Final.Table = [Final.Payload, Final.TurbineMass,Final.CompMass,Final.HXMass,Final.OTMMass, Final. SOFCMass,Final.HTSMMass,Final.BatteryMass,Final.MassPenalty, Final.Mass,Final.TSFC_to,Final.TSFC_cruise,Final.UseableFuel]
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