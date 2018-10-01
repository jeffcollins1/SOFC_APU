OD.SOFC_area = param.FC_area_fixed*ones(10,10); 
Pair_int_max = param.P_int_max;
OD.air_in = param.air_in*ones(10,10); 
OD.P_int = Pair_int*ones(10,1)*linspace(100,Pair_int_max,10); 

OD.height = ones(10,1)*linspace(0,10000,10); %Altitude, meters
OD.thrust_climb = 12726*ones(n1,n2); %Required thrust, lbs 
 
OD.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
OD.asr = 0.2*ones(n1,n2); % Area specific resistance, ohm-cm^2
OD.P_fc = 1000*ones(n1,n2); %Operating pressure for SOFC
OD.T_fc = 1023*ones(n1,n2); %Inlet temperature for SOFC
OD.T_otm = OD.T_fc; %Operating temperature for OTM
OD.T_oxygen_pump = 323*ones(n1,n2); %Inlet temperature of vacuum pump
OD.T_motor = 77*ones(n1,n2); %temperture of H2 gas after cooling superconducting motors
OD.spu = 0.1*ones(n1,n2); 
OD.steamratio = 0.05*ones(n1,n2); %Percentage of humidification at fuel inlet
OD.velocity_climb = (ones(n1,1)*linspace(0,70,n2)); %Velocity, m/s 
OD.P_non_perm = ones(n1,1)*linspace(600,3000,n2); %Range of intake pressures for OTM, kPa
OD.C1_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 1
OD.airflow = 1*ones(n1,n2);
OD.P_perm = 75*ones(n1,n2); %linspace(50,250,n1)'*ones(1,n2); %Pressure of OTM oxygen stream, kPa; 
OD.OTM_perc_theoretical = 0.7*ones(n1,n2); %Actual percentage of theoretical O2 recovered
OD.T1_eff = 0.88*ones(n1,n2); %Mechanical efficiency of turbine
OD.C2_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 2
OD.Blower_eff = 0.5*ones(n1,n2); %efficiency of blower
OD.Blower_dP = 20*ones(n1,n2); %Pressure rise in blower in kPa
OD.T0= -0.0065*OD.height + 14.987 + 273.1; %Ambient Temperature Kelvin;
OD.P0 = 107*exp(-0.0001*OD.height)- 10; %Ambient Pressure as a function of altitude
mission = [];
tic
param_od = run_cycle_od(OD,param_od);
toc
