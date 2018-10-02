%%  Test Cycle
%Create options to test
n1 = 10; % number of points in test dimension 1
n2 = 10; % number of points in test dimension 2
options.OTM_area = 2e3*ones(n1,n2); %membrane area in m^2 per kmol airflow
options.SOFC_area = linspace(1e3,5e3,n1)'*ones(1,n2); %membrane area in m^2 per kmol airflow
options.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
options.asr = 0.15*ones(n1,n2); % Area specific resistance, ohm-cm^2
options.P_fc = 1000*ones(n1,n2); %Operating pressure for SOFC
options.T_fc = 1023*ones(n1,n2); %Inlet temperature for SOFC
options.T_otm = options.T_fc; %Operating temperature for OTM
options.j0_otm = 7*ones(n1,n2); %Nominal oxygen flux through OTM NmL/cm^2*min
options.P0_otm = 2.1*ones(n1,n2); %Nominal oxygen pressure ratio across OTM (total pressure ratio *.21)
options.T_oxygen_pump = 323*ones(n1,n2); %Inlet temperature of vacuum pump
options.T_motor = 77*ones(n1,n2); %temperture of H2 gas after cooling superconducting motors
options.spu = 0.2*ones(n1,n2); 
options.steamratio = 0.01*ones(n1,n2); %Percentage of humidification at fuel inlet
options.PR_comp = ones(n1,1)*linspace(15,60,n2); %Range of intake pressures for OTM, kPa
options.P_non_perm = 50*ones(n1,n2); %Range of intake pressures for OTM, kPa
options.C1_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 1
options.airflow = ones(n1,n2); %Initial airflow, kmol/s
options.P_perm = 50*ones(n1,n2); %Pressure of OTM oxygen stream, kPa; 
options.T1_eff = 0.88*ones(n1,n2); %Mechanical efficiency of turbine
options.C2_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 2
options.Blower_eff = 0.5*ones(n1,n2); %efficiency of blower
options.Blower_dP = 20*ones(n1,n2); %Pressure rise in blower in kPa
options.TO_weight = 820000*ones(n1,n2);%take-off weight in lbs
options.Lift_2_Drag = 13.5*ones(n1,n2);%lift to drag ratio
options.prop_eff = 0.95*ones(n1,n2);%propulsor efficiency
options.motor_eff = 0.984*ones(n1,n2);%motor efficiency

band = [0;500;4500;9500;19500;29500;39500;]/3.1;%altitude converted to m
mission.alt = [(band(2:end)+band(1:end-1))/2;band(end)];
mission.duration = [(band(2:end)+band(1:end-1))/3.54;11];%assume a steady climb rate of 3.54m/s = climb to 39500 ft in 1 hr, then 11 hours cruise
mission.mach_num = min(0.8,0.5+0.3*mission.alt/9000);
tic
param = run_cycle(options,mission);
toc
