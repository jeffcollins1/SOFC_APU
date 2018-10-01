%%  Test Cycle
%Create options to test
n1 = 10; % number of points in test dimension 1
n2 = 10; % number of points in test dimension 2
options.height = 10000*ones(n1,n2); %Altitude, meters
options.OTM_area = 3e3*ones(n1,n2); %membrane area in m^2 per kmol airflow
options.SOFC_area = 3e3*ones(n1,n2); %membrane area in m^2 per kmol airflow
options.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
options.asr = 0.15*ones(n1,n2); % Area specific resistance, ohm-cm^2
options.P_fc = 1000*ones(n1,n2); %Operating pressure for SOFC
options.T_fc = 1023*ones(n1,n2); %Inlet temperature for SOFC
options.T_otm = options.T_fc; %Operating temperature for OTM
options.j0_otm = 7*ones(n1,n2); %Nominal oxygen flux through OTM NmL/cm^2*min
options.P0_otm = 10*ones(n1,n2); %Nominal oxygen pressure ratio across OTM
options.T_oxygen_pump = 323*ones(n1,n2); %Inlet temperature of vacuum pump
options.T_motor = 77*ones(n1,n2); %temperture of H2 gas after cooling superconducting motors
options.spu = 0.2*ones(n1,n2); 
options.steamratio = 0.01*ones(n1,n2); %Percentage of humidification at fuel inlet
options.P_non_perm = ones(n1,1)*linspace(1350,3500,n2); %Range of intake pressures for OTM, kPa
options.C1_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 1
options.airflow = ones(n1,n2); %Initial airflow, kmol/s
options.P_perm = logspace(-1.3,-.01,n1)'*(0.21*options.P_non_perm(1,:)); %Pressure of OTM oxygen stream, kPa; 
options.T1_eff = 0.88*ones(n1,n2); %Mechanical efficiency of turbine
options.C2_eff = 0.80*ones(n1,n2); %Mechanical efficiency of compressor 2
options.Blower_eff = 0.5*ones(n1,n2); %efficiency of blower
options.Blower_dP = 20*ones(n1,n2); %Pressure rise in blower in kPa

mission = [];
tic
[Weight,param] = run_cycle(options,mission);
toc
