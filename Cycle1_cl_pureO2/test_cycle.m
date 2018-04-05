%%  Test Cycle
%Create options to test
options.height = 1000; %Altitude, meters
options.OTM_area = 2.5e3; %membrane area in m^2
options.SOFC_area = 2e3; %SOFC active area
options.dT_fc = 50; %Maximum temperature differential, Kelvin
options.asr = 0.2; % Area specific resistance, ohm-cm^2
options.P_fc = 1000; %Operating pressure for SOFC
options.T_fc = 1023; %Inlet temperature for SOFC
options.T_otm = options.T_fc+25; %Operating temperature for OTM
options.T_oxygen_pump = 323; %Inlet temperature of vacuum pump
options.spu = 0.1; 
options.steamratio = 0.05; %Percentage of humidification at fuel inlet
options.height = 100; %Altitude, meters
options.velocity = 100; %Velocity, m/s 
options.P_non_perm = 1000; %Constant permeate pressure for OTM, kPa
options.C1_eff = 0.80; %Mechanical efficiency of compressor 1
options.airflow = 1;
options.P_perm = 50; %Pressure of OTM oxygen stream, kPa
options.OTM_perc_theoretical = 0.5; %Actual percentage of theoretical O2 recovered
options.T1_eff = 0.88; %Mechanical efficiency of turbine
options.C2_eff = 0.80; %Mechanical efficiency of compressor 2
mission = [];

Param = run_cycle(options,mission);
