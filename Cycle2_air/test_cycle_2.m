%%  Test Cycle
%Create options to test
options.height = 100; %Altitude, meters
options.OTM_area = 2.5e3; %membrane area in m^2
options.SOFC_area = 4e3; %SOFC active area
options.dT_fc = 50; %Maximum temperature differential, Kelvin
options.asr = 0.2; % Area specific resistance, ohm-cm^2
options.P_fc = 1000; %Operating pressure for SOFC
options.T_fc = 1023; %Inlet temperature for SOFC
options.T_otm = options.T_fc+25; %Operating temperature for OTM
options.T_oxygen_pump = 323; %Inlet temperature of vacuum pump
options.T_motor = 77; %temperture of H2 gas after cooling superconducting motors
options.spu = 0.1; 
options.steamratio = 0.05; %Percentage of humidification at fuel inlet
options.velocity = 100; %Velocity, m/s 
options.P_non_perm = 3000; %Constant permeate pressure for OTM, kPa
options.C1_eff = 0.80; %Mechanical efficiency of compressor 1
options.airflow = 1;
options.P_perm = 100; %Pressure of OTM oxygen stream, kPa
options.OTM_perc_theoretical = 0.7; %Actual percentage of theoretical O2 recovered
options.T1_eff = 0.88; %Mechanical efficiency of turbine
options.C2_eff = 0.80; %Mechanical efficiency of compressor 2
options.Blower_eff = 0.5; %efficiency of blower
options.Blower_dP = 20; %Pressure rise in blower in kPa
options.T0= 298; %Ambient Temperature Kelvin;
options.P0 = 100; %Ambient Pressure as a function of altitude
options.O2ref = 0.2062; % molar fraction of oxygen in reference atmosphere
options.O2Xref = 3.9; %kJ/mol
options.N2ref = 0.7651;
options.H2Oref = 1.3;% kJ/mol
options.H2Oref = 0.0190;
options.rel_hum_ref = 0.6;
mission = [];

param = run_cycle_2(options,mission);
