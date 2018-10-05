%%  Test Cycle
%Create options to test
n1 = 10; % number of points in test dimension 1
n2 = 10; % number of points in test dimension 2
options.OTM_area = 2e3*ones(n1,n2); %membrane area in m^2 per kmol airflow
options.SOFC_area = linspace(1e3,5e3,n1)'*ones(1,n2); %membrane area in m^2 per kmol airflow
options.dT_fc = 50*ones(n1,n2); %Maximum temperature differential, Kelvin
options.asr = 0.3*ones(n1,n2); % Area specific resistance, ohm-cm^2
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
options.TO_weight = 871200/2.2*ones(n1,n2);%take-off mass in kg
options.Lift_2_Drag = 13.5*ones(n1,n2);%lift to drag ratio
options.prop_eff = 0.95*ones(n1,n2);%propulsor efficiency
options.motor_eff = 0.984*ones(n1,n2);%motor efficiency
options.air_frame_weight = options.TO_weight - 101000 - 112760 - 4*9670/2.2;%airframe mass in kg: Max Fuel = 101000 kg.  Max payload = 112760 kg.  4 engines, each 9670lb
options.propulsor_weight = 0.3*4*(9670/2.2)*ones(n1,n2); %Weight propulsor portion (30%) of 4 RR RB-211 engines
options.motor_power_den = 24*ones(n1,n2); %Power density of HTSM
options.OTM_specific_mass = 0.048907*10000/81*ones(n1,n2); %Weight per m^2 OTM membrane, kg:  assumes 0.048907kg/ 81cm^2 cell
options.sofc_specific_mass = 0.05508*10000/81*ones(n1,n2); %Weight per m^2, kg:  assumes 0.05508kg/ 81cm^2 cell
%% should edit to be fcn of minimum temperature difference
options.heat_exchange_power_den = 15*ones(n1,n2); %Power density of heat exchangers, kW/kg

band = [0;500;4500;9500;19500;29500;39500;]/3.1;%altitude converted to m
mission.alt = [(band(2:end)+band(1:end-1))/2;band(end)];
mission.duration = [(band(2:end)-band(1:end-1))/7.5/3600;14.38];%assume a steady climb rate of 7.5m/s , then 14.38 hours cruise
mission.mach_num = min(0.8,0.5+0.3*mission.alt/9000);
tic
param = run_cycle(options,mission);
toc
payload = param.weight.payload;
payload(payload<0.8*mean(mean(param.weight.payload))) = nan;
ax1 = surf(options.PR_comp,param.i_den,payload);