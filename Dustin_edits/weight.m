function weight = weight(options,FC,OTM,HL)
T = 245196;% Thrust in N:  Thrust = Cd/Cl*Weight; cruise Cd = 0.036; Cl = .48; weight = 820,000lbs; thus 55,000lbs cruise thrust. Makes sens because max thrust at cruise is 100,000lbs.Range = 6100nautical miles, velocity = 635knts
P = T.*V./options.prop_eff;%propellor momentum theory at cruise http://164.100.133.129:81/econtent/Uploads/09-%20Ducted%20Fans%20and%20Propellers%20%5BCompatibility%20Mode%5D.pdf ,  http://web.mit.edu/16.unified/www/SPRING/systems/Lab_Notes/airpower.pdf
scale = P./options.motor_eff./(FC.Power + OTM.net_work + HL.blower_work);

weight.sofc = scale.*FC.Cells.*0.05508; %Weight per cell, kg
weight.otm = scale.*options.OTM_area.*0.048907.*10000./81; %Weight per cell of each 81 cm^2 OTM membrane, kg
weight.battery = battery_kJ/1260; %battery weight required to assist with takeoff assuming battery energy storage of 1260 kJ/kg;

mass_flow = 28.84*scale; 
A1 = std_atmosphere(options.height);%Ambient conditions as a function of altitude
PratioComp = sqrt(options.P_non_perm./A1.P); %Pressure ratio required for a two stage turbine and two stage compressor
PratioTurb = sqrt(options.P_non_perm./A1.P);
relation = 0.2131*PratioComp.^2 -2.508.*PratioComp + 16.901; %Compressor mass based on pressure ratio and mass flow rate of 1.16 kg/s, from NASA paper
intake_mass = mass_flow.*relation./1.16; %Two compressor stages with mass scale 
blower_mass = ((0.231*(1000/980)^2 -2.508*(1000/980) + 16.901)/1.16).*(FC.H2_supply + FC.H2O_supply)*2; %Extrapolated blower mass based on pressure ratio and mass flow of h2 into fuel cell 
weight.comp = 2*intake_mass + blower_mass; %Compressor mass

turbine_mass = -0.381*PratioTurb.^2 + 5.5*PratioTurb + 1.9167; %Relation for Turbine mass as a function of pressure ratio, 
weight.turb = 2*turbine_mass.*mass_flow./1.16; 
Pden_hx = 15; %Power density of heat exchangers, kW/kg
weight.hx = scale.*(OTM.Q_out + OTM.heat_added + HL.Qremove_fuel + HL.Q_preheat + HL.Q_removed + OTM.Q_oxygen_HX)/Pden_hx; %Heat exchanger weight
Pden_motor = 24; %Power density of HTSM
weight.motor = P./Pden_motor; %Weight of propulstion motors

weight.propulsor = 0.3*4*(9670/2.2); %Weight penalty of replacing power plants of 4 RR RB-211 engines while keeping the propulsor
end