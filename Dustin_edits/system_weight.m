function weight = system_weight(options,FC,OTM,HL,A1)
weight.sofc = options.SOFC_area.*options.sofc_specific_mass;
weight.otm = options.OTM_area.*options.OTM_specific_mass;

mass_flow = 28.84*net_flow(A1);
PratioComp = sqrt(options.P_non_perm./A1.P); %Pressure ratio required for a two stage turbine and two stage compressor
relation = 0.2131*PratioComp.^2 -2.508.*PratioComp + 16.901; %Compressor mass based on pressure ratio and mass flow rate of 1.16 kg/s, from NASA paper
intake_mass = mass_flow.*relation./1.16; %Two compressor stages with mass scale 
blower_mass = ((0.231*(1000/980)^2 -2.508*(1000/980) + 16.901)/1.16).*(FC.H2_supply + FC.H2O_supply)*2; %Extrapolated blower mass based on pressure ratio and mass flow of h2 into fuel cell 
weight.comp = 2*intake_mass + blower_mass; %Compressor mass

PratioTurb = sqrt(options.P_non_perm./A1.P);
turbine_mass = -0.381*PratioTurb.^2 + 5.5*PratioTurb + 1.9167; %Relation for Turbine mass as a function of pressure ratio, 
weight.turb = 2*turbine_mass.*mass_flow./1.16; 

weight.hx = (OTM.Q_out + OTM.heat_added + HL.Qremove_fuel + HL.Q_preheat + HL.Q_removed + OTM.Q_oxygen_HX)./options.heat_exchange_power_den; %Heat exchanger weight

weight.motor = options.motor_eff.*(FC.Power + OTM.net_work + HL.blower_work)./options.motor_power_den; %Weight of propulstion motors
weight.propulsor = options.propulsor_weight; 
end