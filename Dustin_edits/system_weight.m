function weight = system_weight(options,FC,turbo,OTM,HL,HX,HX2)
weight.sofc = options.SOFC_area.*options.sofc_specific_mass;
net_power = FC.Power;
weight.comp = 0;
weight.turb = 0;
for i = 1:1:length(turbo)
    net_power = net_power + turbo{i}.work;
    if turbo{i}.work(1,1)>0
        weight.turb = weight.turb + turbo{i}.mass_flow/1.16.*(-0.381*(turbo{i}.pressure_ratio.^.5).^2 + 5.5*(turbo{i}.pressure_ratio.^.5) + 1.9167);
    else
        weight.comp = weight.comp + turbo{i}.mass_flow/1.16.*(0.2131*(turbo{i}.pressure_ratio.^.5).^2 -2.508.*(turbo{i}.pressure_ratio.^.5) + 16.901); %Compressor mass based on pressure ratio and mass flow rate of 1.16 kg/s, from NASA paper
    end
end
weight.hx = (HL.Qremove_fuel + HL.Q_preheat + HL.Q_removed)./options.heat_exchange_power_den; %Heat exchanger weight

if ~isempty(OTM)
    weight.otm = options.OTM_area.*options.OTM_specific_mass;
    weight.hx = weight.hx + (OTM.Q_out + OTM.heat_added + OTM.Q_oxygen_HX)./options.heat_exchange_power_den; %Heat exchanger weight
end


weight.motor = options.motor_eff.*net_power./options.motor_power_den; %Weight of propulstion motors

% PratioTurb = sqrt(options.PR_comp);
% turbine_mass = -0.381*PratioTurb.^2 + 5.5*PratioTurb + 1.9167; %Relation for Turbine mass as a function of pressure ratio, 
% weight.turb = 2*turbine_mass.*mass_flow./1.16; 
weight.hx_fuel = HX2.fuel; 
weight.hx = HX2.oxycompressor + HX2.fuel + HX2.condenser + HX.oxygen + HX2.HP; %Heat exchanger weight

weight.propulsor = options.propulsor_weight; 
end