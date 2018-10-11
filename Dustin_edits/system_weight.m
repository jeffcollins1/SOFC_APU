function weight = system_weight(options,FC,turbo,OTM,HX)
weight.sofc = options.SOFC_area.*options.sofc_specific_mass;
net_power = FC.Power;
weight.comp = 0;
weight.turb = 0;
for i = 1:1:length(turbo)
    net_power = net_power + turbo{i}.work;
    if turbo{i}.work(1,1)>0
        weight.turb = weight.turb + 2*turbo{i}.mass_flow/1.16.*(-0.381*(turbo{i}.pressure_ratio.^.5).^2 + 5.5*(turbo{i}.pressure_ratio.^.5) + 1.9167);
    else
        weight.comp = weight.comp + 2*turbo{i}.mass_flow/1.16.*(0.2131*(turbo{i}.pressure_ratio.^.5).^2 -2.508.*(turbo{i}.pressure_ratio.^.5) + 16.901); %Compressor mass based on pressure ratio and mass flow rate of 1.16 kg/s, from NASA paper
    end
end
weight.hx_fuel = HX.fuel.mass; 
weight.hx =  HX.fuel.mass + HX.condenser.mass; %Heat exchanger weight
if ~isempty(OTM)
    weight.otm = options.OTM_area.*options.OTM_specific_mass;
    weight.hx = weight.hx + HX.oxycompressor.mass + HX.oxygen.mass + HX.HP.mass; %Heat exchanger weight
end

weight.motor = options.motor_eff.*net_power./options.motor_power_den; %Weight of propulstion motors
weight.propulsor = options.propulsor_weight; 
end