function [outlet,prop] = expander(Inlet,P_out,eff)
R = 8.314*ones(10,10); %kJ/kmol*K
cp = SpecHeat(Inlet);%property(Inlet,'C','kJ/(kmol K)');
gamma = cp./(cp-R);
H_in = enthalpy(Inlet); %property(Inlet,'h','kJ');
S_in = entropy(Inlet); %property(Inlet,'s','kJ/K');
Ideal = Inlet;
Ideal.P = P_out;
Ideal.T = Inlet.T.*(P_out./Inlet.P).^((gamma-1)./gamma);

% Ideal_2 = Ideal;
% Ideal_2.T = Ideal.T+1;
% dT_dS = 1./(entropy(Ideal_2) - entropy(Ideal));%(property(Ideal_2,'s','kJ/K')-property(Ideal,'s','kJ/K'));
% 
% error = 1;
% while max(max(abs(error)))>2.5e-2
%     error = (S_in - property(Ideal))./S_in;  %property(Ideal,'s','kJ/K'))./S_in;
%     Ideal.T = Ideal.T + dT_dS.*error.*S_in;
% end
H_ideal = enthalpy(Ideal); %property(Ideal,'h','kJ');
H_out = H_in - eff.*(H_in - H_ideal);
Cp = SpecHeat(Ideal); %property(Ideal,'C','kJ/K');
outlet = Ideal;
outlet.T = Ideal.T + (H_out - H_ideal)./Cp;
outlet.T = find_T(outlet, H_out);
prop.work = H_in - H_out;
prop.pressure_ratio = P_out./Inlet.P;
prop.eff = eff;
prop.mass_flow = mass_flow(outlet);
end

