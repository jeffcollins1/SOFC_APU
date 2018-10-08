function [outlet,prop] = compressor(Inlet,P_out,eff)
gamma = property(Inlet,'C','kJ/(kmol K)')./property(Inlet,'O','kJ/(kmol K)');
H_in = property(Inlet,'h','kJ');
S_in = property(Inlet,'s','kJ/K');
Ideal = Inlet;
Ideal.P = P_out;
Ideal.T = (Inlet.T.*(P_out./Inlet.P).^((gamma-1)./gamma));

Ideal_2 = Ideal;
Ideal_2.T = Ideal.T+1;
dT_dS = 1./(property(Ideal_2,'s','kJ/K')-property(Ideal,'s','kJ/K'));

error = 1;
while max(abs(error))>1e-3
    error = (S_in - property(Ideal,'s','kJ/K'))./S_in;
    Ideal.T = Ideal.T + dT_dS.*error.*S_in;
end
H_ideal = property(Ideal,'h','kJ');
H_out = H_in + 1./eff.*(H_ideal - H_in);
Cp = property(Ideal,'cp','kJ/K');
outlet = Ideal;
outlet.T = Ideal.T + (H_out - H_ideal)./Cp;
outlet.T = find_T(outlet, H_out);
prop.work = H_in - H_out;
prop.pressure_ratio = P_out./Inlet.P;
prop.eff = eff;
prop.mass_flow = mass_flow(outlet);
end

