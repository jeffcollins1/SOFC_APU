function [Outlet,Work] = compressor(Inlet,P_out,eff)
gamma = property(Inlet,'C','kJ/(kmol K)')./property(Inlet,'O','kJ/(kmol K)');
H_in = property(Inlet,'h','kJ');
S_in = property(Inlet,'s','kJ/K');
Ideal = Inlet;
Ideal.P = P_out;
Ideal.T = Inlet.T.*(P_out./Inlet.P).^((gamma-1)./gamma);

Ideal_2 = Ideal;
Ideal_2.T = Ideal.T+1;
dT_dS = 1./(property(Ideal_2,'s','kJ/K')-property(Ideal,'s','kJ/K'));

error = 1;
while abs(error)>1e-3
    error = (S_in - property(Ideal,'s','kJ/K'))./S_in;
    Ideal.T = Ideal.T + dT_dS.*error.*S_in;
end
H_ideal = property(Ideal,'h','kJ');
H_out = H_in + 1./eff.*(H_ideal - H_in);
Cp = property(Ideal,'cp','kJ/K');
Outlet = Ideal;
Outlet.T = Ideal.T + (H_out - H_ideal)./Cp;
Outlet.T = find_T(Outlet, H_out);
Work = H_in - H_out;


