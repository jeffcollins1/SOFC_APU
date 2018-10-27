function [Outlet,Prop] = compressor(Inlet,P_out,eff)
[n1,n2] = size(Inlet.T); 

%MM = molarmass(Inlet)*ones(n1,n1);

Ru = 8.314*ones(n1,n2); %Universal Gas constant, kJ/(kg*K)
%gamma = property(Inlet,'C','kJ/(kmol K)')./property(Inlet,'O','kJ/(kmol K)');
[~,H_in] = enthalpy(Inlet); %property(Inlet,'h','kJ');
S_in = entropy(Inlet); %property(Inlet,'s','kJ/K');
cp = SpecHeat(Inlet);
% H_in = ones(10,1)*H_1';
% S_in = ones(10,1)*S_1';
% cp = ones(10,1)*cp1';
gamma = cp./(cp - Ru); 
Ideal = Inlet;
Ideal.P = P_out;
Ideal.T = (Inlet.T.*(P_out./Inlet.P).^((gamma-1)./gamma));

% Ideal_2 = Ideal;
% Ideal_2.T = Ideal.T+1;
% dT_dS = 1./(entropy(Ideal_2) - entropy(Ideal));  %(property(Ideal_2,'s','kJ/K')-property(Ideal,'s','kJ/K'));
% %dT_dS = ones(10,1)*dT_dS1';
% error = 1;
% while max(abs(error))>1e-3
%    
%     error = (S_in - entropy(Ideal))./S_in; %property(Ideal,'s','kJ/K'))./S_in;
%     Ideal.T = Ideal.T + 2*error.*dT_dS; %.*S_in./MM; % .*error.*S_in;
% end
[~,H_ideal] = enthalpy(Ideal); %property(Ideal,'h','kJ');
H_out = H_in + 1./eff.*(H_ideal - H_in);
Outlet = Ideal;
Outlet.T = Ideal.T + (H_out - H_ideal)./cp;
Outlet.T = find_T(Outlet, H_out);
Prop.work = H_in - H_out;
Prop.pressure_ratio = P_out./Inlet.P;
Prop.eff = eff;
Prop.mass_flow = mass_flow(Inlet);
end

