function S = enthalpy(A)
S = 0;
if isfield(A,'O2')
    S = S + refpropm('H','T',A.T,'P',A.P,'OXYGEN')/1000; %Enthalpy of permeate O2 stream
end
if isfield(A,'N2')
    S = S + refpropm('H','T',A.T,'P',A.P,'NITROGEN')/1000; %Enthalpy of permeate O2 stream
end
if isfield(A,'H2')
    S = S + refpropm('H','T',A.T,'P',A.P,'HYDROGEN')/1000; %Enthalpy of permeate O2 stream
end
if isfield(A,'H2O')
    S = S + refpropm('H','T',A.T,'P',A.P,'WATER')/1000; %Enthalpy of permeate O2 stream
end