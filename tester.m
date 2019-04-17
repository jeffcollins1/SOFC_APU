A1.T = 400*ones(10,10);
A1.P = 100*ones(10,10); 
A1.N2 = 0.79*ones(10,10);
A1.O2 = 0.21*ones(10,10); 
MM = molarmass(A1)
specheat = SpecHeat(A1);
[~,H] = enthalpy(A1);
ss = entropy(A1); 