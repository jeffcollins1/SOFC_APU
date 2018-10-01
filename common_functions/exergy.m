function A = exergy(A,T0,P0)
A.H = property(A,'h','kJ');
A.S = property(A,'s','kJ/K');
A0 = A;
A0.T = T0;
A0.P = P0;
A0.H = property(A0,'h','kJ');
A0.S = property(A0,'s','kJ/K');

A.X = A.H - A0.H - A0.T.*(A.S-A0.S)  ; % - A0.T*(A.S - A0.S) %A0.H - A0.T*(A.S - A0.S);
end