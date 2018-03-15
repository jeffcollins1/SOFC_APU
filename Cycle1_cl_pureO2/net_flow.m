function n = net_flow(A)
n = 0;
if isfield(A,'O2')
    n = n + A.O2;
end
if isfield(A,'N2')
    n = n + A.N2;
end
if isfield(A,'H2')
    n = n + A.H2;
end
if isfield(A,'H2O')
    n = n + A.H2O;
end