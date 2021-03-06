function n = mass_flow(A)
CESInom = {'O2';'N2';'H2';'H2O';};
mass = [32;28;2;18;];
n = 0;
F = fieldnames(A);
for i = 1:1:length(F)
    s = nonzeros((1:1:length(CESInom))'.*strcmp(F{i},CESInom));
    if ~isempty(s)
        n = n + A.(F{i})*mass(i);
    end
end