function Val = property(A,out,unit)
%usable for cv, cp, h, s, 
CESInom = {'O2';'N2';'H2';'H2O';};
REFPROPnom = {'OXYGEN';'NITROGEN';'HYDROGEN';'WATER'};
MM = [32,28,2,18]; %molar mass
Val = 0;
NetFlow = net_flow(A);
F = fieldnames(A);
for i = 1:1:length(F)
    s = nonzeros((1:1:length(CESInom))'.*strcmp(F{i},CESInom));
    if ~isempty(s)
        X = A.(F{i})./NetFlow;
        n = refproparray(out,'T',A.T,'P',X.*A.P,REFPROPnom{s});
        switch unit
            case {'kJ/kg' , 'kJ/(kg K)'}
                Val = Val + n./1000*A.(F{i})./NetFlow;
            case {'kJ/kmol' , 'kJ/(kmol K)'}
                Val = Val + n./1000*MM(s)*X;
            case {'kJ' , 'kJ/K'}
                Val = Val + n./1000*MM(s)*A.(F{i});
        end
    end
end
