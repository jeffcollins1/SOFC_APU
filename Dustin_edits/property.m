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
        X(A.(F{i})==0) = 1;
        n = zeros(size(A.T));
        for j = 1:1:length(A.T(1,:))
            n(:,j) = refproparray(out,'T',A.T(:,j),'P',X(:,j).*A.P(:,j),REFPROPnom{s});
        end
        switch unit
            case {'kJ/kg' , 'kJ/(kg K)'}
                new_val = n./1000*A.(F{i})./NetFlow;
            case {'kJ/kmol' , 'kJ/(kmol K)'}
               new_val = n./1000*MM(s).*X;
            case {'kJ' , 'kJ/K'}
               new_val = n./1000*MM(s).*A.(F{i});
        end
        new_val(A.(F{i})==0) = 0;
        Val = Val + new_val;
    end
end
end%Ends function property
