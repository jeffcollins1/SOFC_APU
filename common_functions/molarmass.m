function Val = molarmass(A)
%usable for cv, cp, h, s, 
CESInom = {'O2';'N2';'H2';'H2O';};
%REFPROPnom = {'OXYGEN';'NITROGEN';'HYDROGEN';'WATER'};
MM = [32,28,2,18]; %molar mass
Val = 0;
NetFlow = net_flow(A);
F = fieldnames(A);
mass_species = zeros(1,length(F)); 

for i = 1:1:length(F)
    s = nonzeros((1:1:length(CESInom))'.*strcmp(F{i},CESInom));
    if ~isempty(s)
        X = A.(F{i})./NetFlow;
        X(A.(F{i})==0) = 1;
        n = zeros(length(F),length(F));
        for j = 1:1:length(F)
            n(i,j) = X(1)*MM(s);
           mass_species(i,j) = n(i,j); 
        end
       
    end
    
mass_mixture = mass_species(1:end,1)
end
 Val = sum(mass_mixture); 
end%Ends function property
