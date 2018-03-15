function obj = refprop2(prop1, obj1, prop2, obj2, fluid, varargin)

%% Get nval

nval1 = obj1.nominalValue;
nval2 = obj2.nominalValue;

%% Get aeval

aeval1 = obj1.absoluteError;
aeval2 = obj2.absoluteError;
lval1 = nval1 - aeval1/2;
uval1 = nval1 + aeval1/2;
lval2 = nval2 - aeval2/2;
uval2 = nval2 + aeval2/2;

indicUncert1 = sum(abs(aeval1(:))) == 0;
indicUncert2 = sum(abs(aeval2(:))) == 0;
indicUncert3 = isempty(aeval1);
indicUncert4 = isempty(aeval2);

if or(and(indicUncert1, indicUncert2), and(indicUncert3, indicUncert4));
    uncertaintyMode = 0;
else
    uncertaintyMode = 1;
end


obj = ComputeProperties(prop1, nval1, prop2, nval2, fluid);

propertyNames = fieldnames(obj);
propertyNames = propertyNames(not(cellfun(@(x) strcmp(x, 'dotM'), propertyNames)));

if uncertaintyMode == 1
    
    tmp01 = ComputeProperties(prop1, lval1, prop2, lval2, fluid);
    tmp02 = ComputeProperties(prop1, lval1, prop2, uval2, fluid);
    tmp03 = ComputeProperties(prop1, uval1, prop2, lval2, fluid);
    tmp04 = ComputeProperties(prop1, uval1, prop2, uval2, fluid);
    
    tmp01sat = ComputeProperties('P', tmp01.P.nominalValue, 'x', obj.x.nominalValue, fluid);
    tmp02sat = ComputeProperties('P', tmp02.P.nominalValue, 'x', obj.x.nominalValue, fluid);
    tmp03sat = ComputeProperties('P', tmp03.P.nominalValue, 'x', obj.x.nominalValue, fluid);
    tmp04sat = ComputeProperties('P', tmp04.P.nominalValue, 'x', obj.x.nominalValue, fluid);
    
    for k = 1:size(propertyNames)
        tmp01.(propertyNames{k}).nominalValue(tmp01.x.nominalValue ~= obj.x.nominalValue) = tmp01sat.(propertyNames{k}).nominalValue(tmp01.x.nominalValue ~= obj.x.nominalValue);
        tmp02.(propertyNames{k}).nominalValue(tmp02.x.nominalValue ~= obj.x.nominalValue) = tmp02sat.(propertyNames{k}).nominalValue(tmp02.x.nominalValue ~= obj.x.nominalValue);
        tmp03.(propertyNames{k}).nominalValue(tmp03.x.nominalValue ~= obj.x.nominalValue) = tmp03sat.(propertyNames{k}).nominalValue(tmp03.x.nominalValue ~= obj.x.nominalValue);
        tmp04.(propertyNames{k}).nominalValue(tmp04.x.nominalValue ~= obj.x.nominalValue) = tmp04sat.(propertyNames{k}).nominalValue(tmp04.x.nominalValue ~= obj.x.nominalValue);
        values = [tmp01.(propertyNames{k}).nominalValue ...
            tmp02.(propertyNames{k}).nominalValue ...
            tmp03.(propertyNames{k}).nominalValue ...
            tmp04.(propertyNames{k}).nominalValue];
        
        propLowerValues = min(values, [], 2);
        propUpperValues = max(values, [], 2);
        propAbsoluteError = abs(propUpperValues - propLowerValues);
        obj.(propertyNames{k}).absoluteError = propAbsoluteError;
    end
    
    obj.(prop1).absoluteError = aeval1;
    obj.(prop2).absoluteError = aeval2;
end

%% Base function

    function results = ComputeProperties(prop1, value1, prop2, value2, fluid)
        
        %% Conversion of the prop1 and prop2 strings

        % Formatting of the refpropm properties required

        conv_table = { ...
            'rho', 'D'; ... % Density [kg/m3];
            's', 'S'; ... % Specific entropy [J/(kg/K)];
            'h', 'H'; ... % Specific enthalpy [J/kg];
            'x', 'Q'; ... % Quality (vapor fraction) (kg/kg);
            'mu', 'V'; ... % Dynamic viscosity [Pa*s];
            'gamma', 'I'; ... % Surface tension [N/m];
            'lambda', 'L'; ... % Thermal conductivity [W/(m K)];
            'u', 'U'; ... % Specific internal energy [J/kg];
            'cp', 'C'; ... % Specific heat at constant pressure [J/(kg K)];
            'cv', 'O'; ... % Specific heat at constant volume [J/(kg K)];
            'kappa', 'K' ... % Ratio of specific heats (Cp/Cv) [-]
            };
        % Not changed:
        %   - P: Pressure [Pa];
        %   - T: Temperature [K];
        %   - A: Speed of sound [m/s];
        % Not implemented:
        %   - X: liquid and gas composition as the format of the output is not
        %        compliant with the current code base. It can be added by
        %        some additionnal work (transpose ?)
        
        % Conversion of the property strings prop1 and prop2
        for j = 1:size(conv_table, 1)
            if strcmp(prop1, conv_table{j, 1})
                prop1 = conv_table{j, 2};
            end
            if strcmp(prop2, conv_table{j, 1})
                prop2 = conv_table{j, 2};
            end
        end
        
        % Units conversion for refpropm input
        if strcmp(prop1, 'P')
            value1 = value1 ./ 1e3;
        end
        
        if strcmp(prop2, 'P')
            value2 = value2 ./ 1e3;
        end
        
        if numel(value1) ~= numel(value2)
            error('Not the same number of elements.') %%% FIXME
        end
        
        if size(fluid, 2) == 2
            fluid = [fluid(:, 1) repmat(fluid(:,2), 1, size(value1, 1))];
        end
        
        
        resultsStruct = struct;
        for i = 1:size(value1, 1)
            % Generation of the fluid string
            compositionString = cell2mat(transpose(cellfun(@(x) ...
                [num2str(x) ', '], fluid(:,i+1), 'UniformOutput', 0)));
            compositionString = ['[' compositionString(1:end-2) ']'];
            components = fluid(:, 1);
            componentsString = cell2mat(transpose(cellfun(@(x) ...
                ['''' x ''', '], components, 'UniformOutput', 0)));
            fluidString = [componentsString, compositionString];
            
            refpropCommand = ['[resultsStruct.P(i, 1), ' ...
                'resultsStruct.T(i, 1), resultsStruct.rho(i, 1), ' ...
                'resultsStruct.A(i, 1), resultsStruct.s(i, 1), resultsStruct.h(i, 1), ' ...
                'resultsStruct.x(i, 1), resultsStruct.mu(i, 1), ' ...
                'resultsStruct.lambda(i, 1), ' ...
                'resultsStruct.u(i, 1), resultsStruct.cp(i, 1), ' ...
                'resultsStruct.cv(i, 1), resultsStruct.kappa(i, 1)] = ' ...
                'refpropm(''PTDASHQVLUCOK'', ''' prop1 ''', value1(i, 1), ''' ...
                prop2 ''', value2(i, 1), ' fluidString ');'];
            eval(refpropCommand);
        end
        
        props = fieldnames(resultsStruct);
        for l = 1:size(props, 1)
            if strcmp(props{l}, 'P')
                value = resultsStruct.(props{l}) * 1e3;
            else
                value = resultsStruct.(props{l});
            end
            results.(props{l}) = PhysicalValue();
            results.(props{l}).nominalValue = value;
            results.(props{l}).absoluteError = zeros(size(value));
        end
        
%         results.P.nominalValue = resultsStruct.P.*1e3;
%         results.T.nominalValue = resultsStruct.T;
%         results.rho.nominalValue = resultsStruct.rho;
%         results.s.nominalValue = resultsStruct.s;
%         results.h.nominalValue = resultsStruct.h;
%         results.x.nominalValue = resultsStruct.x;
%         results.A.nominalValue = resultsStruct.A;
%         results.mu.nominalValue = resultsStruct.mu;
%         % results.gamma.nominalValue = resultsStruct.gamma;
%         results.lambda.nominalValue = resultsStruct.lambda;
%         results.u.nominalValue = resultsStruct.u;
%         results.cp.nominalValue = resultsStruct.cp;
%         results.cv.nominalValue = resultsStruct.cv;
%         results.kappa.nominalValue = resultsStruct.kappa;
        
        results.x.nominalValue(results.x.nominalValue < 0) = 0;
        results.x.nominalValue(results.x.nominalValue > 1) = 1;
        
    end


end
