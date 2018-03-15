function results = refprop(spec1, value1, spec2, value2, ...
    fluidstruct, varargin)
%REFPROP computes the thermophysical fluid properties of various pure
%    fluids or mixtures of pure fluids. This tool is a backend to the NIST
%    REFPROP MATLAB function 'refpropm'. The later needs to be present in
%    your MATLAB PATH and be functional for this function to work properly.
%    
%    USAGE: results = refprop(spec1, value1, spec2, value2, fluidstruct);
%
%    Available properties:
%    - 'A': Speed of sound [m/s];
%    - 'cp': Specific heat at constant pressure [J/(kg K)];
%    - 'cv': Specific heat at constant volume [J/(kg K)];
%    - 'gamma': Surface tension [N/m];
%    - 'h': Specific enthalpy [J/kg];
%    - 'kappa': Ratio of specific heats (Cp/Cv) [-]
%    - 'mu': Dynamic viscosity [Pa*s];
%    - 'P': Pressure [Pa];
%    - 'rho': Density [kg/m3];
%    - 's': Specific entropy [J/(kg/K)];
%    - 'T': Temperature [K];
%    - 'u': Specific internal energy [J/kg];
%    - 'x': Quality (vapor fraction) (kg/kg);
%    - 'lambda': Thermal conductivity [W/(m K)];
%
%    spec1: 'T', 'P', 'h' or 'rho';
%    spec2: 'P', 'rho', 'h', 's', 'u' or 'x';
%
%    fluidstruct is a cell formatted as follow:
%    - first column: fluid name (strings)
%    - all the next columns: mass fractions in the mixture [0:1];
%    
%    Example: fluidstruct = [{'R134a';'R125';'R32'},{0.52;0.25;0.23}]
%
%    OPTIONS:
%    - 'Uncertainties' ['yes'/'no'] - This option depends on LimitsVal1 and 
%      limitsval2.
%    - 'LimitsVal1' [[lower_value(s) upper_value(s)]];
%    - 'LimitsVal2' [[lower_value(s) upper_value(s)]];
%    - 'Properties' string containing any tag from the properties list
%      above, separated with commas (no space). Example: 'P,T,h,s'
%    - 'QuantityQualifier' ['nval' (nominal value) is the default, but you
%      can choose the quantifier you want.
%    - 'LowerValueQualifier' ['lval' (lower value) is the default, but you
%      can choose the quantifier you want.
%    - 'UpperValueQualifier' ['uval' (upper value) is the default, but you
%      can choose the quantifier you want.
%
%    MAIN QUANTIFIERS:
%    - 'mval': Measured value(s), not calibrated
%    - 'nval': Nominal value(s), calibrated
%    - 'lval': Lower value(s), calibrated
%    - 'uval': Upper value(s), calibrated
%    - 'aeval': absolute error(s), calibrated
%    - 'reval': relative error(s), calibrated
%
%    SIMPLE EXAMPLES:
%    R407C = [{'R134a';'R125';'R32'}, {0.52;0.25;0.23}];
%    comb_gas = [{'Oxygen';'Nitrogen';'CO2'}, {0.19;0.80;0.01}, ...
%        {0.18;0.80;0.02}, {0.17;0.80;0.03}];
%    R134a = {'R134a', 1};
%    test01 = refprop('T', [310;300], 'P', [2e5;1.8e5], R407C);
%    test02 = refprop('T', 310, 'P', [2e5;1.8e5], R407C);
%    test03 = refprop('T', [310;300], 'P', 2e5, R407C);
%    test04 = refprop('T', [310;300;350], 'P', 2e5, comb_gas);
%    test05 = refprop('T', [310;300], 'P', [2e5;1.8e5], R407C, ...
%        'Properties', 'P,T,h');
%    
%    EXAMPLES WITH OPTIONS:
%    comb_gas = [{'Oxygen';'Nitrogen';'CO2'},{0.19;0.80;0.01}, ...
%        {0.18;0.80;0.02},{0.17;0.80;0.03}];
%    test06 = refprop('T',[310;300;350],'P',[2e5;1.8e5;1.5e5], ...
%        comb_gas,'Uncertainties', 'yes', ...
%        'LimitsVal1', [309,311;299,301;349,351], ...
%        'LimitsVal2',[1.95e5,2.05e5;1.75e5,1.85e5;1.45e5,1.55e5]);
%

%% Options

% Default options
options = struct( ...
    'uncertainties', 'no', ...
    'limitsval1', [value1 value1], ...
    'limitsval2', [value2 value2], ...
    'properties', 'P,T,A,rho,s,h,x,mu,lambda,u,cp,cv', ...
    'quantityqualifier', 'nval', ...
    'lowervaluequalifier', 'lval', ...
    'uppervaluequalifier', 'uval', ...
    'absoluteerrorqualifier', 'aeval', ...
    'relativeerrorqualifier', 'reval');

% Read the acceptable names
optionNames = fieldnames(options);

% Count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('propertyName/propertyValue pairs needed.')
end

for pair = reshape(varargin,2,[]) % Pair is {propName;propValue}
   inpName = lower(pair{1}); % Make case insensitive

   if any(strcmp(inpName,optionNames))
      % Some overwriting of options are possible here
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
   
end

%% Some tests

if iscell(fluidstruct) ~= 1
    error('fluidstruct needs to be a cell object.')
end

% Are the value1 and value2 inputs column arrays?
are_col_arrays = sum([ ...
    size(value1,2)==1 && size(value1,1)>=1; ...
    size(value2,2)==1 && size(value2,1)>=1]) == 2;


if are_col_arrays ~= 1
    error('value1 and value2 need to be column arrays.')
end

% Is the first column of the fluidstruct a cell full of strings?
is_first_fluids_col_str = sum(cellfun(@(x) ischar(x), ...
    fluidstruct(:,1)))==size(fluidstruct(:,1),1);

if is_first_fluids_col_str ~= 1
    error(['The first column of fluidstruct needs to be filled ' ...
        'with strings.'])
end

% Is the columns 2:end of the fluidstruct are full of numerics?
is_num_next_cols_fluids = min(cellfun(@(x) isnumeric(x), ...
    fluidstruct(:, 2:end)));

if is_num_next_cols_fluids ~= 1
    error('Column 2:end of fluidstruct need to be filled numerics.')
end

% Are the vectors the same lengths?
are_lengths_egual_one = [ ...
    size(value1,1) == 1 ; ...
    size(value2,1) == 1 ; ...
    size(fluidstruct(:,2:end),2) == 1 ...
    ];

% Adapt the vectors sizes, when possible
if sum(are_lengths_egual_one == [0;0;1]) == 3
    fluidstruct(:, 2:size(value1, 1)+1) = repmat(fluidstruct(:, 2), 1, ...
        size(value1, 1));
    elseif sum(are_lengths_egual_one == [0;1;0]) == 3
    value2 = repmat(value2, size(value1, 1), 1);
elseif sum(are_lengths_egual_one == [1;0;0]) == 3
    value1 = repmat(value1, size(value2, 1), 1);
elseif sum(are_lengths_egual_one == [0;1;1]) == 3
    value2 = repmat(value2, size(value1, 1), 1);
    fluidstruct(:, 2:size(value1, 1)+1) = repmat(fluidstruct(:, 2), 1, ...
        size(value1, 1));
elseif sum(are_lengths_egual_one == [1;1;0]) == 3
    value1 = repmat(value1, size(fluidstruct(:, 2:end), 2), 1);
    value2 = repmat(value2, size(fluidstruct(:, 2:end), 2), 1);
elseif sum(are_lengths_egual_one == [1;0;1]) == 3
    value1 = repmat(value1, size(value2, 1), 1);
    fluidstruct(:, 2:size(value1, 1)+1) = repmat(fluidstruct(:, 2), 1, ...
        size(value2, 1));
end

are_lengths_the_same = sum([ ...
    size(value1, 1) == size(value2, 1) ;
    size(value1, 1) == size(fluidstruct(:, 2:end), 2);
    size(value2, 1) == size(fluidstruct(:, 2:end), 2) ...
    ]);

% If the vectors sizes are still not fine -> error
if are_lengths_the_same ~= 3
    error(['You have provided at least two inputs that have ' ...
        'incompatible sizes.']);
end

% If x is not included in the required properties and if uncertainties are
% required, we need to add it. Indeed, quality is need for the
% uncertainties calculations

if and(isempty(regexp(options.properties, 'x', 'once')), ...
        strcmp(options.uncertainties,'yes'))
    options.properties = [options.properties, ',x'];
    del_x = 1;
else
    del_x = 0;
end

% debug
% disp(value1)
% disp(value2)
% disp(fluidstruct)

%% Base function

    function results = compute_properties(spec1, value1, spec2, ...
            value2, fluidstruct)
        
        % Define the properties to be computed
        
        options.required_props = textscan(options.properties, '%s', ...
            'Delimiter', ',');
        options.required_props = unique(options.required_props{1});
        
        % Formatting of the output string used by refpropm
        % (c.f. requirements)
        
        options.outputs_str = cellfun(@(x)['results.' x '.' ...
            options.quantityqualifier '(i, 1), '], ...
            options.required_props, 'UniformOutput', 0);
        options.outputs_str = cell2mat(transpose(options.outputs_str));
        options.outputs_str = options.outputs_str(1:end-2);
        
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
        %        compliant with the current code base.
        
        options.refpropm_props = options.required_props;
        
        for i = 1:size(conv_table, 1)
            pos = cell2mat(cellfun(@(x)strcmp(x, conv_table{i, 1}), ...
                options.required_props, 'UniformOutput',0));
            if sum(pos) ~= 0
                options.refpropm_props{pos, 1} = conv_table{i, 2};
            end
            if strcmp(spec1, conv_table{i, 1})
                spec1 = conv_table{i, 2};
            end
            if strcmp(spec2, conv_table{i, 1})
                spec2 = conv_table{i, 2};
            end
        end
        
        options.refpropm_props = ...
            cell2mat(transpose(options.refpropm_props));
                
        % Fomatting of fluidstruct cell for refropm input
        fluids = cell2mat(cellfun(@(x) ['''' x ''', '], ...
            transpose(fluidstruct(:, 1)), 'UniformOutput', 0));
        
        % Units conversion for refpropm input
        if strcmp(spec1, 'P')
            value1 = value1 ./ 1e3;
        end
        
        if strcmp(spec2, 'P')
            value2 = value2 ./ 1e3;
        end
        
        % Initialization of the output structure
        results = struct;
        
        % Computations
        for i = 1:size(value1, 1)
            % Formatting of the composition vector
            composition = cell2mat(cellfun(@(x) [num2str(x) ', '], ...
                transpose(fluidstruct(:, i+1)), 'UniformOutput', 0));
            composition = ['[' composition(1:end-2) ']'];
            % generation of the refpropm command
            refpropm_cmd = ['[' options.outputs_str '] = refpropm(''' ...
                options.refpropm_props ''', ''' spec1 ''', ' ...
                'value1(i, 1)' ', ''' spec2 ''', ' 'value2(i, 1)' ', ' ...
                fluids composition ');'];
            % disp(refpropm_cmd); % debug
            % evaluation of the refpropm command generated above
            eval(refpropm_cmd); % the number of the outputs change
            % depending on the inputs, so an eval is unfortunately
            % necessary.
        end
        
        % Values adaptations
        
        % SI Units are Pa, not kPa
        if sum(ismember(options.refpropm_props, 'P'))==1
            results.P.(options.quantityqualifier) = ...
                results.P.(options.quantityqualifier) .* 1e3;
        end
        
        % 0<=x<=1
        if sum(ismember(options.refpropm_props, 'Q'))==1
            results.x.(options.quantityqualifier) ...
                (results.x.(options.quantityqualifier) < 0) = 0;
            results.x.(options.quantityqualifier) ...
                (results.x.(options.quantityqualifier) > 1) = 1;
        end
        
        % debug
        % disp(options.refpropm_props)
        % disp(options.outputs_str)
        % disp(refpropm_cmd)
        
    end

%% Computations

results = compute_properties(spec1, value1, spec2, value2, fluidstruct);

%% Uncertainties

if strcmp(options.uncertainties, 'yes')
    
    lvalue1 = options.limitsval1(:, 1);
    uvalue1 = options.limitsval1(:, 2);
    lvalue2 = options.limitsval2(:, 1);
    uvalue2 = options.limitsval2(:, 2);
    
    uncert01 = compute_properties(spec1, lvalue1, spec2, lvalue2, ...
        fluidstruct);
    uncert02 = compute_properties(spec1, uvalue1, spec2, uvalue2, ...
        fluidstruct);
    uncert03 = compute_properties(spec1, lvalue1, spec2, uvalue2, ...
        fluidstruct);
    uncert04 = compute_properties(spec1, uvalue1, spec2, lvalue2, ...
        fluidstruct);
    uncert01_sat = compute_properties(...
        'P',uncert01.P.(options.quantityqualifier), ...
        'Q',uncert01.x.(options.quantityqualifier), ...
        fluidstruct);
    uncert02_sat = compute_properties(...
        'P',uncert02.P.(options.quantityqualifier), ...
        'Q',uncert02.x.(options.quantityqualifier), ...
        fluidstruct);
    uncert03_sat = compute_properties(...
        'P',uncert03.P.(options.quantityqualifier), ...
        'Q',uncert03.x.(options.quantityqualifier), ...
        fluidstruct);
    uncert04_sat = compute_properties(...
        'P',uncert04.P.(options.quantityqualifier), ...
        'Q',uncert04.x.(options.quantityqualifier), ...
        fluidstruct);
    
    properties = fieldnames(results);
    for j = 1:size(properties, 1)
        uncert_mat = [ ...
            uncert01.(properties{j, 1}).(options.quantityqualifier), ...
            uncert02.(properties{j, 1}).(options.quantityqualifier), ...
            uncert03.(properties{j, 1}).(options.quantityqualifier), ...
            uncert04.(properties{j, 1}).(options.quantityqualifier)];
        uncert_sat_mat = [ ...
            uncert01_sat.(properties{j, 1}).(options.quantityqualifier), ...
            uncert02_sat.(properties{j, 1}).(options.quantityqualifier), ...
            uncert03_sat.(properties{j, 1}).(options.quantityqualifier), ...
            uncert04_sat.(properties{j, 1}).(options.quantityqualifier)];
        abs_err = 2*max(abs(uncert_mat - ...
            repmat(results.(properties{j, 1}).(options.quantityqualifier), ...
            1, 4)), [], 2);
        abs_err_sat = 2*max(abs(uncert_sat_mat - ...
            repmat(results.(properties{j, 1}).(options.quantityqualifier), ...
            1, 4)), [], 2);
        % the goal of this is to not cross over the saturation line, if we
        % are in the gas or the liquid state
        if ~isempty(find(or(results.x.(options.quantityqualifier)==0, ...
                results.x.(options.quantityqualifier)==1), 1))
            try
                abs_err(or(results.x.(options.quantityqualifier) == 0, ...
                        results.x.(options.quantityqualifier) == 1)) = ...
                        min([abs_err, abs_err_sat], [], 2);
            catch err
                disp(err)
            end
        end
        results.(properties{j, 1}).(options.absoluteerrorqualifier) = abs_err;
        results.(properties{j, 1}).(options.relativeerrorqualifier) = ...
            abs_err ./ results.(properties{j, 1}).(options.quantityqualifier);
        results.(properties{j, 1}).(options.lowervaluequalifier) = ...
            results.(properties{j, 1}).(options.quantityqualifier) - ...
            results.(properties{j, 1}).(options.absoluteerrorqualifier)/2;
        results.(properties{j, 1}).(options.uppervaluequalifier) = ...
            results.(properties{j, 1}).(options.quantityqualifier) + ...
            results.(properties{j, 1}).(options.absoluteerrorqualifier)/2;       
    end
end

if del_x == 1
    results = rmfield(results, 'x');
else
    results.x.(options.absoluteerrorqualifier)(results.x.(options.quantityqualifier)==1) = NaN;
    results.x.(options.relativeerrorqualifier)(results.x.(options.quantityqualifier)==1) = NaN;
    results.x.(options.lowervaluequalifier)(results.x.(options.quantityqualifier)==1) = NaN;
    results.x.(options.uppervaluequalifier)(results.x.(options.quantityqualifier)==1) = NaN;
    results.x.(options.absoluteerrorqualifier)(results.x.(options.quantityqualifier)==0) = NaN;
    results.x.(options.relativeerrorqualifier)(results.x.(options.quantityqualifier)==0) = NaN;
    results.x.(options.lowervaluequalifier)(results.x.(options.quantityqualifier)==0) = NaN;
    results.x.(options.uppervaluequalifier)(results.x.(options.quantityqualifier)==0) = NaN;
end

end
