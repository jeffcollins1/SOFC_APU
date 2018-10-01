function varargout = refpropm2D(PROPERTY,PARAMETER1,VALUE1,PARAMETER2,VALUE2,varargin)
% propertyMatrix = refpropm2D('PROPERTY','PARAMETER1',VALUE1,'PARAMETER2',VALUE2,'SUBSTANCE')
%
% This function returns a matrix of thermophysical properties as a function
% of two independent state variables. Properties are determined by calling
% the NIST REFPROP database. The property returned by the function is
% defined by the input PROPERTY. The independent state variable types are
% defined by inputs PARAMETER1 and PARAMETER2.
% Type ">> help refpropm" for the list of variables and properties.
%
% The vector values for the two independent state variables are defined by
% inputs VALUE1 and VALUE2. If VALUE1 and VALUE2 are of length M and N
% respectively, then propertyMatrix will be M by N.
%
% For a pure substance, input SUBSTANCE should be a string with the
% corresponding fluid name. For a mixture of substances, SUBSTANCE should
% be a cell array of the substance names plus a final vector element
% representing the corresponding mass fractions. A maximum of 20 substances
% can be handled. Valid substance names are given by the file names with
% extension "*.FLD" in the ~\fluids\ directory.
%
% Example 1: Specific heat at constant pressure for a pure methane
%   cpMatrix = refpropm2D('C','T',[200,290,300],'P',[100,200],'methane');
%
% Example 2: Density for mixture of CF3CH2F and CH2F2
%   substances = {'R134a','R32', [0.8,0.2]};
%   rhoMatrix = refpropm2D('D','T', [280,290,300],'P',[100,200],substances);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Programmed by Xiaochuan Yuan and Steve Miller, MathWorks, Dec 2013


% Pre-assign property matrix for efficiency
M = length(VALUE1);
N = length(VALUE2);
propertyMatrix = zeros(M, N);

% Determine if substance is pure or not
if (nargin==6)
    if iscell(varargin{1})
        substanceList = varargin{1};
    else
        substanceList = {varargin{1}};
    end
else
    % ASSEMBLE SUBSTANCE FROM INPUTS
    substanceList = {varargin{1:end-1},varargin{end}};
end

% Calculate property values
for i = 1:length(PROPERTY)
    for iRow = 1:M
        for jCol = 1:N
            varargout{i}(iRow, jCol) = refpropm(PROPERTY(i), PARAMETER1, ...
                VALUE1(iRow), PARAMETER2, VALUE2(jCol), substanceList{:});
        end
    end
end
end