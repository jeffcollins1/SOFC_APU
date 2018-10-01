function varargout = refproparray( varargin ) 
%refproparray Interface for refpropm which allows arrays to be entered %instead of individual data points % 
%Function calls are very similar to refpropm (identical if only a single %value is to be evaluated) except that vectors can be entered % 
%To use custom mixtures, the substances and concentrations must be arranged %in rows vectors %Lengths of the values, substances, concentrations must be the same length %or of length 1. % 
%Programmed by Scott Wujek, Creative Thermal Solutions, Fall 2013 
%Updated 3/15/14 to allow for an array of outputs to 'x' function call and %to allow for '<' and '>' forcing %
% convert known numeric input variables to numbers 
propVal1 = cell2mat(varargin(3)); 
propVal2 = cell2mat(varargin(5)); %
% determine vector lengths 
L1=length(propVal1); 
L2=length(propVal2); 
L3=size(char(varargin(6)),1); 
if nargin>6 
    concentrations=cell2mat(varargin(nargin)); 
    L4=size(concentrations,1); 
else L4=NaN; concentrations=[]; 
end
maxsize=nanmax([L1, L2, L3, L4]);
%% Make vectors to all be same length (expand singular entries) 
if maxsize>1 && (L1==maxsize || L1==1) && (L2==maxsize || L2==1) && (L3==maxsize || L3==1) && (L4==maxsize || L4==1 || isnan(L4)) 
    if L1==1 
        propVal1=repmat(propVal1, [maxsize,1]); %#ok 
    end 
    if L2==1 
        propVal2=repmat(propVal2, [maxsize, 1]); %#ok 
    end 
    if L3==1 
        if nargin==6 
            varargin{6}=repmat(varargin(6), [maxsize,1]); 
        else
            for inarg=6:nargin-1 
                varargin{inarg}=repmat(varargin(inarg), [maxsize,1]); 
            end
        end
    end
    if L4==1 
        concentrations=repmat(concentrations, [maxsize,1]); %#ok 
    end
elseif maxsize~=1 
    error('Input value vectors are of different lengths or not of length 1') 
end
%% determine number of properties requested, make output variable, and left side of equation that replicates format of refpropm 
propReq = lower(char(varargin(1))); 
numoutputvars=sum(propReq~='>'&propReq~='<'); 
if propReq=='x' 
    varargout = cell([1,2]); 
    variables=genvarname(repmat({'temp'}, [1,2])); 
    for i=1:length(variables) 
        eval([char(variables(i)) '=NaN(maxsize,2);']); 
    end
    string1='[temp(:,arraynum), temp1(:,arraynum)]'; 
else
    varargout = cell(1, numoutputvars);
    variables=genvarname(repmat({'temp'}, [1,numoutputvars])); 
    string1='['; 
    for i=1:numoutputvars 
        eval([char(variables(i)) '=NaN(size(propVal1));']); 
        string1=[string1, char(variables(i))]; %#ok 
        if i~=length(variables) 
            string1=[string1, '(arraynum),']; %#ok 
        else
            string1=[string1, '(arraynum)]']; %#ok
        end
    end
end
if propReq=='x' 
    eval('temp=NaN(size(concentrations,2), maxsize);') 
    eval('temp1=NaN(size(concentrations,2), maxsize);') 
end
%% create right side of refprop function call, and evaluate function 
for arraynum=1:maxsize 
    if nargin==6 %this is case for pure substances or .mix files 
        string2='=refpropm(char(varargin(1)),char(varargin(2)),propVal1(arraynum),char(varargin(4)),propVal2(arraynum),char(varargin{6}(arraynum,:)));'; 
    else %this is the case for custom mixtures 
        string2='=refpropm(char(varargin(1)),char(varargin(2)),propVal1(arraynum),char(varargin(4)),propVal2(arraynum),'; 
        for inarg=6:nargin-1 
            string2=strcat(string2, 'char(varargin{', num2str(inarg), '}(arraynum,:)), '); 
        end
        string2=strcat(string2, 'concentrations(arraynum,:));');
    end
    eval([string1, string2])
end
%% fill in output vector(s) 
if propReq=='x' 
    eval(['varargout(1) = {' char(variables(1)) '''};']); 
    eval(['varargout(2) = {' char(variables(2)) '''};']); 
else
    for varnum=1:length(varargout) 
        eval(['varargout(varnum,1) = {' char(variables(varnum)) '};']); 
    end
end