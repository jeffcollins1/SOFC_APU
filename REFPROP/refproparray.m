function varargout = refproparray( varargin )
%refproparray Interface for refpropm which allows arrays to be entered
%instead of individual data points
%
%Function calls are very similar to refpropm (identical if only a single
%value is to be evaluated) except that vectors can be entered
%
%To use custom mixtures, the substances and concentrations must be arranged
%in rows vectors
%Lengths of the values, substances, concentrations must be the same length
%or of length 1.
%
%Programmed by Scott Wujek, Creative Thermal Solutions, Fall 2013
%Updated 3/15/14 to allow for an array of outputs to 'x' function call and
%to allow for '<' and '>' forcing

%% convert known numeric input variables to numbers
propVal1 = cell2mat(varargin(3));
propVal2 = cell2mat(varargin(5));

%% determine vector lengths
L1=length(propVal1);
L2=length(propVal2);
L3=size(char(varargin(6)),1);
if nargin>6
    concentrations=cell2mat(varargin(nargin));
    L4=size(concentrations,1);
else
    L4=NaN;
    concentrations=[];
end
maxsize=nanmax([L1, L2, L3, L4]);

%%  Make vectors to all be same length (expand singular entries)
if maxsize>1  && (L1==maxsize || L1==1) && (L2==maxsize || L2==1) && (L3==maxsize || L3==1) && (L4==maxsize || L4==1 || isnan(L4))
    if L1==1
        propVal1=repmat(propVal1, [maxsize,1]); %#ok<NASGU>
    end
    if L2==1
        propVal2=repmat(propVal2, [maxsize, 1]); %#ok<NASGU>
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
        concentrations=repmat(concentrations, [maxsize,1]); %#ok<NASGU>
    end
end