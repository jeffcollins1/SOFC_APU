%% Main function to generate tests
function tests = refpropm_basic_test
tests = functiontests(localfunctions);
end

%% Test Functions

function example_codes_num1_test(testCase)
% #1 refpropm example code
actSolution_case1 = refpropm('P','T',373.15,'Q',0,'water');
expSolution_case1 = 1.014179966562579e+02;
% Check equality
verifyEqual(testCase, actSolution_case1, expSolution_case1, ...
    'AbsTol', 1e-10);
end

function example_code_num2_test(testCase)
% #2 refpropm example code
[actSolution_case2_a, actSolution_case2_b]  = ...
    refpropm('SC','T',373.15,'Q',1,'water');
expSolution_case2_a = 7.354119145712214e+03;
expSolution_case2_b = 2.080041256512282e+03;
% Check equality
verifyEqual(testCase, actSolution_case2_a, expSolution_case2_a, ...
    'AbsTol', 1e-10);
verifyEqual(testCase, actSolution_case2_b, expSolution_case2_b, ...
    'AbsTol', 1e-10);
end

function example_codes_num3_test(testCase)
% #3 refpropm example code
actSolution_case3 = refpropm('D','T',323.15,'P',1e2,'water','ammonia',[0.9 0.1]);
expSolution_case3 = 9.463790683846696e+02;
% Check equality
verifyEqual(testCase, actSolution_case3, expSolution_case3, ...
    'AbsTol', 1e-10);
end

function example_codes_num4_test(testCase)
% #4 refpropm example code
[actSolution_case4_a, actSolution_case4_b]  = ...
    refpropm('X','P',5e2,'Q',0.4,'R134a','R32',[0.8, 0.2]);
expSolution_case4_a = [0.857408717090157;0.142591282909843];
expSolution_case4_b = [0.713889524601922;0.286110475398078];
% Check equality
verifyEqual(testCase, actSolution_case4_a, expSolution_case4_a, ...
    'AbsTol', 1e-10);
verifyEqual(testCase, actSolution_case4_b, expSolution_case4_b, ...
    'AbsTol', 1e-10);
end

function example_codes_num5_test(testCase)
% #5 refpropm example code
actSolution_case5 = refpropm('T','C',0,' ',0,'water');
expSolution_case5 = 6.470960000000000e+02;
% Check equality
verifyEqual(testCase, actSolution_case5 ,expSolution_case5, ...
    'AbsTol', 1e-10);
end

function example_codes_num6_test(testCase)
% #6 refpropm example code
actSolution_case6 = refpropm('T','M',0,' ',0,'r410a.mix');
expSolution_case6 = 4.602980582532530e+02;
% Check equality
verifyEqual(testCase, actSolution_case6, expSolution_case6, ...
    'AbsTol', 1e-10);
end

%% File fixtures  
function setupOnce(testCase)  % do not change function name
    refpropm_function_path = which('refpropm');
    if isempty(refpropm_function_path)
        error('refpropm function is not loaded in your MATLAB path. Abording...')
    else
        fprintf('refpropm function being used: %s\n', refpropm_function_path);
    end
end

function teardownOnce(testCase)  % do not change function name
end

%% Fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end