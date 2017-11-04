function code = generateSFunction(varargin)
%% generateSFunction Generate S-Function for simulink
% 
% code = generateSFunction( nx,nu,nw,N, <options>, <value> )
%
% Returns and compile the .c file containing the S-Function and create a Simulink
% block interfacing the algorithm
%
% required inputs:
%           .nu                         - size of input vector
%           .nx                         - size of state vector
%           .nw                         - size of known disturbance vector
%           .N                          - horizon length
% options: 
%           name                       - .c file name. Default: 'simulink_my_code'
%           function_name              - name of the function called by the S-Function
%                                        Default: 'proposed_algorithm'
%           model_name                 - name of the generated Simulink block
%                                        Default: 'FalcOpt_lib'
%           maskImage_name             - name of the simulink mask image
%           maxIt                      - max number of iterations. Default: 0
%           trackReference             - if 'true' add input for reference
%           terminal                   - if 'true' add input for 'c' terminal constant
%           contractive                - if 'true' add input for 'c' contractive constant
%           indent                     - indentation to be used in code generation
%           precision                  - 'single' or 'double'
%           type                       - data type, either 'float', 'single' or 'double',
%           debug                      - level of debug code that is generated. Default is 1.
%
% Copyright (c) 2017 Tommaso Robbiani <falcopt@tommasorobbiani.com>
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
defaultNames = struct('maximumIterations', 'maximumIterations','fun', 'myfun');

p = inputParser;
p.addRequired('nx', @isnumeric);
p.addRequired('nu', @isnumeric);
p.addRequired('nw', @isnumeric);
p.addRequired('N', @isnumeric);
p.addParameter('maxIt', 0, @(x)(isnumeric(x) && x>=0 && mod(x,1)==0));
p.addParameter('name','simulink_my_code', @ischar);
p.addParameter('names', defaultNames, @isstruct);
p.addParameter('model_name','FalcOpt_lib', @ischar);
p.addParameter('maskImage_name','mask_image',@ischar);
p.addParameter('trackReference',false, @islogical);
p.addParameter('terminal',false, @islogical);
p.addParameter('contractive',false, @islogical);
p.addParameter('indent', '\t', @(x)(ischar(x)));
p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double'}))));
p.addParameter('type', 'double', @(x)(ischar(x) && any(strcmp(x, {'float', 'single', 'double'}))));
p.addParameter('debug', 1, @(x)(isnumeric(x) && x >= 0));
p.addParameter('extra_output', false, @islogical);
p.parse(varargin{:});
o = p.Results;

code = [];
s_name = o.name; % name of .c file
function_name = o.names.fun;   % name of function called inside wrapper
model_name = o.model_name;

% Names
fields = fieldnames(defaultNames);
for i=1:length(fields)
    if ~isfield(o.names, fields{i})
        o.names.(fields{i}) = defaultNames.(fields{i});
    end
end

% Type
if strcmp(o.type, 'single')
    warning('falcopt:generateMEX:InvalidType', 'The type "single" is depricated and will not be allowed in a future version.');
    o.type = 'float';
end
% Precision (due to legacy usage)
if isempty(o.precision)
    warning('falcopt:generateMEX:MissingPrecision', 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
    switch o.type
        case 'float'
            o.precision = 'single';
        case 'single'
            o.precision = 'single';
        case 'double'
            o.precision = 'double';
        otherwise
            warning('falcopt:generateMEX:InvalidPrecision', ['Precision was not specified, could not infer proper precision from given type "' o.type '". Setting to "double".']);
            o.precision = 'double';
    end
end
switch o.type
    case 'double'
        real_T = 'double';
    case 'float'
        real_T = 'float';
end
switch o.precision
    case 'double'
        real_s = 'SS_DOUBLE';
        real_c = 'real_T';
    case 'single'
        real_s = 'SS_SINGLE';
        real_c = 'real_T';
end


%% display informations
fprintf('\n\ngenerating Simulink S-Function and block ...')

%% parameter settings

code = [ code, sprintf(['#define S_FUNCTION_NAME ' s_name '\n'...
        '#define S_FUNCTION_LEVEL 2\n\n'])];
    
% input ports parameters
% for each port must specify: 
%   .name: the name that will appear in simlink block
%   .size: the size of the port 
%   .data_s: data type in simulink
%   .data_c: data type will be passed to algorithm (c-type)

% u_pred
in.u_pred.name = 'u_pred';
in.u_pred.size = o.nu*o.N;
in.u_pred.data_s = real_s;
in.u_pred.data_c = real_c;

if o.nw > 0
    % noise w
    in.w.name = 'w';
    in.w.size = o.nw;
    in.w.data_s = real_s;
    in.w.data_c = real_c;
end
if o.trackReference
    % reference x
    in.xref.name = 'xref';
    in.xref.size = o.nx*o.N;
    in.xref.data_s = real_s;
    in.xref.data_c = real_c;
    
    % reference x
    in.uref.name = 'uref';
    in.uref.size = o.nu*o.N;
    in.uref.data_s = real_s;
    in.uref.data_c = real_c;
end
if o.terminal||o.contractive
    if o.contractive
        in.c_ind.name = 'contraction_idx';
        in.c_ind.size = 1;
        in.c_ind.data_s = 'SS_UINT32';
        in.c_ind.data_c = 'uint32_T';
    end
    in.c.name = 'contraction_c';
    in.c.size = 1;
    in.c.data_s = real_s;
    in.c.data_c = real_c;
end
% x0
in.x0.name = 'x0';
in.x0.size = o.nx;
in.x0.data_s = real_s;
in.x0.data_c = real_c;

in_names = fieldnames(in);
nInputs = length(in_names); %nuber of inputs
code = [code, sprintf('#define NUM_INPUTS %i\n',nInputs)];

% output ports parameters
out.u_opt.name = 'u_opt';
out.u_opt.size = o.nu;
out.u_opt.data_s = real_s;
out.u_opt.data_c = real_c;

out.u.name = 'u_out';
out.u.size = o.nu*o.N;
out.u.data_s = real_s;
out.u.data_c = real_c;


if o.debug > 1
    % fval
    out.fval.name = 'fval';
    out.fval.size = 1;
    out.fval.data_s = real_s;
    out.fval.data_c = real_c;
    
    % x
    out.x.name = 'x';
    out.x.size = o.nu*o.N;
    out.x.data_s = real_s;
    out.x.data_c = real_c;   
end

if o.debug > 2
    
    % optimval
    out.optimval.name = 'optimval';
    out.optimval.size = 1;
    out.optimval.data_s = real_s;
    out.optimval.data_c = real_c;
    
    % feasval
    out.feasval.name = 'feasval';
    out.feasval.size = 1;
    out.feasval.data_s = real_s;
    out.feasval.data_c = real_c;
    
    % meritval
    out.meritval.name = 'meritval';
    out.meritval.size = 1;
    out.meritval.data_s = real_s;
    out.meritval.data_c = real_c;
end

% number iteration
out.iter.name = 'iter';
out.iter.size = 1;
out.iter.data_s = 'SS_UINT32';
out.iter.data_c = 'uint32_T';

% exitflag
out.flag.name = 'flag';
out.flag.size = 1;
out.flag.data_s = 'SS_INT32';
out.flag.data_c = 'int32_T';

out_names = fieldnames(out);
nOutputs = length(out_names); % number outputs
code = [code, sprintf('#define NUM_OUTPUTS %i\n',nOutputs)];

%% include     
code = [ code, sprintf('#include "simstruc.h"\n\n')];

code = [ code, sprintf(['#if defined(MATLAB_MEX_FILE)\n'...
                        '#include "tmwtypes.h"\n'...
                        '#include "simstruc_types.h"\n'...
                        '#else\n'...
                        '#include "rtwtypes.h"\n'...
                        '#endif\n\n'])];                   
 
                    
%% algorithm wrapper function
code = [code, sprintf('void algorithm_wrapper(')];
for k=1:nInputs
    code = [code, sprintf([' const %s *%s,'],in.(in_names{k}).data_c,in.(in_names{k}).name)]; %#ok
end
for k=1:nOutputs-1
    code = [code, sprintf([' %s *%s,'],out.(out_names{k}).data_c,out.(out_names{k}).name)]; %#ok
end
code = [code, sprintf(' %s *%s){\n\n',out.(out_names{end}).data_c,out.(out_names{end}).name)];

% internal variables ( max number iteration and vector u )
code = [code, sprintf([o.indent 'int i=0;\n' ...
                       o.indent 'unsigned int lineSearch_nit[' falcopt.internal.toDefineName(o.names.maximumIterations) '];\n' ...
                       o.indent real_T ' u[%i];\n'], o.nu*o.N)];

% Prepare inputs
code = [code, sprintf(['\n' o.indent 'for( i=0;i<%i;i++){ u[i] = %s[i];}\n'],o.nu*o.N,in.u_pred.name)];

% call algorithm function
code = [code, sprintf([o.indent '%s[0] = %s( %s, u, '],out.flag.name,function_name,in.x0.name)]; % flag = fct_name(x0,u,
if o.nw>0
    code = [code, sprintf('%s, ',in.w.name)]; % w
end
if o.trackReference
    code = [code, sprintf('%s, %s, ',in.xref.name,in.uref.name)];  % xref, uref
end
if o.terminal||o.contractive
    if o.contractive
        code = [code, sprintf('*%s, ',in.c_ind.name)]; % contraction_idx
    end
    code = [code, sprintf('*%s, ',in.c.name)]; % contraction_c
end

if o.debug > 1
    code = [code, sprintf('%s, %s, ', out.x.name, out.fval.name)]; % x, fval
end
if o.debug > 2
    code = [code, 'iter, lineSearch_nit, ']; %TODO tommaso check
    code = [code, sprintf('%s, %s, %s);\n', out.optimval.name, out.feasval.name, out.meritval.name)]; % optimval, feasval, meritval
else
    code = [code, sprintf('iter, lineSearch_nit);\n ')];%TODO tommaso check
end

% prepare outputs
code = [code, sprintf([o.indent 'for( i=0;i<%i;i++){ %s[i] = u[i];}\n'],o.nu*o.N,out.u.name)];
code = [code, sprintf([o.indent 'for( i=0;i<%i;i++){ %s[i] = u[i];}\n'],o.nu,out.u_opt.name)];

code = [code, sprintf('}\n\n')];



%% SIMULINK S-FUNCTION

code = [code, sprintf(['/*====================*\n'...
                       ' * S-function methods *\n'...
                       ' *====================*/\n\n'])];
code = [ code, sprintf('static void mdlInitializeSizes(SimStruct *S){\n\n')];

code = [code, sprintf([o.indent 'DECL_AND_INIT_DIMSINFO(inputDimsInfo);\n\tDECL_AND_INIT_DIMSINFO(outputDimsInfo);\n\n'])];

code = [code, sprintf([o.indent 'ssSetNumSFcnParams(S, 0);\n' ...
                o.indent 'if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {\n' ...
                o.indent o.indent 'return; /* Parameter mismatch reported by the Simulink engine*/\n\t}\n\n'])];

%% initialize inputs port
code = [code, sprintf([o.indent '/* total number of inputs port: %i */\n'...
                o.indent 'if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;\n\n'],nInputs)];
for k=1:nInputs
    code = [code, sprintf(['\n' o.indent '/* Input port %i (%s) */\n'...
                o.indent 'ssSetInputPortWidth(S,  %i, %i);\n'...
                o.indent 'ssSetInputPortDirectFeedThrough(S, %i, 1);\n'...
                o.indent 'ssSetInputPortComplexSignal(S, %i, COMPLEX_NO);\n'...
                o.indent 'ssSetInputPortRequiredContiguous(S, %i, 1);\n'...
                o.indent 'ssSetInputPortDataType(S, %i, %s);\n'],k-1,in.(in_names{k}).name,k-1,in.(in_names{k}).size,k-1,k-1,k-1,k-1,in.(in_names{k}).data_s)]; %#ok
end


% outputs port
code = [code, sprintf(['\n' o.indent '/* total number of outputs port: %i */\n'...
                o.indent 'if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;\n\n'],nOutputs)];
            
for k=1:nOutputs
    code = [code, sprintf(['\n' o.indent '/* Output port %i (%s) */\n'...
        o.indent 'ssSetOutputPortWidth(S, %i, %i);\n'...
        o.indent 'ssSetOutputPortComplexSignal(S, %i, COMPLEX_NO);\n'...
        o.indent 'ssSetOutputPortDataType(S, %i, %s);\n'],k-1,out.(out_names{k}).name,k-1,out.(out_names{k}).size,k-1,k-1,out.(out_names{k}).data_s)]; %#ok
end

% initialize sampling time
code = [code, sprintf(['\n' o.indent '/* initialize sampling time */\n'...
                o.indent 'ssSetNumSampleTimes(S, 1);\n'],nOutputs,nOutputs)];
            
code = [code, sprintf(['\n' o.indent '/* Take care when specifying exception free code - see sfuntmpl.doc */\n'...
                        o.indent 'ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE | SS_OPTION_WORKS_WITH_CODE_REUSE));\n'])];
                    
code = [code, sprintf([o.indent 'ssSetNumRWork(S, 0);\n'...
                o.indent 'ssSetNumIWork(S, 0);\n'...
                o.indent 'ssSetNumPWork(S, 0);\n'...
                o.indent 'ssSetNumModes(S, 0);\n'...
                o.indent 'ssSetNumNonsampledZCs(S, 0);\n}\n\n'])];                    
 %% Sampling time
 code = [code, sprintf('static void mdlInitializeSampleTimes(SimStruct *S){\n')];
 code = [code, sprintf([o.indent 'ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);\n'...
                        o.indent 'ssSetModelReferenceSampleTimeDefaultInheritance(S);\n'...
                        o.indent 'ssSetOffsetTime(S, 0, 0.0);\n}\n\n'])];
                    
%% Output function
code = [code, sprintf('static void mdlOutputs(SimStruct *S, int_T tid){\n\n')];

% create pointer to input signals
code = [code, sprintf([o.indent '/* Pointers to input ports */\n'])];
for k=1:nInputs
    code = [code, sprintf([o.indent 'const real_T *%s = (const real_T*) ssGetInputPortSignal(S,%i);\n'],in.(in_names{k}).name,k-1)]; %#ok
end

% create pointer to output singnals
code = [code, sprintf(['\n' o.indent '/* Pointers to output ports */\n'])];
for k=1:nOutputs
    code = [code, sprintf([o.indent '%s *%s = (%s*) ssGetOutputPortRealSignal(S,%i);\n'],out.(out_names{k}).data_c,out.(out_names{k}).name,out.(out_names{k}).data_c,k-1)]; %#ok
end

% call algorithm function
code = [code, sprintf(['\n' o.indent '/* call algorithm wrapper */\n'])];
code = [code, sprintf([o.indent 'algorithm_wrapper('])];
for k=1:nInputs
    code = [code, sprintf([' %s,'],in.(in_names{k}).name)]; %#ok
end
for k=1:nOutputs-1
    code = [code, sprintf([' %s,'],out.(out_names{k}).name)]; %#ok
end
code = [code, sprintf(' %s);\n',out.(out_names{end}).name)];

code = [code, sprintf('}\n\n')];

code = [code, sprintf('static void mdlTerminate(SimStruct *S){}\n\n')];

code = [code, sprintf(['#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */\n'...
                       '#include "simulink.c" /* MEX-file interface mechanism */\n'...
                       '#else\n'...
                       '#include "cg_sfun.h" /* Code generation registration function */\n'...
                       '#endif\n'])];

%% generate simulink model 
block_name = sprintf([model_name '/FalcOpt_solver']);
% open/create model
try 
    open_system( model_name,'loadonly');
catch
    new_system( model_name,'Library');
end

% delete existing block
try  %#ok
    set_param( model_name,'Lock','off');
    delete_block( block_name);
end

% create new block
add_block('simulink/User-Defined Functions/S-Function', block_name,'MakeNameUnique', 'on');
set_param( block_name,'FunctionName',s_name);
set_param( block_name,'Position',[50 50 400 300])
set_param( block_name,'Mask','on');

% create Mask
mask = [];
% image Mask
try
    ud.logo = imread(o.maskImage_name);
    [h_im,w_im,~] = size(ud.logo);
    set_param(gcb,'position',[5,5,w_im,h_im]);
    set_param(gcb,'UserDataPersistent','on');
    set_param(gcb,'UserData',ud);
    mask = [mask, sprintf('ud = get_param(gcb,''UserData'');\n')];
    mask = [mask, sprintf('image(ud.logo,''center'');\n')];
catch
    warning(['image ' o.maskImage_name ' not found. Ignoring simulink block logo']);
    mask = [mask,'disp(''FalcOpt\n %% no logo %%'');',sprintf('\n')];
end

% input ports label
for k=1:nInputs
    mask = [mask, sprintf(['port_label(''input'',%i,''%s'')\n'],k,in.(in_names{k}).name)]; %#ok
end
% outputs port label
for k=1:nOutputs
    mask = [mask, sprintf(['port_label(''output'',%i,''%s'')\n'],k,out.(out_names{k}).name)]; %#ok
end

set_param( block_name,'MaskDisplay', mask);

set_param( model_name,'Lock','on');
save_system(sprintf([model_name '.mdl']));
close_system( model_name);

end