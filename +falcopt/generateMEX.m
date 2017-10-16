%% generateMEX Generate MEX function
% 
% [code, info] = generateMEX(N, dims)
% or
% [code, info] = generateMEX(N, dims, ...) [with options]
% where 'dims' is a struct with the fields .x, .u and .w
%  the dimension of the disturbance 'w' is optional and set to zero if omitted
%
% Returns a string 'code' containing code for the MEX function
%
% the following options are available:
%  .names    - The name of the MEX function
%  .ref     - A boolean, if true then x and u references are expected
%  .precision - A string that determines the precision of computation and stored data in the generated code. 
%                Needs to be either 'single' or 'double'. Default: '' (unspecified, will try to extract from type);
%  .type    - The data type, either 'float', 'single' or 'double',
%  .indent  - Indentation to be used in code generation
%  .inline  - Inline keyword to be uised. Default is inline
%  .verbose - Level of procedural output of this function. Default is 0.
%  .debug   - Level of debug code that is generated. Default is 1.
% If fields of options structs are ommitted, the default is used.
% 

% Copyright (c) 2017 Damian Frick <falcopt@damianfrick.com>
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
function code = generateMEX(varargin)
    indentTypes = {'generic', 'code'};
    defaultNames = struct('mex', 'mymex', 'fun', 'myfun');
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('N', @(x)(isnumeric(x) && x > 0 && mod(x,1) == 0));
    p.addRequired('dims', @(x)(isstruct(x) && all(isfield(x, {'x','u'}))) ...
                                           && (isnumeric(x.x) && x.x > 0 && mod(x.x,1) == 0) ...
                                           && (isnumeric(x.u) && x.u > 0 && mod(x.u,1) == 0) ...
                                           && (~isfield(x, 'w') || (isnumeric(x.w) && x.w >= 0 && mod(x.w,1) == 0)));
    p.addParameter('names', defaultNames, @isstruct);
    p.addParameter('ref', false, @islogical);
    p.addParameter('terminalContraction', false, @islogical);
    p.addParameter('indexContraction', false, @islogical);
    p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double'}))));
    p.addParameter('type', 'double', @(x)(ischar(x) && any(strcmp(x, {'float', 'single', 'double'}))));
    p.addParameter('maxIt', 0, @(x)(isnumeric(x) && x>=0 && mod(x,1)==0));
    p.addParameter('indent', '    ', @(x)(ischar(x) || (isstruct(x) && isfield(x, 'generic') && all(cellfun(@(y)(~isfield(x, y) || ischar(x.(y))), indentTypes)))));
    p.addParameter('inline', 'inline', @ischar);
    p.addParameter('timing', false, @islogical);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('debug', 1, @(x)(isnumeric(x) && x >= 0));
    p.parse(varargin{:});
    options = p.Results;
    N = options.N;
    dims = options.dims;
    
    if options.verbose >= 2
        fprintf('MEX code generation started...\n');
        fprintf('. processing options\n');
    end
    
    %% Processing options
    % Dimensions
    if ~isfield(dims, 'w')
        dims.w = 0;
    end
    % Names
    fields = fieldnames(defaultNames);
    for i=1:length(fields)
        if ~isfield(options.names, fields{i})
            options.names.(fields{i}) = defaultNames.(fields{i});
        end
    end
    % Type
     if strcmp(options.type, 'single')
        warning('falcopt:generateMEX:InvalidType', 'The type "single" is depricated and will not be allowed in a future version.');
        options.type = 'float';
    end
    % Precision (due to legacy usage)
    if isempty(options.precision)
        warning('falcopt:generateMEX:MissingPrecision', 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
        switch options.type
            case 'float'
                options.precision = 'single';
            case 'single'
                options.precision = 'single';
            case 'double'
                options.precision = 'double';
            otherwise
                warning('falcopt:generateMEX:InvalidPrecision', ['Precision was not specified, could not infer proper precision from given type "' options.type '". Setting to "double".']);
                options.precision = 'double';
        end
    end
    switch options.type
        case 'double'
            dataType = 'double';
        case 'float'
            dataType = 'float';
    end
    % Indentation
    if ~isstruct(options.indent)
        indent = options.indent;
        options.indent = struct();
        for i=1:length(indentTypes)
            options.indent.(indentTypes{i}) = indent;
        end
    end
    for i=1:length(indentTypes)
        if ~isfield(options.indent, indentTypes{i})
            options.indent.(indentTypes{i}) = options.indent.generic;
        end
        options.indent.(indentTypes{i}) = options.indent.(indentTypes{i})(:)'; % Make sure is row vector
    end
    
    %% Generate code
    if options.verbose == 1
        fprintf(['Generating code for ' options.names.mex '()\n']);
    end
    
    code = '';
    if options.timing
        if options.verbose >= 2
            fprintf('. generating timer\n');
        end
        
        code = [code, sprintf(['\n' ...
                           options.indent.code '/*********' '\n' ...
                           options.indent.code ' * Timer *' '\n' ...
                           options.indent.code ' *********/' '\n'])];
        [c] = falcopt.generateTimer('indent', options.indent, 'header', false);
        code = [code, c, sprintf('\n')];
        timerName = [options.names.mex '_T'];
        code = [code, sprintf([options.indent.code 'timer ' timerName ';' '\n'])];
    end

    if options.verbose >= 2
        fprintf('. generating code\n');
    end
    code = [code, sprintf(['\n' options.indent.code '#include "mex.h"' '\n'])];
    code = [code, sprintf(['\n' options.indent.code '#include "string.h"' '\n'])];
    code = [code, sprintf(['\n' ...
                           options.indent.code '/****************' '\n' ...
                           options.indent.code ' * MEX Function *' '\n' ...
                           options.indent.code ' ****************/' '\n'])];
    % Generate conversion functions
    if strcmp(options.type, 'float')
        % from MATLAB (convert to single)
        code = [code, sprintf(['\n' ...
                               options.indent.code 'static ' options.inline ' void ' options.names.mex '_fromMATLAB(const double* v, ' dataType '* vt, unsigned int n) {' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<n; i++) { vt[i] = (' dataType ')v[i]; }' '\n'])];
        code = [code, sprintf([options.indent.code '};' '\n'])];
        % to MATLAB (convert to double)
        code = [code, sprintf([options.indent.code 'static ' options.inline ' void ' options.names.mex '_toMATLAB(const ' dataType '* v, double* vt, unsigned int n) {' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<n; i++) { vt[i] = (double)v[i]; }' '\n'])];
        code = [code, sprintf([options.indent.code '};' '\n'])];
    else
        code = [code, sprintf([options.indent.code 'static ' options.inline ' void ' options.names.mex '_copy(const double* v, double* vt, unsigned int n) {' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<n; i++) { vt[i] = v[i]; }' '\n'])];
        code = [code, sprintf([options.indent.code '};' '\n'])];
    end
    code = [code, sprintf([options.indent.code 'static ' options.inline ' void ' options.names.mex '_UINTtoMATLAB(const unsigned int* v, double* vt, unsigned int n) {' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<n; i++) { vt[i] = (double)v[i]; }' '\n'])];
    code = [code, sprintf([options.indent.code '};' '\n'])];
    % Generate MEX function
    nInputs = [1 1];
    if dims.w > 0
        nInputs = nInputs+1;
    end
    if options.ref
        nInputs = nInputs+2;
    end
    if options.terminalContraction
        nInputs = nInputs+1;
    end
    if options.indexContraction
        nInputs = nInputs+1;
    end
    nInputs = nInputs+[0,1];
    code = [code, sprintf(['\n' ...
                           options.indent.code  'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {' '\n'])];
    % Define variables
    index = 0;
    if strcmp(options.type, 'float')
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' x0[%i];' '\n'], dims.x)];
        if dims.w > 0
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic dataType ' w[%i];' '\n'], N*dims.w)];
        end
        if options.ref
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic dataType ' xref[%i];' '\n'], N*dims.x)];
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic dataType ' uref[%i];' '\n'], N*dims.x)];
        end
    else
        code = [code, sprintf([options.indent.code options.indent.generic 'const ' dataType '* x0 = mxGetPr(prhs[' num2str(index) ']);' '\n'])];
        if dims.w > 0
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic 'const ' dataType '* w = mxGetPr(prhs[' num2str(index) ']);' '\n'])];
        end
        if options.ref
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic 'const ' dataType '* xref = mxGetPr(prhs[' num2str(index) ']);' '\n'])];
            index = index+1;
            code = [code, sprintf([options.indent.code options.indent.generic 'const ' dataType '* uref = mxGetPr(prhs[' num2str(index) ']);' '\n'])];
        end
    end
    if options.indexContraction
        index = index+1;
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int contraction_idx = (unsigned int)(*mxGetPr(prhs[' num2str(index) ']));' '\n'])];
    end
    if options.terminalContraction
        index = index+1;
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic dataType ' contraction_c = (' dataType ')(*mxGetPr(prhs[' num2str(index) ']);' '\n'])];
        else
            code = [code, sprintf([options.indent.code options.indent.generic dataType ' contraction_c = *mxGetPr(prhs[' num2str(index) ']);' '\n'])];
        end
    end
    code = [code, sprintf([options.indent.code options.indent.generic dataType ' u[%i];' '\n'], N*dims.u)];
    code = [code, sprintf([options.indent.code options.indent.generic 'int* flag;' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* nitOut = mxCreateDoubleMatrix(1,1, mxREAL);' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int nit[1];' '\n'])];
    if options.maxIt > 0
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* lineSearch_nitOut = mxCreateDoubleMatrix(1,%i, mxREAL);' '\n'], options.maxIt)];
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int lineSearch_nit[%i];' '\n'], options.maxIt)];
    end
    if options.debug > 1
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' x[%i];' '\n'], N*dims.x)];
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' fval[1];' '\n'])];
    end
    if options.debug > 2
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' optimval[1];' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' feasval[1];' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic dataType ' meritval[1];' '\n'])];
    end
    if options.debug > 1
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* xOut = mxCreateDoubleMatrix(%i,%i,mxREAL);' '\n'], dims.x, N)];
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* fvalOut = mxCreateDoubleMatrix(1,1,mxREAL);' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* uOut = mxCreateDoubleMatrix(%i,%i,mxREAL);' '\n'], dims.u, N)];
    end
    if options.debug > 2
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* optimvalOut = mxCreateDoubleMatrix(1,1,mxREAL);' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* feasvalOut = mxCreateDoubleMatrix(1,1,mxREAL);' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* meritvalOut = mxCreateDoubleMatrix(1,1,mxREAL);' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic 'char** fieldnames;' '\n'])];
    if options.timing
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* t = mxCreateDoubleMatrix(1,1,mxREAL);' '\n'])];
    end
    if options.maxIt > 0
        code = [code, sprintf([options.indent.code options.indent.generic 'mxArray* lsStruct;' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
    
    % Check and process inputs
    if options.debug > 0
        code = [code, sprintf([options.indent.code options.indent.generic '/* Checking and processing inputs */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'if(nrhs < %i || nrhs > %i) {' '\n' ...
                               options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidInput", "Requires bewteen %i and %i inputs."); }' '\n'], nInputs(1), nInputs(2), nInputs(1), nInputs(2))];
    else
        code = [code, sprintf([options.indent.code options.indent.generic '/* Processing inputs */' '\n'])];
    end
    index = 0;
    % Initial state
    if options.debug > 0
        code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != %i) {' '\n' ... 
                               options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Initial state x0 has invalid dimension. Needs to be a column vector of size %i."); }' '\n'], dims.x, dims.x)];
    end
    if strcmp(options.type, 'float')
        code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_fromMATLAB(mxGetPr(prhs[' num2str(index) ']), x0, %i);' '\n'], dims.x)];
    end
    % Disturbance
    if dims.w > 0
        index = index+1;
        if options.debug > 0
            code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != %i) {' '\n' ... 
                                   options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Disturbance vector w has invalid dimension. Needs to be a column vector of size %i*%i."); }' '\n'], N*dims.w, N, dims.w)];
        end
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_fromMATLAB(mxGetPr(prhs[' num2str(index) ']), w, %i);' '\n'], N*dims.w)];
        end
    end
    
    if options.ref
        % References (state)
        index = index+1;
        if options.debug > 0
            code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != %i) {' '\n' ... 
                                   options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "State reference vector xref has invalid dimension. Needs to be a column vector of size %i*%i."); }' '\n'], N*dims.x, N, dims.x)];
        end
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_fromMATLAB(mxGetPr(prhs[' num2str(index) ']), xref, %i);' '\n'], N*dims.x)];
        end
        % References (input)
        index = index+1;
        if options.debug > 0
            code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != %i) {' '\n' ... 
                                   options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Input reference vector uref has invalid dimension. Needs to be a column vector of size %i*%i."); }' '\n'], N*dims.u, N, dims.u)];
        end
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_fromMATLAB(mxGetPr(prhs[' num2str(index) ']), uref, %i);' '\n'], N*dims.u)];
        end
    end
    
    if options.indexContraction
        index = index+1;
        if options.debug > 0
            code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != 1) {' '\n' ... 
                                   options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Input contraction index idx has invalid dimension. Needs to be a scalar."); }' '\n'])];
        end
    end
    if options.terminalContraction
        index = index+1;
        if options.debug > 0
            code = [code, sprintf([options.indent.code options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != 1) {' '\n' ... 
                                   options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Input contraction factor c has invalid dimension. Needs to be a scalar."); }' '\n'])];
        end
    end
    
    % Warmstart
    index = index+1;
    code = [code, sprintf([options.indent.code options.indent.generic 'if(nrhs == %i) { /* If a warmstart is supplied */' '\n'], index+1)];
    if options.debug > 0
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'if(!mxIsDouble(prhs[' num2str(index) ']) || mxIsComplex(prhs[' num2str(index) ']) || mxGetN(prhs[' num2str(index) ']) != 1 || mxGetM(prhs[' num2str(index) ']) != %i) {' '\n' ... 
                               options.indent.code options.indent.generic options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidDimension", "Warmstart vector uinit has invalid dimension. Needs to be a column vector of size %i*%i."); }' '\n'], N*dims.u, N, dims.u)];
    end
    if strcmp(options.type, 'float')
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic options.names.mex '_fromMATLAB(mxGetPr(prhs[' num2str(index) ']), u, %i);' '\n'], N*dims.u)];
    else
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic options.names.mex '_copy(mxGetPr(prhs[' num2str(index) ']), u, %i);' '\n'], N*dims.u)];
    end
    code = [code, sprintf([options.indent.code options.indent.generic '} else {' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'for(i=0; i<%i; i++) { u[i] = 0.0; } /* Initialize to zero */' '\n'], N*dims.u)];
    code = [code, sprintf([options.indent.code options.indent.generic '}' '\n'])];
    % Initialize timer
    if options.timing
        code = [code, sprintf([options.indent.code options.indent.generic '/* Initialize timer */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'timer_init(&' timerName ');' '\n'])];
    end
    % Prepare outputs
    code = [code, sprintf(['\n' ...
                           options.indent.code options.indent.generic '/* Prepare outputs */' '\n'])];
    if options.debug > 0
        code = [code, sprintf([options.indent.code options.indent.generic 'if(nlhs < 1 || nlhs > 3){' '\n' ...
                               options.indent.code options.indent.generic options.indent.generic 'mexErrMsgIdAndTxt("' options.names.mex ':InvalidOutput", "Requires bewteen 1 and 3 outputs."); }' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic 'if(nlhs >= 2) {' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'plhs[1] = mxCreateNumericMatrix(1,1, mxINT32_CLASS, mxREAL);' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'flag = (int*)mxGetData(plhs[1]);' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic '} else { ' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'flag = (int*)mxMalloc(sizeof(int));' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic '}' '\n'])];
    % Call function
    code = [code, sprintf(['\n' ...
                           options.indent.code options.indent.generic '/* Run optimization algorithm */' '\n'])];
    if options.timing
        code = [code, sprintf([options.indent.code options.indent.generic 'timer_start(&' timerName '); /* Start timer */' '\n'])];
    end
    
    code = [code, sprintf([options.indent.code options.indent.generic 'flag[0] = ' options.names.fun '(x0, u'])];
    if dims.w > 0
        code = [code, sprintf(', w')];
    end
    if options.ref
         code = [code, sprintf(', xref, uref')];
    end
    if options.indexContraction
        code = [code, sprintf(', contraction_idx')];
    end
    if options.terminalContraction
        code = [code, sprintf(', contraction_c')];
    end
    if options.debug > 1
        code = [code, sprintf(', x, fval')];
    end
    code = [code, sprintf(', nit')];
    if options.maxIt > 0
        code = [code, sprintf(', lineSearch_nit')];
    end
    if options.debug > 2
        code = [code, sprintf(', optimval, feasval, meritval')];
    end
    code = [code, sprintf([');' '\n'])];
    if options.timing
        code = [code, sprintf([options.indent.code options.indent.generic 'timer_stop(&' timerName '); /* Stop timer */' '\n'])];
    end
    % Process outputs
    code = [code, sprintf(['\n' ...
                           options.indent.code options.indent.generic '/* Process outputs */' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'plhs[0] = mxCreateDoubleMatrix(%i,1,mxREAL);' '\n'], dims.u)];
    if strcmp(options.type, 'float')
        code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(u, mxGetPr(plhs[0]), %i);' '\n'], dims.u)];
    else
        code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(u, mxGetPr(plhs[0]), %i);' '\n'], dims.u)];
    end
    code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_UINTtoMATLAB(nit, mxGetPr(nitOut), 1);' '\n'])];
    if options.maxIt > 0      
        code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_UINTtoMATLAB(lineSearch_nit, mxGetPr(lineSearch_nitOut), %i);' '\n'], options.maxIt)];
    end
    if options.debug > 1
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(x, mxGetPr(xOut), %i);' '\n'], N*dims.x)];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(fval, mxGetPr(fvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(u, mxGetPr(uOut), %i);' '\n'], N*dims.u)];
        else
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(x, mxGetPr(xOut), %i);' '\n'], N*dims.x)];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(fval, mxGetPr(fvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(u, mxGetPr(uOut), %i);' '\n'], N*dims.u)];
        end
    end
    if options.debug > 2
        if strcmp(options.type, 'float')
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(optimval, mxGetPr(optimvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(feasval, mxGetPr(feasvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_toMATLAB(meritval, mxGetPr(meritvalOut), 1);' '\n'])];
        else
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(optimval, mxGetPr(optimvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(feasval, mxGetPr(feasvalOut), 1);' '\n'])];
            code = [code, sprintf([options.indent.code options.indent.generic options.names.mex '_copy(meritval, mxGetPr(meritvalOut), 1);' '\n'])];
        end
    end
    code = [code, sprintf([options.indent.code options.indent.generic 'if(nlhs >= 3) {' '\n'])];
    infoFields = {'iterations'};
    if options.maxIt > 0
        infoFields = [infoFields {'lineSearch'}];
    end
    if options.timing
        infoFields = [infoFields {'time'}];
    end
    % Store additional stuff in info struct
    if options.debug > 1
        infoFields = [infoFields {'u', 'x', 'fval'}];
    end
    if options.debug > 2
        infoFields = [infoFields {'optimval', 'feasval', 'meritval'}];
    end    
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'fieldnames = mxCalloc(%i, sizeof(char*)); /* Allocate memory for storing pointers */' '\n'], length(infoFields))];
    for i=1:length(infoFields)
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'fieldnames[%i] = (char*)mxMalloc(sizeof(char)*%i);' '\n'], i-1, max(cellfun(@length, infoFields))+1)];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'memcpy(fieldnames[%i],"%s",sizeof("%s"));' '\n'], i-1, infoFields{i}, infoFields{i})];
    end
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'plhs[2] = mxCreateStructMatrix(1,1,%i,(const char**)fieldnames);' '\n'], length(infoFields))];
    for i=1:length(infoFields)
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxFree(fieldnames[%i]);' '\n'], i-1)];
    end
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,0,nitOut); /* Number of iterations */' '\n'])];
    if options.maxIt > 0
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'fieldnames[0] = (char*)mxMalloc(sizeof(char)*%i);' '\n'], length('iterations')+1)];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'memcpy(fieldnames[0],"iterations",sizeof("iterations"));' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'lsStruct = mxCreateStructMatrix(1,1,1,(const char**)fieldnames);' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxFree(fieldnames[0]);' '\n'], i-1)];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(lsStruct,0,0,lineSearch_nitOut);' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,1,lsStruct);' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxFree(fieldnames);' '\n'])];
    if options.timing
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxGetPr(t)[0] = timer_getTime(&' timerName ');' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,2,t); /* Runtime */' '\n'])];
    end
    if options.debug > 1
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,3,uOut); /* u */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,4,xOut); /* x */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,5,fvalOut); /* fval */' '\n'])];
    end
    if options.debug > 2
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,6,optimvalOut); /* optimval */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,7,feasvalOut); /* feasval */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic 'mxSetFieldByNumber(plhs[2],0,8,meritvalOut); /* meritval */' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic '}' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'if(nlhs < 2) { mxFree(flag); }' '\n'])];
    code = [code, sprintf([options.indent.code '}' '\n'])];
    
    if options.verbose >= 2
        fprintf('...done\n');
    end
    if options.verbose >= 1
        fprintf('Code successfully generated\n');
    end

end

