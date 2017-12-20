%% vec2strjoin Internal helper function to transform vectors into compact, readable strings for C
% s = falcopt.internal.vec2strjoin(v, sep, [options])
%
% Given a vector v returns a string representation of use in C code.
% The numbers in v are separated by a string provided in "sep".
%
% the following options are available:
% .precision - The precision (optional) determines how the number if printend.
%               Default: 'integer'. See also falcopt.internal.num2str
%
%
% Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory 
%                    Damian Frick <falcopt@damianfrick.com>
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
function s = vec2strjoin(varargin)
    precisions = {'double', 'single', 'integer'};
    defaultPrecision = 'integer';

    p = inputParser;
    p.addRequired('v', @isnumeric);
    p.addRequired('sep', @ischar);
    p.addOptional('precision', defaultPrecision, @(s)(ischar(s) && any(strcmp(s, precisions))));
    p.parse(varargin{:});
    options = p.Results;

    s = strjoin(cellfun(@(s)(falcopt.internal.num2str(s, options.precision)), num2cell(options.v), 'UniformOutput', false), options.sep);

end
