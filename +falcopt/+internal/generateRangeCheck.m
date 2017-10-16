%% generateRangeCheck Internal helper function to generate multiple range checks
%
% s = falcopt.internal.generateRangeCheck(range, name)
%  or
% s = falcopt.internal.generateRangeCheck(range, name, precision)
%
% Generates a string
% Given a scalar "num" returns a string representation for C.
% The precision (optional) determines how the number if printend.
%  'integer' means integer precision (default)
%  'double' means double-precision floating-point format
%   e.g. num = 101 becomes '101.0',
%        num = 3.141592653589793116... becomes s = '3.141592653589793' and
%        num = 33e-6 becomes s = '3.3e-05'.
%  'single' means single-precision floating-point format, additionally an "f" is appended
%   e.g. num = 101 becomes '101.0f',
%        num = 3.141592653589793116... becomes s = '3.14159265f' and
%        num = 33e-6 becomes s = '3.3e-05f'.
%

% Copyright (c) 2017 Damian Frick <falcopt@damianfrick.com>
%                    Giampaolo Torrisi <giampaolo.torrisi@gmail.com>
%                    Tommaso Robbiani <tommasro@student.ethz.ch>
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
function check = generateRangeCheck(varargin)
    precisions = {'integer', 'unsigned integer'};
    defaultPrecision = 'integer';

    p = inputParser;
    p.addRequired('range', @isnumeric);
    p.addRequired('name', @ischar);
    p.addParameter('precision', defaultPrecision, @(c)(ischar(c) && any(strcmp(c, precisions))));
    p.parse(varargin{:});
    options = p.Results;
    range = options.range;
    
    if strcmp(options.precision, 'unsigned integer') && any(range) < 0
        throw(MException());
    end
    % TODO Make sure range is sorted in ascending order
    
    % Get equality ranges
    eq = [true options.range(1:end-1)+1 < range(2:end)] & [range(1:end-1)+1 < range(2:end) true];
    i = 1;
    ranges = {};
    while i <= length(range)
        if eq(i)
            ranges{end+1} = range(i);
            i = i +1;
        else
            idx = i-1+find(eq(i+1:end), 1, 'first');
            if isempty(idx)
                idx = length(range);
            end
            ranges{end+1} = [range(i) range(idx)];
            i = idx+1;
        end
    end
    
    % Transform numerical ranges into checks
    for i=1:length(ranges)
        if length(ranges{i}) == 1
            ranges{i} = sprintf([options.name ' == ' falcopt.internal.num2str(ranges{i}, 'integer') ]);
        elseif strcmp(options.precision, 'unsigned integer') && ranges{i}(1) == 0
            ranges{i} = sprintf([options.name ' <= ' falcopt.internal.num2str(ranges{i}(2), 'integer')]);
        else
            if strcmp(options.precision, 'unsigned integer') && ranges{i}(1) == 0
                ranges{i} = sprintf([options.name ' <= ' falcopt.internal.num2str(ranges{i}(2), 'integer')]);
            else
                ranges{i} = sprintf(['(' options.name ' >= ' falcopt.internal.num2str(ranges{i}(1), 'integer') ') && (' options.name ' <= ' falcopt.internal.num2str(ranges{i}(2), 'integer') ')']);
            end
        end
    end
    
    check = strjoin(ranges, ' || ');
end
