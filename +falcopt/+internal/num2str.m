%% num2str Internal helper function to transform numbers into compact, readable strings for C
%
% s = falcopt.internal.NUM2STR(num)
%  or
% s = falcopt.internal.NUM2STR(num, precision)
%
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
function s = num2str(varargin)
    precisions = {'double', 'single', 'integer'};
    defaultPrecision = 'integer';
    
    p = inputParser;
    p.addRequired('num', @(x)(isnumeric(x) && numel(x) == 1));
    p.addOptional('precision', defaultPrecision, @(s)(ischar(s) && any(strcmp(s, precisions))));
    p.parse(varargin{:});
    options = p.Results;
    num = options.num;
    
    minDigits = 3;
    if any(strcmp(options.precision, {'single', 'double'}))
        switch options.precision
            case 'double'
                nDigits = 16; % If we want 17, then we need more than double precision num
            case 'single'
                nDigits = 9;
        end
         % Determine number of needed digits
        dstr = num2str(abs(num)*10^(nDigits-fix(log10(abs(num)))), nDigits+1);
        d = regexp(dstr(1:end-1), '([0]+)$', 'start');
        if isempty(d)
            d = nDigits;
        else
            d = max(d-1, 1); % Max only has an effect if num = 0
        end
        if d <= minDigits && fix(log10(abs(num))) <= minDigits && fix(log10(abs(num))) > 0
           d = d+fix(log10(abs(num)));
        end
        s = sprintf(['%0.' num2str(d) 'g'], num);
        if isempty(strfind(s, '.'))
            pieces = strsplit(s, 'e');
            if length(pieces) == 1
                s = [pieces{1} '.0'];
            else
                s = [pieces{1} '.0e' pieces{2}];
            end
        end
        if strcmp(options.precision, 'single')
            s = [s 'f'];
        end
    elseif strcmp(options.precision, 'integer')
        % TODO incorporate length
        s = sprintf('%i', num);
    end

end
