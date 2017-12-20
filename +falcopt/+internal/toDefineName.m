%% toDefineName Transform string 'name' into all upper case (for defines),
%   where '_' are appropriately added
% 
% str = falcopy.internal.toDefineName(name)

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
function str = toDefineName(varargin)
    p = inputParser;
    p.addRequired('name', @ischar);
    p.parse(varargin{:});
    options = p.Results;
    name = options.name;
    
    str = regexprep(name, '([A-Z]+|[0-9]+)', '_$1'); % Add '_' in front of continguous upper case characters or numbers
    if str(1) == '_' && name(1) ~= '_' % Remove '_' if it was added to the start of the name
        str = str(2:end);
    end
    str = strrep(str, '__', '_'); % Get rid of duplicated '_'
    str = upper(str); % Make all upper case

end