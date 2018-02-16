%% findEqual
%

% Copyright (c) 2018, ETH Zurich, Automatic Control Laboratory 
%                    Damian Frick <falcopt@damianfrick.com>
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

function [K, indices] = findEqual(x)
    if isempty(x)
        indices = [];
        K = {};
    elseif ~iscell(x) || all(size(x) ~= [1,1])
        error('Only takes 1d cell arrays.');
    end

    idx = 1;
    indices = zeros(1,length(x));
    
    while ~isempty(idx)
        select = diag(indices == 0);
        select = double(select(:,indices == 0));
        indices(select*cellfun(@(c)(isequal(x{idx},c)), x(indices == 0))' == 1) = max(indices)+1;
        idx = find(indices == 0, 1, 'first');
    end
    
    if any(indices == 0)
        error('Something went wrong (this should not be possible).');
    end
    
    K = cell(1,max(indices));
    for i=1:max(indices)
        K{i} = find(indices == i);
    end
end