%% multFlops multiply each element of struct a with scalar b.
%
% c = falcopt.internal.multFlops( a, b)
%
% Inputs:
% a    - struct
% b    - scalar
% Outputs:
% c    - struct
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
function c = multFlops(a, b)
    if(~isscalar(b))
        error('second term must be scalar')
    end
    fields = fieldnames(a);
    for f=1:length(fields)
        c.(fields{f}) = b*a.(fields{f});
    end
end