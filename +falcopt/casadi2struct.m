%% casadi2struct(sx) given a casadi.SX matrix return the structure of the sparse unique matrix
% Inputs:
% - a casadi.SX matrix
% 
% Options:
% - structure: Char, can be 'unique', 'sparse' or 'dense'. Default: 'unique'. 
%
% Outputs:
% - const: A boolean. 1 if the matrix is composed of constant values
% - struct: structure of the matrix
%       .rows   the row indeces of nonzeros elements
%       .cols   the column indeces of nonzeros elements
%       .values the nonzeros elements
%       .mat    the matrix structure
%       .num    number of nonzeros elements
%
% Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory 
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
%
function [const, struct] = casadi2struct(varargin)
    p = inputParser;
    p.addRequired('sx_');
    p.addParameter('structure', 'unique');
    p.addParameter('errorName','function');
    p.parse(varargin{:});
    sx = p.Results.sx_;

    try 
        sx.is_constant();
    catch
        if isnumeric(sx)
            struct.stored.values = sx;
            const = 1;
            return;
        else
            error(['Error evaluating ' p.Results.errorName ': sx must be a casadi.SX matrix'])
        end
    end

    sx = sx.sparsify(1e-10);
    const = sx.is_constant();           % 1 if sx is no symbolic
    coli = sx.sparsity().get_col();     % vector of indices of nonzeros elements     
    nnz = size(coli,2);                 % number of nonzero elements
    [n,m] = size(sx);                   % size of sx matrix



    Mx = densify(sx); % transform sx to dense matrix


    % extract indeces of nonzeros elements
    st.rows = [];
    st.cols = [];
    st.values = [];
    for( i= 1:n)
        for( j = 1:m)
            if( ~is_equal(Mx(i,j),casadi.SX(0),2))||strcmp(p.Results.structure,'dense')
                st.rows = [ st.rows; i];
                st.cols = [ st.cols; j];
                st.values = [ st.values; simplify(Mx(i,j))];  
            end
        end
    end

    struct.access.rows = st.rows;
    struct.access.cols = st.cols;
    struct.access.values = st.values;
    struct.access.num = size(st.rows,1);
    
    if strcmp(p.Results.structure,'dense')
        nnz = length(st.values);
    end
    
    % build matrix structure
    nz = st.values;
    val = num2cell(st.values);
    M = zeros(n,m);
    ind = 1;

    if( nnz~= 0)
        st.values = [nz(1)];
    end
    sp = 0;
    for( i = 1:nnz)
        if( ~isnumeric( val{i}))
            switch p.Results.structure
                case {'sparse','dense'}
                    val{i} = ind;
                    ind = ind+1;
                otherwise
                    val{i} = ind;
                    for( k = i+1:nnz)
                        if( is_equal( nz(k),nz(i),2))
                            val{k} = ind;
                        end
                    end
                    ind = ind+1;
            end
        end
        M(st.rows(i),st.cols(i)) = val{i};

        % transform to unique matrix 
        if( i~=1)
            for( j = 1:size(st.values,1))
                if( is_equal( nz(i),st.values(j),2))
                    sp = 1;
                    break;
                end
            end
            if(( ~sp)||( ~strcmp(p.Results.structure, 'unique')))
                st.values = [st.values; nz(i)] ;
            end
            sp = 0;
        end
    end

    % if constant transform in double matrix
    if( const)
        nz = st.values;
        tmp = struct.access.values;
        st.values = zeros(size(nz,1),1);
        struct.access.values = zeros(size(nz,1),1);
        for( i= 1:size(st.values,1))
            st.values(i) = to_double(nz(i));
            struct.access.values(i) = to_double(tmp(i));
        end
    end

    st.mat = M;
    st.num = size(st.values,1);
    struct.stored = st;
end