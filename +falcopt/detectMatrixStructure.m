%% detectMatrixStructure Detects structure in a square matrix M
% 
% [elements, info] = detectMatrixStructure(M)
% or
% [elements, info] = detectMatrixStructure(M, ...) [with options]
% 
% Returns TODO
%
% the following options are available:
%  .structure - A string that determines the structure of M. Default is 'sparse'.
%                only affects how the structure of M is interpreted and utilized.
%                the possibilities are:
%                'dense' for no structure,
%                'sparse' for ommitted zero elements (only non-zero elements of M are considered)
%                'unique' for ommitted zero and reccurring elements (recurring elements of M are assumed to be the same)
%                'ordered' for ommitted zero and reccurring elements including permutation of elements (elements of M are stored in ascending order of their values)
%                'indexed' for ommitted zero and recurring elements including permutation of elements (elements of M are stored according to their values as indices starting from 1)
%  .symmetric - A boolean, if true then M is considered symmetric and code is generate to exploit this. Default: false.
%  .verbose   - Level of procedural output of this function
% If fields of options structs are ommitted, the default is used.
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
function info = detectMatrixStructure(varargin)
    structures = {'dense', 'sparse', 'unique', 'ordered', 'indexed'};
    defaultStructure = 'sparse';
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('M', @(x)(isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x)))));
    p.addParameter('structure', defaultStructure, @(x)((ischar(x) && any(strcmp(x, structures)))));
    p.addParameter('symmetric', false, @islogical);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.parse(varargin{:});
    options = p.Results;
    M = options.M;
    
    if any(options.symmetric)
        warning('falcopt:MissingImplemementation:SymmetricMatrix', 'The special treatment of symmetric matrices has not yet been implemented. Ignoring...');
    end
    
    %% Check dimension of M
    dims.n = size(M,1);
    
    
    %% Pre-process data
    [c, r, v] = find(M'); % Find non-zero elements, transposed to switch from column-major to row-major
    elements.M.row = r(:)';
    elements.M.col = c(:)';
    elements.M.values = v(:)';
    elements.M.num = length(elements.M.values);
    switch options.structure
        case 'dense'
            throw(MException('falcopt:detectMatrixStructure:MissingImplementation', 'Structure ''dense'' has not yet been implemented. Use ''sparse''.'));
        case 'sparse'
            elements.M.data.row = elements.M.row;
            elements.M.data.col = elements.M.col;
            elements.M.data.num = elements.M.num;
            elements.M.data.indices = 1:elements.M.data.num;
        case 'unique'
            throw(MException('falcopt:detectMatrixStructure:MissingImplementation', 'Structure ''unique'' has not yet been implemented. Use ''sparse''.'));
        case 'ordered'
            throw(MException('falcopt:detectMatrixStructure:MissingImplementation', 'Structure ''ordered'' has not yet been implemented. Use ''sparse''.'));
        case 'indexed'
            throw(MException('falcopt:detectMatrixStructure:MissingImplementation', 'Structure ''indexed'' has not yet been implemented. Use ''sparse''.'));
    end
    % Identify the individual blocks
    if dims.n == size(M,2)
        idx = 1; % Running index of block
        elements.M.blocks.indices{1} = elements.M.data.indices(1);
        elements.M.blocks.numel(1) = 1;
        I = 2:elements.M.data.num;
        while ~isempty(I)
            % Add all elements with same row or column
            elements.M.blocks.indices{idx} = unique([elements.M.blocks.indices{idx}, find(any(repmat(elements.M.row,elements.M.blocks.numel(idx),1) == repmat(elements.M.row(elements.M.blocks.indices{idx})',1,elements.M.num),1) ...
                                                                                    | any(repmat(elements.M.col,elements.M.blocks.numel(idx),1) == repmat(elements.M.col(elements.M.blocks.indices{idx})',1,elements.M.num),1))]);
            if length(elements.M.blocks.indices{idx}) == elements.M.blocks.numel(idx) % If we haven't added any new elements, start a new block
                idx = idx+1;
                elements.M.blocks.indices{idx} = I(1);
                elements.M.blocks.numel(idx) = 1;
                I = I(2:end);
            else % Remove added indices
                elements.M.blocks.numel(idx) = length(elements.M.blocks.indices{idx});
                I = I(all(repmat(I,elements.M.blocks.numel(idx),1) ~= repmat(elements.M.blocks.indices{idx}',1,length(I)),1));
            end
        end
        % Extract dimensions of blocks
        elements.M.blocks.num = length(elements.M.blocks.indices);
        for i=1:elements.M.blocks.num
            elements.M.blocks.size(i) = max(length(unique(elements.M.row(elements.M.blocks.indices{i}))), length(unique(elements.M.col(elements.M.blocks.indices{i}))));
            elements.M.blocks.indices{i} = elements.M.data.indices(elements.M.blocks.indices{i});
        end
    end
    info.elements = elements;
    
    % Store structure in info
    info.structure.M.type = options.structure;
    info.structure.M.mat = M;
    info.structure.M.num = elements.M.num;
    % TODO: always return "ordered" mat (or more)
    if strcmp(options.structure, 'dense') || strcmp(options.structure, 'sparse')
        % Replace all non-zero elements with 1
        info.structure.M.mat(info.structure.M.mat ~= 0) = 1;
    elseif strcmp(options.structure, 'unique')
        % Replace all non-zero elements with 1-to-...
        % TODO, CHECK
        info.structure.M.mat = info.structure.M.mat';
        info.structure.M.mat(info.structure.M.mat ~= 0) = 1:info.structure.M.num;
        info.structure.M.mat = info.structure.M.mat';
    elseif strcmp(options.structure, 'ordered')
        % TODO
    elseif strcmp(options.structure, 'indexed')
        % Nothing needs to be done
    end

end

