function s = vec2strjoin(varargin)
    precisions = {'double', 'single', 'integer'};
    defaultPrecision = 'double';

    p = inputParser;
    p.addRequired('v', @isnumeric);
    p.addRequired('sep', @ischar);
    p.addOptional('precision', defaultPrecision, @(s)(ischar(s) && any(strcmp(s, precisions))));
    p.parse(varargin{:});
    options = p.Results;

    s = strjoin(cellfun(@(s)(falcopt.num2str(s, options.precision)), num2cell(options.v), 'UniformOutput', false), options.sep);

end
