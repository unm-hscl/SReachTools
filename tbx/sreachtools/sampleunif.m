function S = sampleunif(varargin)
% SAMPLEUNIF Generate uniform samples.

samples = fliplr(varargin);

c = cell(size(samples));
[c{:}] = ndgrid(samples{:});

for k = 1:length(c)
    c{k} = reshape(c{k}.', 1, []);
end

S = vertcat(c{:});
