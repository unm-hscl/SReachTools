function prob = fht(varargin)
% FHT Create a first-hitting time problem.
%
%   See also: mht, tht
%
%   Copyright 2019 Adam Thorpe

prob = srt.problems.FirstHitting(varargin{:});
