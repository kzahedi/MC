function [ v ] = binVector( varargin )
%BINMATRIX discretises a vector
%   Required parameter:  data
%   Optional parameters: min, max, bins
%   Ouput: Vector with values in [1, bins]

  p = inputParser;
  p.addRequired('v',        @isvector);
  p.addOptional('min',   0, @isscalar);
  p.addOptional('max',   0, @isscalar);
  p.addOptional('bins', 30, @isscalar);
  p.parse(varargin{:});
  p.Results;

  v    = p.Results.v;
  min_ = p.Results.min;
  max_ = p.Results.max;
  bins = p.Results.bins;

  if min_ == max_
    min_ = min(v);
    max_ = max(v);
  end
  
  v = min(bins, int64(1.0 + floor((v - min_) / (max_ - min_) * bins)));
  
end
