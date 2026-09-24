function v = perRunVal(x, r)
% PERRUNVAL  Run r's value of a fitVesselPatchTimeSeries.m per-parameter field.
%
%   v = perRunVal(x, r)
%
% fitVesselPatchTimeSeries.m stores a value that DIFFERS between runs as a {1 x nRun} cell and a value
% that is the same in every run ONCE, bare (see that file's collapseIdentical) -- so "this run's
% value" is x{r} in the first case and x itself in the second. Every reader that indexes a fit
% result by run goes through here rather than assuming a cell.
%
% See also fitVesselPatchTimeSeries, buildGaussianFitDiag, plotVessels.
    if iscell(x); v = x{r}; else; v = x; end
end
