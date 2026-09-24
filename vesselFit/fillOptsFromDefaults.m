function opts = fillOptsFromDefaults(opts, def)
% FILLOPTSFROMDEFAULTS  Fill every opts field left unset (missing, or [] -- same "use the default"
% convention as fitVessel.m's own fillOptsDefaults) from a fully-populated default struct; anything the
% caller DID set is left untouched. Shared generic-loop core of the self-populating-default-opts /
% no-arg-call convention's own fillOptsDefaults() helpers (see
% .bass/pkm/patterns/self-populating-default-opts.md) -- extracted 2026-07-24 after the identical loop
% was found duplicated near-verbatim across fitMultiVessel.m/fitPatchVessels.m/fitVesselTimeSeries.m/
% movieIrfIm.m/drawVessel.m/fitVesselTimeSeriesDiag.m/showVessel.m. Callers still layer their own
% field-specific type-casts/normalizations/asserts on top of the result -- this only handles the
% generic "missing or empty -> default" fill; it does not know which fields are REQUIRED (no
% unconditional default) or which need special handling (conditional defaults, nested structs).
%
% INPUT
%   opts : the caller-supplied opts struct (any subset of def's fields set, rest missing or []).
%   def  : a fully-populated default opts struct, e.g. defaultOpts()'s own return value.
%
% OUTPUT
%   opts : same struct, with every field of def copied in wherever opts lacked it or it was empty.
%
%   opts = fillOptsFromDefaults(opts, def)
%
% NO-ARG CALL: none -- this is a pure positional-argument plumbing primitive (same convention as
% circlePolygonIntersectionArea.m/gaussianModel.m), not a user-facing opts-driven entry point itself.

fn = fieldnames(def);
for i = 1:numel(fn)
    name = fn{i};
    if ~isfield(opts,name) || isempty(opts.(name)); opts.(name) = def.(name); end
end
end
