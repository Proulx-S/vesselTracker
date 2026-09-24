function [lo, hi] = multiplicativeSeedBound(seedVal, spec)
% MULTIPLICATIVESEEDBOUND  Seed-relative MULTIPLICATIVE bound window -- a parameter bounded by scaling
% an already-converged seed value, rather than by a generic patch-scale heuristic. The multiplicative
% half of the seed-relative bound vocabulary; see additiveSeedBound.m for the other half.
%
% WHICH PARAMETERS -- the ones whose natural bound is a RATIO rather than an offset: amplitude (a),
% size (radius), shape (aspectRatio) and background level (b). A parameter whose natural bound is an
% OFFSET in its own units (position x0/y0 in mm, orientation theta in radians) belongs to
% additiveSeedBound.m instead.
%
% Promoted out of fitPatchVessels.m's own local copy (2026-09-03) alongside additiveSeedBound.m -- see
% that file's own header for the full rationale (second caller, plus this file family's four confirmed
% "KEEP IN SYNC" bound-formula drift incidents). Unchanged in substance from that copy.
%
% SCALAR AND VECTOR SPECS ARE THE SAME FORMULA here, unlike additiveSeedBound.m:
%   scalar s -> v1=1/s, v2=s, i.e. [seedVal/s, seedVal*s]
%   [v1 v2]  -> [seedVal*v1, seedVal*v2]
% so the scalar case genuinely IS the vector case with v1=1/s, v2=s substituted in.
%
% min/max OF THE TWO ENDPOINTS, not v1/v2 in the order given -- required to stay a valid ASCENDING
% interval in two independent situations: (1) seedVal is negative (a and b are genuinely bipolar -- an
% amplitude can be a negative "dip" -- unlike radius/aspectRatio, which are positive by construction),
% which flips the product order; (2) an asymmetric [v1 v2] caller is not required to already know which
% of the two products ends up smaller.
%
% UNBOUNDING NEEDS CALLER-SIDE SPECIAL-CASING for a bipolar parameter, and this function deliberately
% does NOT do it: at SCALAR spec=Inf the general formula divides by Inf, collapsing the lower endpoint
% to +/-0 rather than preserving the seed's own sign (a negative seed gives [-Inf,-0]; a seed of exactly
% 0 degenerates to the single point [0,0]). The correct meaning of "unbounded" is parameter-specific --
% (-Inf,Inf) for a bipolar amplitude, [0,Inf) for a non-negative background/size -- and only the CALLER
% knows which parameter it is holding, so the caller special-cases scalar Inf (see fitPatchVessels.m's
% own fitOnePatchSimultaneous). A [v1 v2] VECTOR spec is never special-cased even when v1/v2 is itself
% Inf -- the min/max formula below already handles that correctly for either sign of seedVal (e.g.
% [-Inf Inf] gives the full real line, by explicit construction rather than by accident).
%
% NO NO-ARG CALL, deliberately -- same documented exemption as additiveSeedBound.m (pure positional
% numeric primitive, no opts struct).
%
% INPUT
%   seedVal : scalar, the already-converged seed value the window is scaled around.
%   spec    : scalar s>1 (symmetric-in-ratio), or 1x2 [v1 v2] (asymmetric, applied directly).
%
% OUTPUT
%   lo, hi  : the bound interval, always ascending (min/max of the two endpoints -- see above).
%
%   [lo, hi] = multiplicativeSeedBound(4, 2);        % ->  2, 8
%   [lo, hi] = multiplicativeSeedBound(-4, 2);       % -> -8, -2   (sign-safe, hence min/max)
%   [lo, hi] = multiplicativeSeedBound(4, [0.5 3]);  % ->  2, 12
%
% See also additiveSeedBound, symmetrizeAspectRatioBound, seedRelativeGaussianBounds, fitPatchVessels.

    if isscalar(spec); v1 = 1/spec; v2 = spec; else; v1 = spec(1); v2 = spec(2); end
    e1 = seedVal*v1;
    e2 = seedVal*v2;
    lo = min(e1, e2);
    hi = max(e1, e2);
end
