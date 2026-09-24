function [lo, hi] = additiveSeedBound(seedVal, spec)
% ADDITIVESEEDBOUND  Seed-relative ADDITIVE bound window -- a parameter bounded by adding an offset to
% an already-converged seed value, rather than by a generic patch-scale heuristic. The additive half of
% the seed-relative bound vocabulary; see multiplicativeSeedBound.m for the other half.
%
% WHICH PARAMETERS -- the ones whose natural bound is an OFFSET in the parameter's own units: position
% (x0/y0, mm from the seed) and orientation (theta, radians from the seed). A parameter whose natural
% bound is a RATIO (amplitude, radius, aspectRatio, background) belongs to multiplicativeSeedBound.m
% instead.
%
% Promoted out of fitPatchVessels.m's own local copy (2026-09-03) once a SECOND caller needed the
% identical formula -- fitVesselPatchTimeSeries.m's own per-component opts.model.<type>(k).seedBounds --
% per this codebase's established promote-once-2+-callers-need-it convention (getNestedField.m,
% phaseScrambleRun.m, seedRelativeGaussianBounds.m, principalAxisProjection.m, ...). Unchanged in
% substance from that copy. Promoting rather than duplicating is deliberate here specifically: this file
% family has been bitten four separate confirmed times by a hand-duplicated bound formula drifting stale
% from the authoritative one -- see [[project_olsfaa2_fitpatchvessels]]'s own "KEEP IN SYNC" history.
%
% SCALAR AND VECTOR SPECS ARE DIFFERENT FORMULAS, deliberately (NOT the same formula with v1==v2):
%   scalar d  -> [seedVal-d, seedVal+d]   (SYMMETRIC around the seed -- d is a half-width)
%   [d1 d2]   -> [seedVal+d1, seedVal+d2] (ASYMMETRIC -- BOTH are ADDED, mirroring
%                multiplicativeSeedBound.m's own [seedVal*v1, seedVal*v2] exactly)
% A caller wanting the scalar case's own symmetric window via the vector form must therefore pass
% [-d, d] explicitly, negative d1 included.
%
% UNBOUNDING -- d (or d1/d2) = Inf works correctly with NO special-casing at all, unlike
% multiplicativeSeedBound.m's own .a/.b: adding +/-Inf always gives -Inf/+Inf regardless of seedVal's
% own sign or magnitude, so additiveSeedBound(seed, Inf) is exactly (-Inf, Inf).
%
% NO NO-ARG CALL, deliberately -- a pure positional-argument numeric primitive with no opts struct,
% the same documented exemption from the self-populating-default-opts convention that gaussianModel.m/
% buildIsochromatGrid.m/reduceIsochromatFull.m/canonicalizeEllipseShape.m already carry.
%
% INPUT
%   seedVal : scalar, the already-converged seed value the window is centred on/offset from.
%   spec    : scalar d (symmetric half-width), or 1x2 [d1 d2] (both added). Inf allowed either way.
%
% OUTPUT
%   lo, hi  : the bound interval. Ascending by construction for the scalar form (d>0 assumed); for the
%             vector form, ascending iff d1<=d2, which is the caller's own responsibility -- unlike
%             multiplicativeSeedBound.m, no min/max reordering happens here, because an additive offset
%             pair has an unambiguous intended order (lower offset first) that a sign flip cannot
%             scramble the way a multiplication by a negative seed can.
%
%   [lo, hi] = additiveSeedBound(2.5, 0.8);        % -> 1.7, 3.3
%   [lo, hi] = additiveSeedBound(2.5, [-0.2 1.0]); % -> 2.3, 3.5
%
% See also multiplicativeSeedBound, seedRelativeGaussianBounds, fitPatchVessels.

    if isscalar(spec)
        lo = seedVal - spec;
        hi = seedVal + spec;
    else
        lo = seedVal + spec(1);
        hi = seedVal + spec(2);
    end
end
