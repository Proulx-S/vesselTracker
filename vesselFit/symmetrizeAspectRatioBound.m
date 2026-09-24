function [lo, hi] = symmetrizeAspectRatioBound(lo, hi)
% SYMMETRIZEASPECTRATIOBOUND  Widens an aspectRatio [lo,hi] bound to be reciprocal-symmetric around 1
% (lo==1/hi) -- CONTAINS the original interval exactly (equal to it already when lo*hi==1), only ever
% widens, never shrinks: k=max(hi,1/lo), lo->1/k, hi->k.
%
% WHY THIS MATTERS: fitVessel.m's own canonicalizeRadiusAspectShape reciprocates aspectRatio whenever
% it fits <1 (aspectRatio -> 1/aspectRatio, radius unchanged, theta+=pi/2 -- the correct radius,
% aspectRatio-space analog of canonicalizeEllipseShape.m's own s1<->s2 exchange for Rx,Ry). That
% canonicalization only actually RUNS when radius, aspectRatio, AND theta are ALL free/estimated
% together (fitVessel.m's own canonicalizeThetaIfFree, its "allEstimated" check) -- in that case ONLY,
% a value legally inside a NON-symmetric [lo,hi] can reciprocate to something OUTSIDE it (x in [lo,hi]
% does not imply 1/x in [lo,hi] unless lo*hi==1), which is exactly how a genuinely bound-respecting
% fit can end up DISPLAYED as "exceeding" its own bound (Seb's own report, 2026-08-17).
%
% CALLER'S OWN RESPONSIBILITY -- deliberately NOT baked into the low-level sx,sy -> radius,aspectRatio
% bound conversion itself (gaussianFitDefaultsRadiusAspect.m / sxSyBoundsToRadiusAspect, fitVessel.m's
% own and its 2 duplicates in fitPatchVessels.m/fitMultiVessel.m). That conversion is a generic,
% mechanical unit change and has no way to know whether symmetrization is even needed -- whether theta
% will actually end up free is a fact only the CALLER knows (e.g. via opts.gaussian.mode.theta /
% opts.gaussian.fixed.theta), not something a shared low-level utility should silently impose as a
% policy (Seb's own call, 2026-08-17, after an earlier version of this fix put the widening inside
% that conversion directly, unconditionally -- correctly flagged as bounds a higher-level caller set
% getting silently changed by lower-level code). Call this ONLY when theta will actually be estimated
% jointly with the shape pair -- if theta is fixed (or in a different mode than radius/aspectRatio),
% canonicalizeRadiusAspectShape's own reciprocal never runs at all (see canonicalizeThetaIfFree's own
% allEstimated guard), so a non-symmetric aspectRatio bound is entirely safe there, and may be the
% caller's own deliberate, tighter choice.
%
%   [lo,hi] = symmetrizeAspectRatioBound(lo, hi)
%   % e.g. lo=0.02, hi=2.2 (asymmetric, lo*hi=0.044) -> lo=1/22.7, hi=22.7 (lo*hi=1)

    k = max(hi, 1/lo);
    lo = 1/k;
    hi = k;
end
