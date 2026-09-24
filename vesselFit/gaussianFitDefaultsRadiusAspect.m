function d = gaussianFitDefaultsRadiusAspect(im2d, voxSz2, seedXY)
% GAUSSIANFITDEFAULTSRADIUSASPECT  gspec.defaultsFun for the 'radius,aspectRatio' parameterization --
% calls gaussianFitDefaults.m (UNTOUCHED, still shared directly by fitPatchVessels.m/fitMultiVessel.m
% for their own sx,sy-only seed-relative bounds) then converts its sx/sy start/lower/upper into
% radius/aspectRatio space. start converts EXACTLY (a single point: radius0=sqrt(sx0*sy0),
% aspectRatio0=sx0/sy0 -- since sx0==sy0==shape0 there, this reduces to radius0=shape0, aspectRatio0=1,
% a sensible isotropic starting guess). lower/upper are an APPROXIMATION: sx/sy's own bounds form an
% axis-aligned BOX ([sxLo,sxHi] x [syLo,syHi]); radius=sqrt(sx*sy) is monotonic increasing in both sx
% and sy, so its own outer bound is exactly the box's own two diagonal corners
% (sqrt(sxLo*syLo), sqrt(sxHi*syHi)); aspectRatio=sx/sy is increasing in sx, decreasing in sy, so ITS
% outer bound is the OTHER diagonal (sxLo/syHi, sxHi/syLo). Both are valid outer bounds for the box's
% own image under each formula, but NOT the tightest possible box in radius/aspectRatio space (the true
% feasible region there isn't axis-aligned either) -- a heuristic STARTING bound only, same spirit as
% gaussianFitDefaults.m's own heuristic, not meant to be exact. NOT used by fitPatchVessels.m/
% fitMultiVessel.m (still sx,sy-only, unaffected) -- see fitVessel.m's own gaussianMethodSpec.
%
% NOT symmetrized here (2026-08-17, see symmetrizeAspectRatioBound.m's own header) -- a purely
% mechanical sx,sy -> radius,aspectRatio unit conversion has no way to know whether theta will end up
% free (the ONE condition under which an asymmetric aspectRatio bound is actually unsafe, since
% canonicalizeRadiusAspectShape's own reciprocal only runs then) -- that's a fact only fitOneRun's own
% caller-supplied paramSpec knows, so the SYMMETRIZATION decision (not this conversion) is made there
% instead, conditionally, right after this function returns.
%
% Promoted out of fitVessel.m's own local function of the same name (2026-09-07) so there is exactly
% ONE copy of this formula -- this codebase's established convention (see gaussianFitDefaults.m's own
% header: it was promoted out of fitVessel.m for precisely this reason, after an earlier hand-copy
% drifted into fitPatchVesselsJoint.m). fitVessel.m itself now calls this shared file directly instead
% of keeping a local copy; fitVesselPatchTimeSeries.m's own 'path.heuristic'/'path(heuristic)' seed
% pathway also calls it directly (via resolveHeuristicSeedPeak), for the identical
% single-peak-heuristic-from-an-image need in the radius/aspectRatio parameterization.
    dSxSy = gaussianFitDefaults(im2d, voxSz2, seedXY);
    d.start = rmfield(dSxSy.start, {'sx','sy'});
    d.lower = rmfield(dSxSy.lower, {'sx','sy'});
    d.upper = rmfield(dSxSy.upper, {'sx','sy'});
    d.start.radius = sqrt(dSxSy.start.sx * dSxSy.start.sy);
    d.start.aspectRatio = dSxSy.start.sx / dSxSy.start.sy;
    d.lower.radius = sqrt(dSxSy.lower.sx * dSxSy.lower.sy);
    d.upper.radius = sqrt(dSxSy.upper.sx * dSxSy.upper.sy);
    d.lower.aspectRatio = dSxSy.lower.sx / dSxSy.upper.sy;
    d.upper.aspectRatio = dSxSy.upper.sx / dSxSy.lower.sy;
end
