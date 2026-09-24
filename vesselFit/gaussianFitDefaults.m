function d = gaussianFitDefaults(im2d, voxSz2, seedXY)
% GAUSSIANFITDEFAULTS  Default StartPoint/bounds heuristic for a single rotated-2D-Gaussian
% (gaussianModel.m) fit against im2d -- same formula fitVesselProfile.m's own Gaussian fit uses,
% generalized to accept an arbitrary (single-frame or time-averaged) im2d rather than always the
% run-pooled time-average. Bounds here are PATCH-CENTER-anchored (x0/y0 bounded by the whole crop, not
% by any particular peak's own seed) -- see seedRelativeGaussianBounds.m for the per-peak, seed-anchored
% variant multi-peak callers (fitPatchVessels.m's sequential loop, fitMultiVessel.m) actually use.
%
% Promoted out of fitVessel.m's own local `gaussianDefaults` so there is exactly ONE copy of this
% formula -- previously hand-copied (with independent, confirmed drift) into fitPatchVesselsJoint.m.
% fitVessel.m itself now calls this directly instead of keeping a local copy.
%
% x0/y0 are center-voxel-relative mm, matching fitVesselProfile.m's UNITS/ORIGIN convention (same
% ctrCol/ctrRow formula) -- REQUIRED for a 'fixed' source pulled from lumenGaussFit.x0/.y0 to mean the
% same physical coordinate as this function's own x0/y0.
%
%   d = gaussianFitDefaults(im2d, voxSz2, seedXY)
%   % d.start/.lower/.upper -- each a struct with fields a/x0/y0/sx/sy/theta/b (gaussianModel.m's own
%   % parameter names/order)

    sigma = 0.3;
    [ny,nx] = size(im2d);
    ctrCol = (nx+1)/2; ctrRow = (ny+1)/2;
    imSmooth = imgaussfilt(im2d, sigma);
    b0 = mean(imSmooth(:));
    % seedXY, ROUNDED, clamped into [1,nx]/[1,ny] before indexing -- a caller passing the ORIGINAL
    % detection seed (fitPatchVessels.m's own sequential loop, x0Seed=Xmm(peakRow,peakCol)) always lands
    % exactly on an integer pixel already in range, so this clamp is a no-op there. A caller passing a
    % FITTED (sub-pixel-precision) x0/y0 back-converted to pixel space -- fitMultiVessel.m's own joint
    % refit, seeded from a PRIOR stage's own converged position, not a raw detection -- has no such
    % guarantee: a peak that converged right at (or within half a pixel of) its own x0/y0 bound can
    % round to pixel 0 or nx+1/ny+1, indexing outside im2d entirely ("Index in position 1 is invalid" --
    % a real crash on real data, 2026-07-22). Clamping to the nearest valid edge pixel is the correct
    % behavior here regardless: a0's own job is just a representative LOCAL amplitude estimate near the
    % seed, and the nearest in-bounds pixel is the best available substitute for one that's fallen just
    % outside the patch.
    seedRowPix = min(max(round(seedXY(2)),1), ny);
    seedColPix = min(max(round(seedXY(1)),1), nx);
    a0 = max(imSmooth(seedRowPix,seedColPix) - b0, 0.01*range(imSmooth(:)));
    x0_0 = (seedXY(1)-ctrCol)*voxSz2(2); y0_0 = (seedXY(2)-ctrRow)*voxSz2(1);
    shapeLb = 0.1*min(voxSz2);   % was 0.25*min(voxSz2) -- Seb's own ask, loosened for narrower vessels
    shapeUb = max(nx*voxSz2(2), ny*voxSz2(1))/3;
    shape0  = 0.5*mean(voxSz2);
    % theta is UNBOUNDED (-Inf,Inf) -- a bounded theta can force gradient descent to traverse the WHOLE
    % parameter range to reach a solution that's actually just beyond the opposite bound, since the true
    % optimum and its periodic image can be numerically far apart despite describing the identical
    % ellipse (theta+pi, and theta+pi/2 with sx/sy swapped, both describe the same ellipse -- see
    % canonicalizeEllipseShape.m). Canonicalization happens ONLY AFTER the fit converges, never during
    % optimization. x0/y0 bounds are SIGNED (symmetric around the center-voxel origin).
    d.start = struct('a',a0, 'x0',x0_0, 'y0',y0_0, 'sx',shape0,  'sy',shape0,  'theta',0,    'b',b0);
    d.lower = struct('a',0,  'x0',(1-ctrCol)*voxSz2(2),  'y0',(1-ctrRow)*voxSz2(1),  'sx',shapeLb, 'sy',shapeLb, 'theta',-Inf, 'b',min(imSmooth(:)));
    d.upper = struct('a',10*range(imSmooth(:)), 'x0',(nx-ctrCol)*voxSz2(2), 'y0',(ny-ctrRow)*voxSz2(1), 'sx',shapeUb, 'sy',shapeUb, 'theta',Inf, 'b',max(imSmooth(:)));
end
