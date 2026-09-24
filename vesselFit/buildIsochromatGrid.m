function [XIso, YIso] = buildIsochromatGrid(nIso, nx, ny, ctrCol, ctrRow, voxSz2)
% BUILDISOCHROMATGRID  Super-resolution ISOCHROMAT GRID position -- nIso x nIso sub-voxel positions per
% voxel (ODD, so one lands exactly at each voxel's own center), in physical mm, center-voxel-relative
% (see fitVessel.m's own POSITION CONVENTION). Position-only (independent of any fit parameter) -- built
% ONCE and reused across every model evaluation.
%
% Promoted out of fitVessel.m's own buildModelEvaluator (previously the grid-construction half of that
% closure) so it has exactly ONE home instead of being hand-copied per caller -- this exact formula had
% independently drifted (stale bound heuristics riding alongside a correct copy of this grid math) in
% fitPatchVesselsJoint.m and fitPatchVesselsDiag.m before the 2026-07-21 rework that introduced this file;
% promoting the position-only math here removes one axis of that duplication. fitVessel.m/
% fitMultiVessel.m/fitPatchVessels.m/fitPatchVesselsDiag.m all call this directly now.
%
%   [XIso, YIso] = buildIsochromatGrid(nIso, nx, ny, ctrCol, ctrRow, voxSz2)
%
% nIso=1 degenerates to exactly a single point at each voxel's own center (see reduceIsochromat.m).

    offs = ((1:nIso)-0.5)/nIso - 0.5;
    xIsoGridPix = (1:nx) + offs(:);
    yIsoGridPix = (1:ny) + offs(:);
    [XIsoPix,YIsoPix] = meshgrid(xIsoGridPix(:), yIsoGridPix(:));
    XIso = (XIsoPix-ctrCol)*voxSz2(2); YIso = (YIsoPix-ctrRow)*voxSz2(1);
end
