function predVox = reduceIsochromatFull(predIso, nIso, ny, nx)
% REDUCEISOCHROMATFULL  Same isochromat-averaging reduction as reduceIsochromat.m, but returns the FULL
% [ny x nx] voxel-resolution grid -- NOT restricted/flattened to a validMask. Needed wherever a peak's
% own contribution is SUBTRACTED from (or composited into) a full 2D image rather than compared only at
% validMask voxels (e.g. fitPatchVessels.m's own sequential peak-subtraction, fitPatchVesselsDiag.m's own
% diagnostic panels) -- validMask(true(ny,nx)) still linearizes to a COLUMN VECTOR in MATLAB (logical
% indexing never preserves 2D shape, even when every element is true), so reduceIsochromat.m itself
% cannot be reused for this "no restriction" case by simply passing an all-true mask.
%
%   predVox = reduceIsochromatFull(predIso, nIso, ny, nx)

    predVox = squeeze(mean(mean(reshape(predIso, nIso,ny,nIso,nx), 1), 3));
end
