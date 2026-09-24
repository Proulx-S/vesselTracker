function predValid = reduceIsochromat(predIso, nIso, ny, nx, validMask)
% REDUCEISOCHROMAT  Average each voxel's own nIso x nIso isochromats (see buildIsochromatGrid.m) down to
% voxel resolution, then restrict to validMask -- in the SAME linear order frame(validMask) produces, so
% a caller can subtract/compare directly against real voxel data.
%
% Promoted out of fitVessel.m's own local copy -- see buildIsochromatGrid.m's own header for why. For
% the "no restriction, full [ny x nx] frame" case (e.g. peak-subtraction against the whole patch, or a
% diagnostic panel needing every voxel, not just the fit region), call reduceIsochromatFull.m instead --
% NOT this function with validMask=true(ny,nx): logical indexing never preserves 2D shape (even an
% all-true mask linearizes to a column vector), so that substitution would silently break any caller
% expecting a 2D result.
%
%   predValid = reduceIsochromat(predIso, nIso, ny, nx, validMask)

    predVox = reduceIsochromatFull(predIso, nIso, ny, nx);
    predValid = predVox(validMask);
end
