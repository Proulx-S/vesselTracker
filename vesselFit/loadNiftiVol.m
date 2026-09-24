function s = loadNiftiVol(f)
% LOADNIFTIVOL  niftiread-based stand-in for MRIread(f).vol, returned as a {.vol,.fspec}-shaped
% struct so callers that expect a loaded-MRI-like struct keep working without a FreeSurfer MATLAB
% library dependency.
%
% NOTE: FreeSurfer's MRIread returns .vol with dims 1/2 (row/col) SWAPPED relative to niftiread's
% native array -- verified empirically (permute(niftiread(f),[2 1 3 4]) matched MRIread(f).vol to
% floating-point equality on a real preprocessed volume from this dataset). The permute below
% applies that same convention so downstream indexing (crop x/y, COM, etc.) behaves identically to
% the original MRIread-based pipeline.
    v = double(niftiread(f));
    permDims = 1:ndims(v);
    permDims(1:2) = [2 1];
    v = permute(v, permDims);
    s = struct('vol', v, 'fspec', f);
