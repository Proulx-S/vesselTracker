function tsIm = loadNiftiTs(fList, opts)
% LOADNIFTITS  Load one or more preprocessed 4D NIfTI timeseries files (ONE FILE PER RUN) into the
% ts-track struct shape the vessel pipeline reads (vessel.tsIm in huMoMain2's own tsImToVessel.m /
% promoteTsIm.m), but for the FULL image rather than a vessel patch:
%
%   tsIm.im     {1 x nRun} cell, each [ny x nx x 1 x T] double -- loadNiftiVol.m orientation, i.e.
%               rows = NIfTI dim 2, cols = NIfTI dim 1 (FreeSurfer MRIread's own .vol convention, which
%               every crop/x/y index in this pipeline assumes -- see loadNiftiVol.m's own NOTE).
%   tsIm.vSize  [3 x 1] voxel size in mm (NIfTI pixdim 1:3, in NIfTI dim order: vSize(1) is the
%               COLUMN spacing, vSize(2) the ROW spacing).
%   tsIm.dt     [1 x nRun] frame interval per run, in SECONDS (NIfTI pixdim 4, unit-converted).
%   tsIm.fName  {1 x nRun} the files, as given.
%
% Every run must share one voxel grid (image size and voxel size); frame counts may differ.
%
% INPUT
%   fList : char (one run) | cellstr (one entry per run) of .nii/.nii.gz paths.
%   opts  : struct --
%       .verbose  true | false -- print one line per file loaded (true).
%
% NO-ARG CALL -- loadNiftiTs() prints this opts listing and returns opts populated with defaults.
%
%   tsIm = loadNiftiTs({'run1.nii.gz','run2.nii.gz'});
%
% See also loadNiftiVol, makeVessel, drawVesselCenter.

    if nargin==0
        printOptsHelp();
        tsIm = defaultOpts();
        return
    end
    if nargin<2 || isempty(opts); opts = struct(); end
    opts = fillMissing(opts, defaultOpts());
    if ischar(fList) || isstring(fList); fList = {char(fList)}; end
    fList = cellfun(@char, fList, 'UniformOutput', false);
    nRun  = numel(fList);
    assert(nRun>0, 'loadNiftiTs:noFile', 'fList is empty.');

    tsIm = struct('im',{cell(1,nRun)}, 'vSize',[], 'dt',nan(1,nRun), 'fName',{fList});
    for r = 1:nRun
        f = fList{r};
        assert(exist(f,'file')==2, 'loadNiftiTs:missingFile', 'file not found: %s', f);
        if opts.verbose; fprintf('loadNiftiTs: run %d/%d <- %s\n', r, nRun, f); end
        info = niftiinfo(f);
        v = loadNiftiVol(f);                      % .vol is [ny x nx x nz x T] after the dim-1/2 swap
        vol = v.vol;
        assert(size(vol,3)==1, 'loadNiftiTs:notSingleSlice', ...
            '%s has %d slices -- this pipeline fits 2D single-slice patches.', f, size(vol,3));
        tsIm.im{r} = vol;
        vSize = info.PixelDimensions(1:3).';
        if isempty(tsIm.vSize)
            tsIm.vSize = vSize;
            sz1 = size(vol, [1 2]);
        else
            assert(all(abs(vSize - tsIm.vSize) < 1e-6), 'loadNiftiTs:voxSizeVaries', ...
                'voxel size differs between runs (run 1: %s, run %d: %s).', mat2str(tsIm.vSize.',4), r, mat2str(vSize.',4));
            assert(isequal(size(vol,[1 2]), sz1), 'loadNiftiTs:imSizeVaries', ...
                'image size differs between runs (run 1: %s, run %d: %s).', mat2str(sz1), r, mat2str(size(vol,[1 2])));
        end
        tsIm.dt(r) = frameIntervalSec(info);
    end
end

% ---------------------------------------------------------------------------
function dt = frameIntervalSec(info)
    % NIfTI pixdim(4) in the header's own time unit -> seconds.
    dt = double(info.PixelDimensions(min(4, numel(info.PixelDimensions))));
    if numel(info.PixelDimensions) < 4; dt = nan; return; end
    switch lower(char(info.TimeUnits))
        case 'second';       % as is
        case 'millisecond';  dt = dt / 1e3;
        case 'microsecond';  dt = dt / 1e6;
        otherwise            % 'None' or anything else -- leave the raw header value, flag it
            warning('loadNiftiTs:timeUnits', 'unrecognized NIfTI TimeUnits ''%s'' -- .dt left in raw header units.', char(info.TimeUnits));
    end
end

% ---------------------------------------------------------------------------
function opts = defaultOpts()
    opts.verbose = true;
end

function opts = fillMissing(opts, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})); opts.(fn{i}) = d.(fn{i}); end
    end
end

function printOptsHelp()
    fprintf('loadNiftiTs opts -- allowed value(s) per field (default in parentheses):\n');
    fprintf('  opts.verbose : true | false -- one line per file loaded   (true)\n\n');
end
