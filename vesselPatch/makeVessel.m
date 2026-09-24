function vessel = makeVessel(tsIm, center, opts)
% MAKEVESSEL  Crop one square patch per vessel center out of a full-image ts-track (loadNiftiTs.m
% output) and return the vessel cell array fitVesselPatchTimeSeries.m fits -- the SAME struct shape
% huMoMain2's own tsImToVessel.m / getVesselRoi2.m / promoteTsIm.m produce, reduced to the fields the
% fit and its visualizers actually read:
%
%   vessel{v}.sId        char, subject/dataset id (opts.sId)
%   vessel{v}.label      char, vessel label (opts.label)
%   vessel{v}.dt         [1 x nRun] frame interval per run, seconds (from tsIm.dt)
%   vessel{v}.center     [1 x 2] the (x, y) pixel this patch was centered on, in full-image pixels
%   vessel{v}.cropXlim   [1 x 2] full-image column range of the patch
%   vessel{v}.cropYlim   [1 x 2] full-image row range of the patch
%   vessel{v}.tsIm.im    {1 x nRun} of [cropSz+1 x cropSz+1 x 1 x T]
%   vessel{v}.tsIm.vSize [3 x 1] mm
%   vessel{v}.tsIm.x/.y  = cropXlim/cropYlim (getVesselRoi2.m's own per-track copy)
%   vessel{v}.tsIm.fName {1 x nRun} source files
%   vessel{v}.tsIm.rois  empty [0 x 0] struct with the standard roi fields (label/mask/poly/im/offset/
%                        parent/mod) -- no ROI is drawn here; the fit uses the whole patch unless
%                        opts.includeMask/.excludeMask name a label on it.
%
% The crop window is getVesselRoi2.m's own: x = round(cx) + [-1 1]*cropSz/2 (same for y), so a patch
% is (cropSz+1) pixels wide and the clicked pixel is its exact center pixel. Requires an EVEN cropSz.
%
% INPUT
%   tsIm   : loadNiftiTs.m output.
%   center : [nVessel x 2] (x, y) pixel centers -- drawVesselCenter.m's output.
%   opts   : struct --
%       .cropSz  even positive integer, patch is cropSz+1 pixels square    (10 -> 11 x 11)
%       .sId     char                                                     ('example')
%       .label   {} = 'vessel01', 'vessel02', ... | cellstr, one per row of center   ({})
%
% NO-ARG CALL -- makeVessel() prints this opts listing and returns opts populated with defaults.
%
%   vessel = makeVessel(tsIm, center);
%   vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', fitOpts);
%
% See also loadNiftiTs, drawVesselCenter, fitVesselPatchTimeSeries.

    if nargin==0
        printOptsHelp();
        vessel = defaultOpts();
        return
    end
    if nargin<3 || isempty(opts); opts = struct(); end
    opts = fillMissing(opts, defaultOpts());
    assert(isstruct(tsIm) && isfield(tsIm,'im') && ~isempty(tsIm.im), 'makeVessel:noIm', ...
        'tsIm must be a loadNiftiTs.m output with a non-empty .im.');
    assert(isnumeric(center) && size(center,2)==2, 'makeVessel:badCenter', ...
        'center must be [nVessel x 2] (x, y) pixel coordinates.');
    cropSz = opts.cropSz;
    assert(isscalar(cropSz) && cropSz>0 && mod(cropSz,2)==0, 'makeVessel:badCropSz', ...
        'opts.cropSz must be an even positive integer (patch is cropSz+1 pixels square).');
    nV   = size(center, 1);
    nRun = numel(tsIm.im);
    [ny, nx] = size(tsIm.im{1}, [1 2]);
    labels = opts.label;
    if isempty(labels); labels = arrayfun(@(v) sprintf('vessel%02d', v), 1:nV, 'UniformOutput', false); end
    if ischar(labels); labels = {labels}; end
    assert(numel(labels)==nV, 'makeVessel:labelCount', 'opts.label has %d entries for %d centers.', numel(labels), nV);

    emptyRois = struct('label',{},'mask',{},'poly',{},'im',{},'offset',{},'parent',{},'mod',{});
    vessel = cell(1, nV);
    for v = 1:nV
        cx = center(v,1); cy = center(v,2);
        x = round(cx) + [-1 1]*cropSz/2;
        y = round(cy) + [-1 1]*cropSz/2;
        assert(x(1)>=1 && y(1)>=1 && x(2)<=nx && y(2)<=ny, 'makeVessel:cropOutOfBounds', ...
            'vessel %d (%s): a %dx%d patch around (x=%g, y=%g) falls outside the %dx%d image.', ...
            v, labels{v}, cropSz+1, cropSz+1, cx, cy, nx, ny);
        s = struct();
        s.sId      = opts.sId;
        s.label    = labels{v};
        s.dt       = tsIm.dt;
        s.center   = [cx cy];
        s.cropXlim = x;
        s.cropYlim = y;
        s.tsIm.im  = cell(1, nRun);
        for r = 1:nRun
            s.tsIm.im{r} = tsIm.im{r}(y(1):y(2), x(1):x(2), :, :);
        end
        s.tsIm.vSize = tsIm.vSize(:);
        s.tsIm.x     = x;
        s.tsIm.y     = y;
        s.tsIm.fName = tsIm.fName;
        s.tsIm.rois  = emptyRois;
        vessel{v} = s;
    end
    fprintf('makeVessel: %d vessel patch(es) of %dx%d pixels, %d run(s)\n', nV, cropSz+1, cropSz+1, nRun);
end

% ---------------------------------------------------------------------------
function opts = defaultOpts()
    opts.cropSz = 10;
    opts.sId    = 'example';
    opts.label  = {};
end

function opts = fillMissing(opts, d)
    fn = fieldnames(d);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}); opts.(fn{i}) = d.(fn{i}); end
    end
end

function printOptsHelp()
    fprintf('makeVessel opts -- allowed value(s) per field (default in parentheses):\n');
    fprintf('  opts.cropSz : even positive integer -- patch is cropSz+1 pixels square   (10)\n');
    fprintf('  opts.sId    : char                                                    (''example'')\n');
    fprintf('  opts.label  : {} = vessel01, vessel02, ... | cellstr, one per center   ({})\n\n');
end
