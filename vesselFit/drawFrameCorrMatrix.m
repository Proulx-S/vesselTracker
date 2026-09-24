function [hIm, hRow, hCol, hRunLines] = drawFrameCorrMatrix(ax, C, tAxis, cLim, cmap, crosshairColor, crosshairWidth, runEdges)
% DRAWFRAMECORRMATRIX  Draws a square frame-by-frame correlation matrix into a GIVEN axes, plus a
% MOVABLE row/column crosshair marking one "current" frame -- the one shared primitive for this
% pipeline's cross-frame similarity panel, in the same "caller owns the axes/figure lifecycle" family as
% plotVesselPatch.m/drawMaskOutline.m/drawModelImage.m (positional arguments, per-argument [] defaults,
% no opts struct, never creates a figure or axes of its own).
%
% This primitive does NOT compute C -- correlating WHAT (raw patch frames, per-frame fit residuals, the
% fitted model itself) over WHICH voxels is a caller-level scientific choice with genuinely different
% meanings, so it stays with the caller (see fitVesselTimeSeriesDiag.m's own CORRELATION MATRIX note)
% rather than being buried in a drawing primitive that would then have to grow an option per choice.
%
% CROSSHAIR -- the whole reason this returns handles rather than just drawing. In an animation the
% matrix itself is STATIC (it is a property of the whole series, not of one frame) and only the
% crosshair moves, so a caller re-renders nothing per frame -- it just reassigns the two lines' own
% constant coordinate:
%     set(hRow, 'YData',[tNow tNow]);   set(hCol, 'XData',[tNow tNow]);
% in tAxis units (see fitVesselTimeSeriesDiag.m's own frame loop). Deliberately NOT wrapped in a helper
% function: a two-property assignment is shorter and clearer at the call site than any wrapper, and a
% local function here would not be callable from another file anyway. Drawn as full-width/full-height
% line objects rather than xline/yline deliberately -- an xline's own position is settable too, but
% xline/yline attach to the axes' own ruler and are drawn ON TOP of every later child regardless of
% stacking order, which on a dark image panel reads as two opaque bars across the data; plain lines
% respect ordinary child order and can be restacked by the caller if it wants.
%
% INPUT
%   ax    : target axes (an existing axes handle).
%   C     : [n x n] correlation (or any square similarity) matrix. NOT required to be exactly symmetric
%           or to have a unit diagonal -- this function only ever DISPLAYS it, and a caller may
%           legitimately pass a matrix built with a masked/partial voxel set. NaNs render transparent
%           (AlphaData), so a frame with no finite data (e.g. an all-NaN censored frame) reads as a
%           genuine gap rather than as a spurious zero correlation.
%   tAxis : [1 x n] monotonic values labelling both axes -- e.g. frame index 1:n, or acquisition time in
%           seconds. Default 1:n. Used for the image's own XData/YData, so the crosshair position a
%           caller assigns to hRow/hCol is in THESE units, not pixel indices.
%   cLim  : [lo hi] color limits. Default [-1 1] -- the full correlation range, SYMMETRIC about zero and
%           FIXED rather than data-driven, so the same color means the same correlation across every
%           vessel/call (a data-driven range would make a structureless matrix look dramatic by
%           autoscaling to its own noise). Pass an explicit range to zoom in.
%   cmap  : [N x 3] colormap. Default: util's own colormap_divergingHue.m called with an explicit
%           blue/red hue pair -- a correlation is a SIGNED quantity about zero, the same reasoning
%           plotGaussianFitPanels.m's own residual panel uses (see that file's own xc note for why
%           colormap_divergingHue.m directly, and NEVER bare/no-arg).
%   crosshairColor : RGB triplet for the crosshair. Default [0 1 0] (green -- reads against both ends of
%           a blue/red diverging map, unlike white/black which collide with the map's own neutral and
%           its saturated extremes respectively).
%   crosshairWidth : crosshair LineWidth. Default 1.0.
%   runEdges : [] (default, nothing drawn) | vector of tAxis values at which to draw a STATIC boundary
%           line on both axes -- for a matrix spanning several concatenated runs, so run structure is
%           visible as such rather than being mistaken for signal. Drawn in a dim gray, dashed, visually
%           subordinate to the crosshair.
%
% OUTPUT
%   hIm       : the image object (for a caller that wants to restack/alpha it).
%   hRow      : the HORIZONTAL crosshair line (constant y = current frame).
%   hCol      : the VERTICAL crosshair line (constant x = current frame).
%   hRunLines : [1 x m] gobjects for the run-boundary lines ([] when runEdges is empty).
%
%   [hIm,hRow,hCol] = drawFrameCorrMatrix(ax, C, 1:nFrame);
%   set(hRow,'YData',[37 37]); set(hCol,'XData',[37 37]);   % per animation frame -- no redraw of C
%
% See also FITVESSELTIMESERIESDIAG, PLOTVESSELPATCH, DRAWMODELIMAGE.

    assert(nargin>=2, 'drawFrameCorrMatrix:tooFewArgs', 'ax and C are both required.');
    assert(isnumeric(C) && ismatrix(C) && size(C,1)==size(C,2) && ~isempty(C), ...
        'drawFrameCorrMatrix:badC', 'C must be a nonempty square 2-D numeric matrix (got %s).', mat2str(size(C)));
    n = size(C,1);

    if nargin<3 || isempty(tAxis); tAxis = 1:n; end
    assert(isnumeric(tAxis) && isvector(tAxis) && numel(tAxis)==n, 'drawFrameCorrMatrix:badTAxis', ...
        'tAxis must be a numeric vector with numel == size(C,1) (%d), got %s.', n, mat2str(size(tAxis)));
    tAxis = tAxis(:).';

    if nargin<4 || isempty(cLim); cLim = [-1 1]; end
    assert(isnumeric(cLim) && numel(cLim)==2 && cLim(2)>cLim(1), 'drawFrameCorrMatrix:badCLim', ...
        'cLim must be [lo hi] with hi>lo.');
    if nargin<5 || isempty(cmap); cmap = defaultDivergingMap(); end
    assert(isnumeric(cmap) && size(cmap,2)==3 && size(cmap,1)>=2, 'drawFrameCorrMatrix:badCmap', ...
        'cmap must be [N x 3] with N>=2.');
    if nargin<6 || isempty(crosshairColor); crosshairColor = [0 1 0]; end
    assert(isnumeric(crosshairColor) && numel(crosshairColor)==3, 'drawFrameCorrMatrix:badColor', ...
        'crosshairColor must be an RGB triplet.');
    if nargin<7 || isempty(crosshairWidth); crosshairWidth = 1.0; end
    assert(isnumeric(crosshairWidth) && isscalar(crosshairWidth) && crosshairWidth>0, ...
        'drawFrameCorrMatrix:badWidth', 'crosshairWidth must be a positive scalar.');
    if nargin<8; runEdges = []; end

    C = double(C);

    % A single-frame matrix has no extent for imagesc to spread over -- give it a nominal one so the
    % panel still renders (degenerate-input guard, same spirit as plotVesselPatch.m's own uniform-image
    % CLim broadening) rather than collapsing to a zero-width axes.
    if n==1
        xd = tAxis + [-0.5 0.5];
        hIm = imagesc(ax, xd, xd, C([1 1],[1 1]));
    else
        hIm = imagesc(ax, tAxis, tAxis, C);
    end
    set(hIm, 'AlphaData', ~isnan(get(hIm,'CData')));   % NaN -> transparent, a real gap not a fake zero
    axis(ax, 'image');
    set(ax, 'CLim', cLim, 'YDir', 'reverse');   % YDir reverse: row 1 at the TOP, so the matrix reads
                                                 % like a matrix (and like every image panel in this
                                                 % pipeline), not like a scatter plot
    colormap(ax, cmap);
    hold(ax, 'on');

    % Run boundaries FIRST, so the crosshair drawn below sits on top of them.
    hRunLines = gobjects(1,0);
    if ~isempty(runEdges)
        assert(isnumeric(runEdges) && isvector(runEdges), 'drawFrameCorrMatrix:badRunEdges', ...
            'runEdges must be a numeric vector of tAxis values.');
        lim = [min(tAxis) max(tAxis)];
        if n==1; lim = tAxis + [-0.5 0.5]; end
        for k = 1:numel(runEdges)
            e = runEdges(k);
            hRunLines(end+1) = plot(ax, lim, [e e], ':', 'Color',[0.55 0.55 0.55], 'LineWidth',0.75); %#ok<AGROW>
            hRunLines(end+1) = plot(ax, [e e], lim, ':', 'Color',[0.55 0.55 0.55], 'LineWidth',0.75); %#ok<AGROW>
        end
    end

    % Crosshair, initialised at the FIRST frame -- full-extent lines whose position the caller then
    % moves per animation frame (see moveFrameCorrCrosshair). XLim/YLim are read back off the axes
    % AFTER axis('image') rather than from tAxis, so the lines span exactly the drawn extent (imagesc
    % pads by half a pixel at each edge and axis('image') keeps that padding).
    xl = xlim(ax); yl = ylim(ax);
    hRow = plot(ax, xl, [tAxis(1) tAxis(1)], '-', 'Color',crosshairColor, 'LineWidth',crosshairWidth);
    hCol = plot(ax, [tAxis(1) tAxis(1)], yl, '-', 'Color',crosshairColor, 'LineWidth',crosshairWidth);
    xlim(ax, xl); ylim(ax, yl);   % the two plot() calls above can nudge the limits; pin them back
end

% ---------------------------------------------------------------------------
function cmap = defaultDivergingMap()
    % Blue/red diverging map for a signed quantity about zero -- colormap_divergingHue.m called with
    % EXPLICIT [] args, NEVER bare: that shared utility's own nargin==0 path runs a self-demo (opens a
    % figure, writes a PNG to pwd, returns a demo-specific parameterization), not its documented
    % positional defaults -- see plotGaussianFitPanels.m's own xc note, which documents this footgun in
    % full. Falls back to MATLAB's own built-in 'turbo'-free diverging option only if that utility is
    % genuinely not on the path, which for this pipeline means a broken environment setup, so it warns
    % rather than substituting silently.
    if exist('colormap_divergingHue','file')
        cmap = colormap_divergingHue({[250 300],[340 40]}, [], [], [], [], 256);
        return
    end
    warning('drawFrameCorrMatrix:noDivergingHue', ...
        ['colormap_divergingHue.m is not on the path (util/ missing from the environment setup) -- ' ...
         'falling back to MATLAB''s own parula, which is NOT diverging and will misrepresent the sign ' ...
         'of a correlation. Fix the path rather than relying on this.']);
    cmap = parula(256);
end
