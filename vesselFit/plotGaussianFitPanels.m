function ax = plotGaussianFitPanels(diag, titleStr, opts)
% PLOTGAUSSIANFITPANELS  Six-panel gaussian-fit example figure (manuscript Fig4, panels xa-xf): time-
% averaged vessel patch, fitted isochromat-resolution model image, voxel-resolution residual (measured -
% predicted) image, TWO-SIDED radial profile of measured data, TWO-SIDED radial profile of the model's
% own isochromat samples, and measured-vs-predicted scatter -- one column each, TOP TO BOTTOM (xa/xb/xc
% down column 1, xd/xe/xf down column 2). xa/xb share ONE zero-floored CLim, ceilinged at the COMBINED
% max of BOTH panels' own data (im for xa, Zfit for xb) rather than xa's im alone -- the fit can overshoot
% the measured peak on a noisy patch, so xb's own data can set the ceiling too; xc is its OWN diverging
% colormap/symmetric CLim about zero (util's own colormap_divergingHue.m by default, a CIE-LCH-built
% diverging map, called directly with an explicit blue/red hue pair -- the most commonly used
% diverging-colormap choice, NOT colormap_blueNeutralRed.m, which is reserved for activation maps
% specifically -- see xc's own inline comment below; caller-overridable via opts.residualColormap, see
% INPUT) -- a residual is signed, not directly comparable to xa/xb's plain intensity scale. All
% three of xa/xb/xc get a colorbar + a 1-mm scale bar (decorateImagePanel, local to this file).
%
% xd/xe are BOTH two-sided by raw image left/right of center (NOT mirrored -- every voxel/isochromat of
% the WHOLE patch shows exactly once, see their own inline comments below) and radius- (not major-axis-)
% normalized, and BOTH share the SAME XLim (opts.xlimScale, see INPUT) and YLim ([0, data max], harmonized
% across the two so they read on one common scale, and starting at 0 since both series are physically
% non-negative intensities) with gridlines on -- xd shows raw MEASURED data, xe shows the model's own
% ISOCHROMAT sample values PLUS a continuous curve of diag.modelFun's own 1D-mode profile (the SAME
% formula the fit itself used, evaluated at tight regular intervals across xe's own two-sided XLim --
% see xe's own inline comment below) (no measured/predicted comparison drawn directly on either; that
% comparison is xf's own job instead). All three of xd/xe/xf use white as their main color; xf's own unity line is
% solid white with the axes grid on (deliberately distinct styling from fitVesselDiag.m's own reuse of
% the SAME shared drawing primitives, whose calls are unaffected by any of this -- see
% drawRadialProfileMeasuredPredicted.m/drawRadialProfileModelContinuous.m/
% drawMeasuredVsPredictedScatter.m's own new optional trailing arguments).
%
% Takes the SAME DIAG STRUCT fitVesselDiag.m's own low-level calling form does (see that file's own file
% header DIAG STRUCT) -- build one via buildGaussianFitDiag.m (a plain fitVessel.m-shaped result) or
% buildGaussianFitDiagFromPeak.m (one peak of a fitPatchVessels.m result), or hand-build one directly.
% Reuses drawModelImage.m/drawRadialProfileMeasuredPredicted.m/drawRadialProfileModelContinuous.m/
% drawMeasuredVsPredictedScatter.m/effectiveDistanceFromFit.m -- all promoted out of fitVesselDiag.m's
% own inline panel-drawing code -- rather than re-deriving any of that math a second time.
%
% Deliberately NOT built as a vessel-mode wrapper, and carries none of fitVesselDiag.m's own
% printIt/stayOpen/multi-vessel machinery: this is a curated MANUSCRIPT EXAMPLE figure (one call, one
% figure, one vessel/peak the caller already chose), not a QC/batch tool. Always returns a single
% visible figure's 6 axes handles, for a caller (see makeFigx3.m, this project's own
% "N axes from one shared tiledlayout -> N independent svg panels" convention -- makeFigx5.m's
% own precedent) to copy/format/insert individually -- the figure itself is never
% closed here.
%
% INPUT
%   diag       : DIAG STRUCT (see fitVesselDiag.m's own file header) -- .im/.voxSz2/.modelFun/
%                .paramNames/.shapeParamNames/.fitVals/.validMask/.predicted/.nIso/.method.
%   titleStr   : figure Name / sgtitle (default '').
%   opts.xlimScale : positive scalar -- xd/xe's shared XLim is [-radius radius]*opts.xlimScale, where
%                    radius = sqrt(s1*s2) is the fit's own equivalent-area radius (default 4).
%   opts.scale     : positive scalar -- xc's own CLim is [-c c]/opts.scale, where c is xa/xb's own shared
%                    CLim ceiling (the combined max of xa's im and xb's Zfit, see xa/xb above), i.e. the
%                    residual color scale is a fixed FRACTION of that shared peak intensity, not of the
%                    residual's own (potentially noisy/outlier-driven) max (default 10).
%   opts.residualColormap : xc's own colormap -- [] (default): this figure's own blue/red diverging map
%                    (colormap_divergingHue({[250 300],[340 40]}, [], [], [], [], 256), see xc's own inline
%                    comment below); char/string: a 'fname' or 'fname(args...)' spec resolved via
%                    resolveColormapFcn.m, N=256 (e.g. 'turbo', or SAME convention plotVessels.m's own
%                    opts.colormap.OL.method uses -- 'colormap_divergingHue([],0,0,0,0.5)' -- to pass
%                    colormap_divergingHue's own hues/wNeutral/wTransition/lNeutral/cOuter positional args
%                    through, see that file's own header for the exact grammar); numeric [N x 3]: used
%                    verbatim as the colormap itself (e.g. a differently-tuned colormap_divergingHue call).
%
% OUTPUT
%   ax : [1x6] axes handles, in xa..xf order:
%          ax(1) xa -- time-averaged patch (plotVesselPatch.m, no contour)
%          ax(2) xb -- fitted isochromat-resolution model image
%          ax(3) xc -- voxel-resolution residual (measured - predicted) image
%          ax(4) xd -- two-sided radial profile of MEASURED data, every voxel of the whole patch
%                      (signed left(-)/right(+) of center by raw image x-position, radius-normalized)
%          ax(5) xe -- two-sided radial profile of the model's own ISOCHROMAT samples, every isochromat
%                      of the whole patch (same signed/radius-normalized convention as xd)
%          ax(6) xf -- measured vs. predicted scatter
%
% NO-ARG CALL -- plotGaussianFitPanels() (zero input arguments) prints opts help and returns the default
% opts struct, same convention as fitVessel.m/plotRoiFitFinal.m's own no-arg call -- see PRINTOPTSHELP.
%
%   [diagS, titleStr] = buildGaussianFitDiagFromPeak(vesselHuman{1}, 'tsImMotionCorrectedUp', 'patchVesselsSimultaneous');
%   ax = plotGaussianFitPanels(diagS, titleStr);
%   ax = plotGaussianFitPanels(diagS, titleStr, struct('xlimScale',6));

    if nargin==0
        printOptsHelp();
        ax = defaultOpts();
        return
    end
    if nargin<2 || isempty(titleStr); titleStr = ''; end
    if nargin<3 || isempty(opts); opts = struct(); end
    opts = fillOptsDefaults(opts);

    im = double(diag.im);
    [ny,nx] = size(im);
    voxSz2 = diag.voxSz2;
    nIso = diag.nIso;
    white = [1 1 1];

    nP = numel(diag.paramNames);
    vals = cell(1,nP);
    for p = 1:nP; vals{p} = diag.fitVals.(diag.paramNames{p}); end

    % 'gaussian's own paramNames{4:5} are radius/aspectRatio (this project's sole parameterization -- see
    % fitVessel.m's own PARAMETERIZATION note) -- gaussianModel.m has no such concept, so convert to
    % literal sx,sy here, ONCE, as a plain local variable (modelVals), used for every diag.modelFun call
    % below. 'parabolic' has no analogous reparameterization (Rx,Ry is already native), so modelVals==vals.
    if strcmpi(diag.method,'gaussian')
        sqrtAR = sqrt(vals{5});
        modelVals = [vals(1:3), {vals{4}*sqrtAR, vals{4}/sqrtAR}, vals(6:7)];
    else
        modelVals = vals;
    end

    figArgs = {'Color','k', 'Visible','on', 'Name',['plotGaussianFitPanels: ' titleStr]};
    try
        fig = figure(figArgs{:}, 'WindowStyle','docked');
    catch
        fig = figure(figArgs{:});
    end
    try; theme(fig,'dark'); catch; end
    t = tiledlayout(fig, 3, 2, 'Padding','tight', 'TileSpacing','tight');

    % --- xa (tile 1): time-averaged patch, no contour -------------------------------------------------
    axA = nexttile(t,1);
    [toX, toY] = plotVesselPatch(axA, im, voxSz2);

    % --- isochromat grid -- SAME construction as fitVessel.m's own buildModelEvaluator ----------------
    offs = ((1:nIso)-0.5)/nIso - 0.5;
    xIsoGridPix = (1:nx) + offs(:); yIsoGridPix = (1:ny) + offs(:);
    [XIsoPix,YIsoPix] = meshgrid(xIsoGridPix(:), yIsoGridPix(:));
    XIso = toX(XIsoPix); YIso = toY(YIsoPix);
    Zfit = diag.modelFun(XIso, YIso, modelVals{:});

    % xa/xb's shared CLim: forced to START AT 0 (not min(im(:))) -- an intensity image's own near-zero
    % background should read as near-black, not a lighter gray offset by whatever small positive baseline
    % the patch happens to have -- and ceilinged at the COMBINED max of xa's own data (im) and xb's own
    % data (Zfit), not xa's im alone, since the fit can overshoot the measured peak on a noisy patch.
    % Degenerate-image guard mirrors plotVesselPatch.m's own (broaden rather than leave a zero-width CLim).
    cLim = [0, max(max(im(:)), max(Zfit(:)))];
    if diff(cLim)==0; cLim(2) = 1; end
    set(axA, 'CLim',cLim, 'FontSize',14);
    title(axA, 'time-averaged patch', 'Color','w', 'FontSize',14);
    decorateImagePanel(axA, 'intensity');

    % --- xb (tile 3): fitted isochromat-resolution model image -- SAME cLim as xa above ------------------
    axB = nexttile(t,3);
    drawModelImage(axB, XIso(1,:), YIso(:,1), Zfit, cLim);
    set(axB, 'FontSize',14);
    title(axB, sprintf('isochromat-resolution model (nIso=%d)', nIso), 'Color','w', 'FontSize',14);
    decorateImagePanel(axB, 'intensity');

    % --- xc (tile 5): voxel-resolution residual (measured - predicted) -- reduceIsochromatFull.m gives
    % the model's own downsampled (voxel-resolution) prediction at EVERY voxel (not just diag.validMask,
    % unlike diag.predicted -- see that file's own header), so the residual is defined over the whole
    % patch, matching xa/xb's own full-patch extent. UNLIKE xa/xb, this is its OWN diverging colormap +
    % SYMMETRIC CLim about zero: a residual is a SIGNED quantity, not directly comparable to xa/xb's own
    % plain intensity scale, so sharing their grayscale/CLim (this panel's original behavior) would have
    % made "no residual" and "dim signal" visually indistinguishable.
    %
    % Colormap: opts.residualColormap ([] default) -- util's own colormap_divergingHue.m (already on path
    % -- doIt_human.m's own environment setup puts util/ there), called DIRECTLY with an explicit blue/red
    % hue pair -- NOT colormap_blueNeutralRed.m (util PR #5, 2026-08-19: that legacy wrapper is now reserved for
    % ACTIVATION MAPS specifically, Seb's own convention; a residual map calls the generalized function
    % directly instead, even though it happens to use the SAME blue/red hues -- blue/red is simply the
    % most commonly used diverging-colormap choice per the earlier dataviz-skill-driven research
    % (matplotlib's own 'coolwarm'/'RdBu' defaults, crystallography's Fo-Fc convention, etc.), not
    % something reserved to this project's own activation-map vocabulary). Built in CIE LCH (searches
    % the sRGB gamut for the max-Delta-E equal-lightness pair within each hue's own sector, then ramps
    % L*/C* linearly in LCH space, converting to sRGB only at the very end), so it never suffers the
    % naive-RGB-interpolation "Mach band" artifact a raw RGB lerp through gray would.
    % CALLED WITH EXPLICIT [] ARGS, NEVER BARE -- colormap_divergingHue() with truly zero arguments does
    % NOT return its own documented positional defaults; nargin==0 instead runs that function's own
    % "self-demo" path (builds and leaves open an invisible demo figure, saves a PNG to pwd, and returns
    % a DIFFERENT, demo-specific parameterization) -- confirmed directly (2026-08-19, on
    % colormap_blueNeutralRed.m before this generalization; colormap_divergingHue.m inherits the same
    % dispatch): the demo's own 'baseline' cmap is NOT bit-equal to the true positional-default cmap.
    % This is a real footgun in that shared utility's own no-arg convention (inconsistent with this
    % codebase's own self-populating-default-opts pattern elsewhere), left AS IS (a shared tool's public
    % no-arg contract is a bigger, cross-project call) -- flagged to Seb rather than silently worked
    % around with no explanation. The interpolation math itself was audited and is correct: a synthetic
    % check confirmed the default (auto lNeutral, "pure chroma onset") L* profile is perfectly monotonic
    % from the near-white outer tip down to the inner max-separation lightness with zero sign reversals,
    % and chroma is exactly 0 across the whole neutral zone.
    %
    % CLim: +-cLim(2)/opts.scale, where cLim(2) is xa/xb's own shared CLim ceiling (the combined max of
    % xa's im and xb's Zfit, see xa/xb above) -- a FIXED FRACTION of that shared peak intensity (Seb's own
    % explicit choice), not of the residual's own max (this panel's original behavior, replaced -- that
    % was vulnerable to a single outlier voxel dominating the whole color range) and not noise-sigma-
    % normalized either (also considered, explicitly declined) -- so the scale stays anchored to xa/xb's
    % own signal amplitude and is directly adjustable via opts.scale.
    axC = nexttile(t,5);
    ZfitVox = reduceIsochromatFull(Zfit, nIso, ny, nx);
    residualVox = im - ZfitVox;
    residLim = cLim(2) / opts.scale; if residLim==0; residLim = 1; end
    cLimResidual = [-residLim residLim];
    drawModelImage(axC, toX(1:nx), toY(1:ny), residualVox, cLimResidual, resolveResidualColormap(opts.residualColormap));
    set(axC, 'FontSize',14);
    title(axC, 'residual (measured - predicted)', 'Color','w', 'FontSize',14);
    decorateImagePanel(axC, 'residual');

    % --- shared xd/xe geometry: BOTH two-sided by raw image left/right of center (NOT mirrored -- every
    % voxel/isochromat of the WHOLE patch shows exactly once) and radius- (not major-axis-) normalized.
    % Deliberately NOT effectiveDistanceFromFit.m (that primitive is unsigned/major-axis-normalized by
    % design, still used unchanged by fitVesselDiag.m's own panel 3) -- a distinct, non-reused metric,
    % computed inline. Sign is taken from each point's own raw IMAGE x-position relative to the fit's x0
    % (negative = left of center, positive = right of center, in the image's own x-axis -- NOT the fit's
    % rotated major axis), so a left-of-center and a right-of-center point are two genuinely different
    % data points here, never one folded onto the other. Magnitude is the SAME elliptical-normalized rho
    % effectiveDistanceFromFit.m itself uses, but scaled by the EQUIVALENT-AREA radius sqrt(s1*s2) (==
    % diag.fitVals.radius exactly, for 'gaussian') instead of that primitive's own major-axis sMajor --
    % "mid-way between long and short axis", not major-axis-biased.
    s1 = modelVals{4}; s2 = modelVals{5};
    radiusEq = sqrt(s1*s2);
    [Xpix,Ypix] = meshgrid(1:nx,1:ny);
    Xmm = toX(Xpix); Ymm = toY(Ypix);

    % xd's own population: every VOXEL of the whole patch (measured intensity).
    xrAll =  (Xmm-diag.fitVals.x0)*cos(diag.fitVals.theta) + (Ymm-diag.fitVals.y0)*sin(diag.fitVals.theta);
    yrAll = -(Xmm-diag.fitVals.x0)*sin(diag.fitVals.theta) + (Ymm-diag.fitVals.y0)*cos(diag.fitVals.theta);
    rhoAll = sqrt((xrAll/s1).^2 + (yrAll/s2).^2);
    sideSign = ones(ny,nx); sideSign(Xmm<diag.fitVals.x0) = -1;
    teffSigned  = sideSign .* rhoAll .* radiusEq;
    measuredAll = im(:);

    % xe's own population: every ISOCHROMAT of the whole patch (the model's own value there) -- SAME
    % sign/radius-normalization convention as xd above, just evaluated on the finer isochromat grid
    % (XIso/YIso/Zfit, already built for xb) instead of the voxel grid.
    xrIso =  (XIso-diag.fitVals.x0)*cos(diag.fitVals.theta) + (YIso-diag.fitVals.y0)*sin(diag.fitVals.theta);
    yrIso = -(XIso-diag.fitVals.x0)*sin(diag.fitVals.theta) + (YIso-diag.fitVals.y0)*cos(diag.fitVals.theta);
    rhoIso = sqrt((xrIso/s1).^2 + (yrIso/s2).^2);
    sideSignIso = ones(size(XIso)); sideSignIso(XIso<diag.fitVals.x0) = -1;
    teffIsoSigned = sideSignIso .* rhoIso .* radiusEq;
    valsIsoAll = Zfit(:);

    % Shared axes limits: XLim is opts.xlimScale radii either side of center (NOT data-driven -- a fixed,
    % adjustable window so xd/xe stay comparable across different vessels/calls); YLim starts at 0 (both
    % series here are intensities -- physically non-negative) and is harmonized across BOTH panels' own
    % data (measuredAll for xd, valsIsoAll for xe) so the two intensity scales read directly against
    % each other.
    xlimShared = [-radiusEq radiusEq] * opts.xlimScale;
    yAll = [measuredAll(:); valsIsoAll(:)];
    ylimShared = [0 max(yAll)];
    if ylimShared(2)<=ylimShared(1); ylimShared(2) = ylimShared(1)+1; end

    % --- xd (tile 2): two-sided radial profile of MEASURED data only (white circle, black edge) ----------
    axD = nexttile(t,2);
    hMeasD = drawRadialProfileMeasuredPredicted(axD, teffSigned(:), measuredAll, [], white);
    hMeasD.SizeData = 80;   % double drawRadialProfileMeasuredPredicted.m's own default marker size
                              % (40) -- xd/xf, top/bottom of the right column (Seb's own explicit ask,
                              % 2026-08-25), overridden here rather than in that shared primitive
                              % (also used unchanged by fitVesselDiag.m) to keep this a THIS-figure-
                              % only styling choice.
    xlim(axD, xlimShared); ylim(axD, ylimShared); grid(axD,'on');
    set(axD, 'PlotBoxAspectRatio',[1 1 1], 'FontSize',14);
    xlabel(axD, 'effective distance from center (mm, radius-normalized, left(-)/right(+))', 'Color','w', 'FontSize',14);
    ylabel(axD, 'intensity', 'Color','w', 'FontSize',14);
    title(axD, 'radial profile: measured (two-sided)', 'Color','w', 'FontSize',14);

    % --- xe (tile 4): two-sided radial profile of ISOCHROMAT samples (white dot marker) PLUS a
    % continuous curve of diag.modelFun's own 1D MODE profile (gaussianModel.m/parabolicModel.m's own
    % 1D MODE, see either file's own header) -- reuses the fit's OWN formula rather than re-deriving
    % the Gaussian/parabolic shape a second time here, evaluated at TIGHT regular intervals (401
    % samples) spanning the SAME two-sided [-radiusEq, radiusEq]*opts.xlimScale range as everything
    % else on this axes (Seb's own explicit ask, 2026-08-25).
    %
    % 1D MODE expects its own X argument as `rho*sMajor` (sMajor=max(s1,s2)), NOT this panel's own
    % radiusEq=sqrt(s1*s2) equivalent-area convention (see xd/xe's own sign/radius-normalization note
    % above) -- tCurve below is built directly on THIS panel's own radiusEq-scaled axis (so it reads
    % against the SAME x-axis the isochromat dots/xd's own measured markers already use), then
    % converted to the 1D MODE's own expected units (rho = tCurve/radiusEq, XforModel = rho*sMajor)
    % purely as an input-unit conversion for the call itself -- algebraically identical either way
    % (both reduce to the SAME underlying a*exp(-rho^2/2)+b / a*(1-rho^2) formula, see either model
    % file's own 1D MODE note), just expressed on a different-but-proportional distance axis.
    sMajor = max(s1,s2);
    tCurve = linspace(xlimShared(1), xlimShared(2), 401);
    curveVals = diag.modelFun((tCurve/radiusEq)*sMajor, [], modelVals{:});

    axE = nexttile(t,4);
    drawRadialProfileModelContinuous(axE, tCurve, curveVals, teffIsoSigned(:), valsIsoAll, white);
    xlim(axE, xlimShared); ylim(axE, ylimShared); grid(axE,'on');
    set(axE, 'PlotBoxAspectRatio',[1 1 1], 'FontSize',14);
    xlabel(axE, 'effective distance from center (mm, radius-normalized, left(-)/right(+))', 'Color','w', 'FontSize',14);
    ylabel(axE, 'intensity', 'Color','w', 'FontSize',14);
    title(axE, 'radial profile: isochromats (two-sided)', 'Color','w', 'FontSize',14);

    % --- xf (tile 6): measured vs. predicted scatter -- white marker, solid white unity line, grid on ----
    axF = nexttile(t,6);
    measured  = im(diag.validMask);
    predicted = diag.predicted(:);
    hMeasF = drawMeasuredVsPredictedScatter(axF, predicted, measured, white, white, '-', true);
    hMeasF.SizeData = 80;   % double drawMeasuredVsPredictedScatter.m's own default marker size (40)
                              % -- xd/xf, top/bottom of the right column (Seb's own explicit ask,
                              % 2026-08-25), same rationale as xd's own override above.
    set(axF, 'FontSize',14);
    xlabel(axF, 'predicted intensity', 'Color','w', 'FontSize',14);
    ylabel(axF, 'measured intensity', 'Color','w', 'FontSize',14);
    title(axF, 'measured vs. predicted', 'Color','w', 'FontSize',14);

    try; sgtitle(t, titleStr, 'Interpreter','none', 'Color','w', 'FontSize',14); catch; end

    ax = [axA axB axC axD axE axF];
end

% ---------------------------------------------------------------------------
function decorateImagePanel(ax, colorbarLabel)
    % Colorbar + 1-mm scale bar, this figure's own fixed styling for xa/xb/xc -- local to this file
    % (unlike drawScaleBar.m itself, a generic reusable primitive) since the exact styling (white, 1 mm,
    % FontSize 14) is specific to THIS manuscript figure, not a general "caller owns the axes" contract.
    cb = colorbar(ax); cb.Color = 'w'; cb.FontSize = 14;
    cb.Label.String = colorbarLabel; cb.Label.Color = 'w';
    drawScaleBar(ax, 1, [1 1 1]);
    % drawScaleBar.m hardcodes its own scale-bar LABEL text at FontSize 7, with no parameter to
    % override -- blanket-overridden here, post-hoc, on every text object under ax (the panel's own
    % title, already FontSize 14 from its own explicit title(...) call above, plus the scale-bar
    % label just created) rather than adding a FontSize argument to that shared, generic primitive
    % for one caller's own styling choice. findall, not findobj: a Title defaults to
    % HandleVisibility='off', which findobj silently skips (see makeFigx3.m's own
    % applyLightThemeChildren for the same, previously-confirmed gotcha).
    set(findall(ax, 'Type','text'), 'FontSize', 14);
end

% ---------------------------------------------------------------------------
function cmap = resolveResidualColormap(spec)
    % opts.residualColormap resolution -- [] (default): this figure's own blue/red diverging colormap
    % (see xc's own inline comment above for why THIS hue pair/THIS function, not
    % colormap_blueNeutralRed.m); char/string: a 'fname' or 'fname(args...)' spec, resolved via
    % resolveColormapFcn.m (SAME grammar/convention plotVessels.m's own opts.colormap.OL.method already
    % uses, e.g. 'colormap_divergingHue([],0,0,0,0.5)' -- see that file's own header for the exact
    % positional-argument rules); numeric [N x 3]: used verbatim, caller's own already-built colormap
    % (e.g. a differently-tuned colormap_divergingHue call, or any other diverging map).
    if isempty(spec)
        cmap = colormap_divergingHue({[250 300],[340 40]}, [], [], [], [], 256);
    elseif ischar(spec) || isstring(spec)
        cmap = resolveColormapFcn(char(spec), '', 256);
    else
        assert(isnumeric(spec) && ismatrix(spec) && size(spec,2)==3, ...
            'plotGaussianFitPanels:badResidualColormap', ...
            'opts.residualColormap must be [] (default), a colormap ''fname''/''fname(args...)'' spec (char/string), or an [N x 3] numeric colormap matrix.');
        cmap = spec;
    end
end

% ---------------------------------------------------------------------------
function printOptsHelp()
    % Reference listing of opts fields/defaults -- printed by the NO-ARG CALL (see file header) only; a
    % real call (any nargin>0) fills in the same defaults silently.
    fprintf('plotGaussianFitPanels parameters -- allowed value(s) (default in parentheses):\n');
    fprintf('  opts.xlimScale : positive scalar -- xd/xe share XLim = [-radius radius]*xlimScale, radius = sqrt(s1*s2)  (4)\n');
    fprintf('  opts.scale     : positive scalar -- xc''s own CLim = [-c c]/scale, c = xa/xb''s shared CLim ceiling  (10)\n');
    fprintf('  opts.residualColormap : [] (colormap_divergingHue blue/red default) | ''fname''/''fname(args...)'' colormap spec (char/string, see resolveColormapFcn.m) | [N x 3] numeric colormap  ([])\n');
end

% ---------------------------------------------------------------------------
function opts = defaultOpts()
    opts.xlimScale = 4;
    opts.scale     = 10;
    opts.residualColormap = [];
end

% ---------------------------------------------------------------------------
function opts = fillOptsDefaults(opts)
    def = defaultOpts();
    fn = fieldnames(def);
    unknown = setdiff(fieldnames(opts), fn);
    assert(isempty(unknown), 'plotGaussianFitPanels:unknownOptsField', ...
        'unrecognized opts field(s): %s -- allowed fields: %s (see plotGaussianFitPanels() (no-arg call)).', ...
        strjoin(unknown, ', '), strjoin(fn, ', '));
    for i = 1:numel(fn)
        name = fn{i};
        if ~isfield(opts, name) || isempty(opts.(name)); opts.(name) = def.(name); end
    end
end
