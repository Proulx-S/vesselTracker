function [vessel, diagOut] = fitVesselTimeSeriesDiag(vessel, fld, opts)
% FITVESSELTIMESERIESDIAG  QA movie for a PER-FRAME vessel fit (fitVesselTimeSeries.m / fitVessel.m with
% at least one mode='perFrame' parameter) -- a 4x5 panel figure ANIMATED over the timeseries frames:
% column 1 is four small per-FRAME panels that change with the movie, and columns 2-5 are ONE spanning
% panel showing the whole series' own frame-by-frame correlation matrix, static, with a crosshair that
% tracks the frame column 1 is currently showing.
%
% The complement to fitVesselDiag.m/fitPatchVesselsDiag.m, which diagnose ONE STILL: fitVesselDiag.m's
% own vessel-mode form already takes opts.run/opts.frame and renders one frame per call, and
% fitVessel.m's own opts.verbose already loops EVERY frame -- as one separate figure per frame (see that
% file's own diagFrames loop), which is exactly why it is unusable for a 350-frame run. This file is that
% loop turned into a single animated figure, plus the cross-frame context (the correlation matrix) that
% no per-frame still can carry by construction.
%
% WHAT IS ACTUALLY BEING QA'd -- a per-frame fit in this pipeline (see doIt_human.m's own motion-estimate
% section) typically frees only x0/y0 per frame and holds a/radius/aspectRatio/theta/b per VESSEL. So the
% fitted SHAPE is frozen and only its POSITION tracks the vessel frame to frame, and the question this
% figure exists to answer is "does that position actually track, on every frame, or does it jump/stick/
% drift?" -- which is a question about a SEQUENCE and cannot be answered from any single frame. Panels 1/2
% answer it spatially (is the contour on the vessel? is the residual structureless?), panels 3/4
% answer it radially, and the correlation matrix answers it for all frames at once.
%
% LAYOUT -- one tiledlayout(4,5), TileSpacing 'none' and Padding 'tight' (space between elements
% minimized, Seb's own ask), inset from the left by opts.leftMargin to hold tiles 1/2's own colorbars
% and tiles 3/4's own y-tick labels. NO PANEL CARRIES A TITLE: the colorbar labels and y-labels name
% every quantity, so titles were pure vertical cost. Tile numbers are the nexttile(t,N) call sites:
%   1  (row 1, col 1) -- the frame's own patch image (plotVesselPatch.m), MASKED to the fit region
%                     (AlphaData; see MASKED TO THE FIT REGION below), with THIS FRAME's fitted outline
%                     overlaid at ONE level: the contour enclosing opts.gaussContourPercIntegral percent
%                     of the Gaussian's integrated signal, 95% by default (gold,
%                     drawGaussianPeakContour.m). Colorbar 'meas', left. XLim/YLim flush to the image.
%   6  (row 2, col 1) -- the frame's own VOXEL-resolution residual image (measured - predicted), same
%                     masking and same flush limits, on its own SYMMETRIC diverging scale about zero --
%                     the SAME colormap the manuscript figure's residual panel uses
%                     (plotGaussianFitPanels.m:322's own colormap_divergingHue call, verbatim).
%                     Colorbar 'meas-pred', left. A well-tracked frame leaves noise here; a
%                     vessel-shaped dipole means the position is off on THIS frame.
%   11 (row 3, col 1) -- MEASURED intensity as small white dots, fit-mask voxels only, against signed
%                     distance from the fitted center; PLUS the PREDICTED profile as a continuous curve
%                     in the SAME GOLD as tile 1's outline, so the outline and the model trace read as
%                     one object seen two ways. The curve is gaussianModel.m's own 1-D mode -- the
%                     formula the optimizer itself used -- not a re-derivation. y-label 'meas,model'.
%                     NO x ruler: it shares tile 16's (see HORIZONTAL EXTENT below).
%   16 (row 4, col 1) -- the same voxels' RESIDUAL as small white dots, with a y=0 reference line and an
%                     ASYMMETRIC y-axis (the data's own range, not +-max -- real residuals are not
%                     symmetric about zero and forcing symmetry throws away half the vertical
%                     resolution). y-label 'meas-pred'. Carries the shared x ruler for tiles 11/16.
%   2  (rows 1-4, cols 2-5, ONE spanning tile) -- the frame-by-frame correlation matrix over the whole
%                     selected series (drawFrameCorrMatrix.m), parula, CLim [0 1], colorbar 'frame xcorr'
%                     along the BOTTOM, with a RED horizontal + vertical crosshair at the frame column 1
%                     is showing,
%                     moving down/right as the movie plays. Every ruler decoration is off except the x
%                     ticks/tick labels -- the matrix is symmetric, so a y ruler duplicating the x one
%                     is clutter. STATIC content: it is a property of the series, not of a frame, so it
%                     is drawn ONCE and only the crosshair moves (see PERFORMANCE below). See
%                     CORRELATION MATRIX for what is actually correlated.
%
%   6  (rows 1-4, col 6, ONE spanning tile) -- PER-FRAME DISSIMILARITY, axes SWAPPED: the metric on x,
%                     frame on y, sharing the correlation matrix's own YLim and YDir so frame N sits at
%                     the same height in both and the two are read across row by row. What is plotted is
%                     the mean (or opts.dissim.stat median) of 1-rho from that frame to every other,
%                     diagonal excluded -- "how unlike the rest of the series is this frame", so a spike
%                     is a candidate frame to censor. 1-rho is Seb's OWN correlation->distance
%                     conversion, taken from bassReg2/QAdendrogram.m:162
%                         linkage(squareform(1-rho), 'average')   % Convert correlation to distance
%                     applied there to the matrix bassReg2/xCorrQA.m:142 builds (cross-frame corr over
%                     masked voxels) -- the same object this figure draws, and the same pipeline that
%                     wrote censor files from the resulting clusters. That old pair is also where this
%                     panel's placement comes from: QAdendrogram drew its dendrogram
%                     'Orientation','right', immediately right of the matrix with frames running
%                     vertically. (xCorrQA.m:151 is also where CLim [0 1] comes from -- imagesc(rho,[0 1]).)
%                     A DENDROGRAM is deliberately not what is drawn: the ask was a metric AS A FUNCTION
%                     OF FRAME, and a tree is not that. The per-frame mean is the pointwise reduction of
%                     the very distances that clustering consumed pairwise.
%
% HORIZONTAL EXTENT -- tiles 11/16 share tiles 1/6's own x extent in mm (opts.xlimScale=1 makes it an
% exact match), so column 1 reads as one vertically-aligned stack rather than four unrelated panels.
% Checked rather than assumed, because the 1-D axis is radius-NORMALIZED and so is not literally the
% image's x coordinate: across vessels 1/5/13/20/30/36 the largest |teff| over MASKED voxels is
% 1.26-2.00 mm against a 2.20 mm image half-extent, so nothing is clipped. (Out-of-mask voxels do reach
% 3.1-5.4 mm -- they are no longer drawn.) The panel WIDTHS line up too, which is why tiles 1/2 get
% drawThinColorbar.m rather than a native colorbar(ax,'westoutside'): a native one shrinks its own axes
% and would leave tiles 1/2 narrower than tiles 3/4.
% Built with the "oversampled columns + spanning nexttile" idiom plotIRFmat.m:83-112 already uses in this
% repo (nexttile(t, idx, [nRow nCol])) -- deliberately NOT a new plotVessels.m route: that dispatcher
% selects ONE route per call and each route creates its own figure/tiledlayout with tiles claimed by bare
% linear index (opts.grid there means a grid of VESSELS, not of panel slots), so hosting a mixed-kind
% spanning layout would mean threading a per-tile route+span concept through most of it. That is exactly
% the call plotVessels.m's own SCOPE note says to flag rather than make silently, and the same call
% plotVesselGroupMetric.m:10-20 and plotGaussianFitPanels.m:37-43 already made for their own figures.
%
% CORRELATION MATRIX -- [nFrame x nFrame] Pearson correlation between frames, each frame contributing one
% observation per voxel. WHAT is correlated is opts.corrMat.source, and the three choices answer genuinely
% different questions, which is why this is a parameter and not a hardcoded choice:
%   'residual' (DEFAULT) -- the per-frame residual images. This is the FIT-QA reading, and the default
%             because that is this file's whole subject: if the fit tracks the vessel on every frame, the
%             residuals are noise and the matrix is structureless (near-zero off-diagonal, a plain dark
%             field with a bright diagonal). Any BLOCK structure means a set of frames share a residual
%             pattern the model does not explain -- i.e. systematic fit failure over that stretch, which
%             is precisely the failure a per-frame still cannot show you.
%   'patch'   -- the raw frame images. This is the DATA-QA reading (the classic frame-similarity
%             "carpet"): it shows motion, dropouts and artifacts in the input, regardless of the fit.
%             Off-diagonal correlation is high and structured here even for a perfect fit, because a
%             vessel patch mostly looks like itself; read the STRUCTURE, not the level.
%   'model'   -- the fitted model images. Mostly a control: with shape held perVessel, this matrix is a
%             pure function of the fitted (x0,y0) trajectory, so it shows what the fit BELIEVES the
%             motion was, with the data removed.
% WHICH voxels contribute is opts.corrMat.voxels, and 'valid' (the fit's own includeMask minus
% excludeMask -- the region the optimizer actually saw) is the default for ALL three sources, because
% that is the only region where "does the model explain the data" is even a meaningful question: outside
% it nothing was ever modelled, so the "residual" there is just the raw image, and including it answers
% a different question than the one asked.
% What that costs, MEASURED rather than assumed (vessel 1 of the human cache, 60 frames of run 1, mask =
% 26 of 121 voxels) -- 'all' DILUTES the correlation rather than inflating it, for every source:
%     source     valid          all
%     residual   +0.117 mean    +0.099 mean
%     patch      +0.874         +0.714
%     model      +0.973         +0.978
% i.e. the out-of-mask background at this patch size is noise-dominated, not static structure, so it
% adds uncorrelated observations. (An earlier version of this comment asserted the OPPOSITE -- that a
% static background would drive every off-diagonal towards 1 -- which is the intuitive guess and is
% simply wrong here; the numbers above are why it is stated as a measurement now.) So 'all' is a
% legitimate choice for source='patch' if whole-patch similarity is what you want; it is just never the
% SHARPER choice for fit QA, since it dilutes the very region the fit was computed on.
% NaN handling: any voxel not finite in EVERY selected frame is dropped from the observation set (dropped
% once, globally, so every pair of frames is correlated over the SAME voxels -- a pairwise-complete
% variant would make different matrix entries incomparable). A frame with no finite data left renders as
% a transparent row/column (drawFrameCorrMatrix.m's own NaN AlphaData), a visible gap rather than a
% fake zero.
%
% CACHING -- a movie already on disk is NOT re-rendered unless opts.force (Seb's own ask, 2026-09-02;
% opts.force is this pipeline's own name for that switch, set from doIt_human.m's own forceThis exactly
% as drawVessel.m's opts.force already is). The whole point is that re-running the doIt over a populated
% cache costs nothing instead of ~57 min for 36 vessels.
% The target path is fully determined by the vessel's identity + fld/outFld + opts.run + format (see
% resolveMovieFile), so it is resolved BEFORE any work happens. Then:
%   file exists + vessel ALREADY carries a .qaMovie whose .fName is that same file
%       -> FULL skip. Nothing is computed, the stored record is passed straight through. This is the
%          steady state.
%   file exists but there is NO matching record (a movie from an earlier session, or a vessel loaded
%   from a cache predating .qaMovie)
%       -> PARTIAL skip: PASS 1 only, so .qaMovie still comes out fully populated. Pass 1 is the cheap
%          half (~2 ms/frame, so ~1-2 s for a 350-frame run against ~94 s to render). Deliberately not
%          left absent: making the OUTPUT contract depend on whether a cache hit occurred is exactly
%          the sort of conditional field that bites later.
%   otherwise (or opts.force) -> render.
% A record whose .fName points ELSEWHERE (opts.movie.dir changed, different run, different format) is
% treated as absent rather than trusted -- it is stale for this call.
% WHAT THIS DOES NOT CHECK, and it matters: the existing file is never validated against the CURRENT
% opts. Change opts.frameStride, opts.corrMat.source, a colormap, the figure size -- anything that
% alters what the movie LOOKS like -- and a cache hit will happily keep the old file while the rebuilt
% .qaMovie describes the new settings. Use opts.force after changing any rendering option. (Same
% limitation drawVessel.m's own roi cache has: it does not re-derive whether a cached roi was drawn
% under the same opts either.) On a FULL skip this cannot bite, since the stored record is passed
% through untouched alongside the file it actually describes.
% opts.movie.format=='none' writes no file, so there is nothing to cache against and it always runs.
%
% MOVIE FORMATS -- opts.movie.format. The governing constraint is that a QA movie is worthless if you
% cannot actually WATCH it where you work, and Seb works in VSCODE, so "opens in VSCode's own built-in
% media preview, no extension required" is the requirement the default has to satisfy:
%   'mp4' (DEFAULT) -- H.264 in mp4. VSCode's bundled media-preview renders it in an Electron
%             (Chromium) webview, so it plays inline, scrubs frame-by-frame, and is by far the smallest.
%             Measured on a full 350-frame run at the default figPosition: 4.8 MB in 94 s. Against the
%             identical run written as avi -- measured back-to-back at the then-default 1800x1000, so
%             read the RATIO rather than the absolutes -- mp4 was 14.5x SMALLER (7.0 vs 101.7 MB) for
%             1.6x the time (112 vs 70 s). That trade is the whole reason this is the default: the extra
%             time is paid once per movie, whereas being unable to open the file where you actually work
%             is paid every time you look at it.
%   'gif'   -- also displays inline in VSCode (as an animated IMAGE, so no scrub bar, no seeking) and
%             needs no external binary at all. The fallback when ffmpeg is unavailable.
%   'avi'   -- Motion JPEG AVI, straight out of VideoWriter. Kept because it needs nothing external, but
%             VSCode CANNOT display it (Chromium has no MJPEG-in-AVI path), so you would have to open it
%             in an external player -- which is exactly the friction this default exists to remove.
%   'none'  -- render only, write nothing.
%
% WHY mp4 GOES THROUGH ffmpeg RATHER THAN VideoWriter -- MATLAB cannot write H.264 in this install:
% VideoWriter.getProfiles() offers only Archival / Motion JPEG 2000 / Motion JPEG AVI / Grayscale AVI /
% Indexed AVI / Uncompressed AVI, and asking for 'MPEG-4' errors outright
% (MATLAB:audiovideo:VideoWriter:profileNotFound; confirmed R2025a/Linux, 2026-09-01). ffmpeg 8.0.1 with
% libx264 is present and, checked specifically because it is the usual trap, WORKS when shelled out to
% from inside MATLAB despite MATLAB's own rewritten LD_LIBRARY_PATH.
%
% Frames are staged as LOSSLESS PNG and encoded in ONE ffmpeg pass at the end, rather than being pushed
% through VideoWriter's MJPEG first. The staging is cheap in SPACE -- ~0.10 MB/frame, so ~35 MB for a
% 350-frame run, because the figure is mostly flat black and PNG compresses that superbly -- and it
% avoids a first generation of lossy compression on precisely the content that shows it worst: thin gold
% contour lines, small white scatter markers, 10 pt axis text. It is NOT free in time: end to end mp4
% runs ~270 ms/frame against avi's ~200 (measured on a full 350-frame run; a per-element micro-benchmark
% suggested PNG writing was cheaper than MJPEG writeVideo, but that benchmark re-wrote one identical
% frame and did not survive contact with varying real frames -- trust the full-run number).
% The staged folder is removed on success and DELIBERATELY KEPT on an ffmpeg failure: the renders are
% the expensive part, and re-running one ffmpeg command by hand is trivial next to re-rendering.
%
% FRAME STAGING (mp4) -- captured frames go to ONE FLAT RAW rgb24 FILE and are encoded in a single
% ffmpeg pass at the end, rather than as per-frame lossless PNG (changed 2026-09-03). Purely a measured
% decision: writing one frame cost 66.8 ms as PNG against 3.7 ms as a raw fwrite -- 18x -- and PNG
% compression was ~45% of the entire per-frame budget, the largest single item in the loop. Still
% lossless, so no quality is traded; what it costs is transient disk, 4.9 MB/frame (~1.7 GB per
% 350-frame movie) against PNG's 0.07 MB/frame. The staging file is deleted on success and KEPT on an
% ffmpeg failure (the renders are the expensive part; re-running one ffmpeg command is not). It never
% goes under /local -- see opts.movie.stageDir.
% -f rawvideo carries no geometry of its own, so the encode states -s WxH from the SAME refSz every
% frame was pinned to. That makes the settle-and-pin step load-bearing rather than merely defensive: a
% single differently-sized frame would shear the rest of the movie instead of erroring the way
% VideoWriter used to.
%
% PARALLELISM -- opts.parallel (default TRUE), via ensureParPool.m like every other parfor user here.
% The axis follows opts.run, which is Seb's own rule:
%   opts.run = [] or a scalar -> parfor across that vessel's RUNS, one vessel at a time.
%   opts.run = 'cat'          -> one job per vessel, so parfor across VESSELS.
% Deliberately NOT one flat (vessel,run) job list, which would parallelise both at once but would have
% to hand a worker its own copy of a vessel struct for every run of that vessel (~18 MB each here).
% Rendering figures on a worker WORKS -- verified end to end, and fitVesselDiag.m already documents
% being called from inside fitVessel.m's own parfor for the same reason -- but scaling is sub-linear:
% measured 2.13x on 4 workers, and NOT because of CPU saturation (getframe is flat from 1 to 32
% threads, 55.7 vs 50.6 ms), so the ceiling is contention inside the graphics pipeline rather than
% arithmetic. opts.stayOpen is forced off with a warning while parallel: a figure on a worker never
% reaches the client. Progress lines interleave out of order, as parfor output always does.
%
% MASKED TO THE FIT REGION -- every one of the four per-frame panels shows ONLY the voxels the fit
% actually used (its own includeMask minus excludeMask): tiles 1/6 by AlphaData, tiles 11/16 by
% selecting the masked voxels. Every color limit and axis limit is likewise computed over those voxels
% alone, so nothing on screen is drawn on a scale derived from data that is not on screen.
% WHY THIS MATTERS, because an earlier version of this figure got it wrong in both directions -- the
% per-frame motion fit is deliberately restricted to a narrow informative region (doIt_human.m's own
% motion-estimate section: opts.includeMask=fullVesselFromGaussianFit_d1, "to focus on the peak position
% and limit the impact of spurious peaks"), 26 of 121 voxels on vessel 1 of the human cache. A large
% residual OUTSIDE that region is territory the optimizer never saw, not something it got wrong.
% Concretely: vessel 13's whole-patch residual panel showed a bright structured corner reaching +950
% against a ~+-100 noise floor, which reads as gross failure and is nothing of the kind. The first fix
% was to OUTLINE the mask (dashed cyan) while still drawing everything; masking outright (2026-09-02,
% Seb's own ask) is strictly better -- the outline told you where the boundary was but still let
% excluded signal set the color scale and still invited the misreading. The mask is frame-invariant, so
% this costs nothing per frame.
% Consequence worth knowing: this figure no longer shows what SURROUNDS the fit, which is a real thing
% to want when judging whether a fit locked onto the right feature. showVessel.m/fitPatchVesselsDiag.m
% on the same fld are the whole-patch views; this one is deliberately the fit's own-eye view. The
% correlation matrix's opts.corrMat.voxels='valid' default now agrees with the panels rather than
% differing from them.
%
% FIXED SCALES ACROSS FRAMES -- every color limit and every axis limit is computed ONCE over the whole
% selected frame set and held constant for the entire movie. This is not a detail: per-frame autoscaling
% would make a frame's apparent brightness/spread depend on its own extremes, so a fit visibly failing
% on one frame would be silently renormalised to look like every other frame, which would defeat the
% entire purpose. The consequence to be aware of is the opposite one -- a single outlier frame sets the
% scale for all of them, so opts.scale (residual panel) and opts.xlimScale (radial panels) are the knobs
% for when that happens.
%
% SIGNED RADIAL AXIS -- tiles 11/16 use the same convention plotGaussianFitPanels.m's own xd/xe panels
% established, and for the same reason: magnitude is the elliptical-normalized rho scaled by the
% EQUIVALENT-AREA radius sqrt(sx*sy), and the SIGN is taken from each voxel's own raw image x-position
% relative to the fitted x0 (negative = left of center, positive = right). Two-sided, never mirrored --
% every voxel appears exactly once, so a left/right asymmetry (the signature of a mis-positioned fit) is
% visible as such instead of being folded on top of itself. Deliberately NOT
% effectiveDistanceFromFit.m, which is unsigned/major-axis-normalized by design and would fold exactly
% the asymmetry this figure is looking for. The x-axis WINDOW is the image panels' own mm extent times
% opts.xlimScale -- see HORIZONTAL EXTENT above -- so it is fixed across frames rather than breathing
% with the fit (see FIXED SCALES above).
%
% INPUT
%   vessel : ONE scalar vessel struct, OR an array/cell of them (e.g. the whole vesselHuman). A cell is
%            normalized to a struct array, then looped by index -- SAME multi-vessel convention as
%            fitPatchVesselsDiag.m's own vessel input. ONE MOVIE PER VESSEL (there is no multi-vessel
%            grid form: this figure is already a 20-tile layout for a single vessel, and a movie's frame
%            axis is per-vessel anyway since runs/frame counts differ between vessels -- see the
%            per-vessel nFrame column in this file's own PR).
%            MULTI-VESSEL SKIPPING -- a vessel with no image, or one never run through the per-frame fit
%            for this exact fld/outFld, is SKIPPED with a warning identifying it by index/sId/label
%            (fitVesselTimeSeriesDiag:vesselSkipped) rather than aborting the whole batch; same
%            convention as fitPatchVesselsDiag.m. Errors (:noEligibleVessels) only if NONE qualify. A
%            single scalar vessel that doesn't qualify errors directly rather than silently doing
%            nothing.
%   fld    : char, '<fld>.<outFld>' -- e.g. 'tsIm.gaussMotionEstimate' reads
%            vessel.tsIm.gaussMotionEstimate, i.e. the SAME fld fitVesselTimeSeries.m was called with,
%            dot-joined with the SAME opts.outFld it was called with. Folded into one dot-path argument
%            rather than two, same reasoning as fitPatchVesselsDiag.m's own fld: this function only ever
%            READS an existing result, so a bare "outFld" argument would misleadingly suggest it
%            produces one.
%            MUST name a fit with at least one mode='perFrame' parameter -- a fully perVessel/perRun fit
%            has nothing to animate, and is rejected (:notPerFrame) rather than rendered as 350 identical
%            frames.
%   opts   : struct, ALL fields optional (any missing/empty field falls back to its default below) --
%            a FLAT struct except for the two nested groups .corrMat/.movie:
%     .run          [] (DEFAULT) | positive integer | 'cat'
%                     []      -- ONE MOVIE PER RUN. A 6-run vessel gets 6 movies, each named ..._runN,
%                                and .qaMovie comes back as a 6-entry struct array in run order. This is
%                                the default because a run is the unit the scanner actually acquired and
%                                the unit motion correction was estimated over; splicing runs together
%                                invents adjacencies that never existed.
%                     integer -- that one run only.
%                     'cat'   -- every run's frames CONCATENATED into one series/one movie, with run
%                                boundaries marked on the correlation matrix (drawFrameCorrMatrix.m's own
%                                runEdges). This is where cross-RUN structure shows up, which is often
%                                where motion-estimate problems live, but it costs the whole run count in
%                                one job. (Was called 'all' before 2026-09-03; 'all' now ERRORS rather
%                                than being silently accepted, since its meaning would have changed
%                                under anyone's existing call.)
%                   Also decides the PARALLEL AXIS -- see PARALLELISM below.
%     .frames       [] (default, every frame of the selected run(s)) | vector of positive integers
%                   indexing INTO the selected series (so with .run='all', frame numbering is
%                   series-global, not per-run). Bounds-checked (:badFrames).
%     .frameStride  positive integer (default 1) -- keep every Nth frame. Applied AFTER .frames. The
%                   cheap way to preview a long run: stride 10 on a 350-frame run renders 35 frames.
%                   NOTE the correlation matrix is computed on the STRIDED set, so its own axes are the
%                   frames actually shown, never a denser hidden set.
%     .corrMat.source   'residual' (default) | 'patch' | 'model' -- see CORRELATION MATRIX above.
%     .corrMat.voxels   'valid' (default) | 'all' -- see CORRELATION MATRIX above.
%     .corrMat.clim     [lo hi] (default [-1 1], the full symmetric correlation range, FIXED so the same
%                   color means the same correlation across vessels/calls).
%     .corrMat.colormap [] (default: drawFrameCorrMatrix.m's own blue/red diverging map) | [N x 3].
%     .scale        positive scalar (default 10) -- the residual panel's own CLim is +-cLim(2)/.scale,
%                   where cLim(2) is the patch panel's own ceiling; SAME convention and default as
%                   plotGaussianFitPanels.m's own opts.scale (a fixed fraction of peak intensity, not of
%                   the residual's own max, which one outlier voxel would otherwise dominate).
%     .xlimScale    positive scalar (default 1) -- tiles 11/16's shared XLim is the IMAGE panels' own mm
%                   extent times this, so 1 (the default) makes the 1-D panels line up exactly with the
%                   images above them (see HORIZONTAL EXTENT above); <1 zooms into the fit, >1 pulls
%                   back. DELIBERATELY NOT plotGaussianFitPanels.m's own meaning for the same-named
%                   option (there it is a multiple of the fitted equivalent-area RADIUS, default 4):
%                   that convention does not transfer here, because these per-frame fits are routinely
%                   SUB-VOXEL -- vessel 1 of the human cache fits radius=0.147 mm against a 0.4 mm
%                   voxel -- so "4 radii" spans well under a single voxel and would leave ~12 of 121
%                   voxels on the axes (verified visually).
%     .gaussContourPercIntegral  scalar in (0,100) (default 95) -- tile 1's outline is the single contour
%                   enclosing this percent of the fitted Gaussian's own INTEGRATED signal. Converted to
%                   an amplitude fraction by plotVessels.m:3558's own formula, u = 1 - XX/100 (same
%                   option name, same meaning, same open interval as there) -- so 95% of the signal is
%                   the 0.05-of-peak contour. ONE level only, deliberately: this panel is a
%                   frame-by-frame position check, and a multi-level contour stack at this patch size is
%                   ink rather than information.
%     .dotSize      positive scalar (default 6) -- MarkerSize for tiles 11/16's dot markers.
%     .colorbarWidth  scalar in (0,0.5) (default 0.010) -- thickness of tiles 1/6's own thin colorbars,
%                   in normalized figure units (drawThinColorbar.m's own convention).
%     .leftMargin   scalar in [0,0.5) (default 0.075) -- normalized-figure strip reserved to the LEFT of
%                   the whole tiledlayout for tiles 1/6's colorbars and tiles 11/16's y-tick labels.
%                   Needed because TileSpacing is 'none' and drawThinColorbar.m places its bar outside
%                   the axes without shrinking it, so there is otherwise nowhere for it to go.
%     .residualColormap  [] (default: blue/red diverging) | char 'fname'/'fname(args...)' resolved via
%                   resolveColormapFcn.m | numeric [N x 3]. SAME convention as
%                   plotGaussianFitPanels.m's own opts.residualColormap.
%     .movie.format 'mp4' (default) | 'gif' | 'avi' | 'none' -- see MOVIE FORMATS below for the full
%                   comparison and for why mp4 is the default. Short version: mp4 and gif both open in
%                   VSCode's own built-in preview, avi does NOT; mp4 is also ~15x smaller than avi and
%                   scrubbable, so it is the default. 'none' renders and (with .stayOpen) leaves the
%                   figure at its last frame without writing anything.
%     .movie.frameRate  positive scalar (default 10) -- mp4/avi FrameRate, or 'gif' 1/DelayTime.
%     .movie.crf    scalar in [0 51] (default 20) -- x264's own constant-rate-factor, LOWER is better
%                   quality and bigger; 'mp4' ONLY. At the default, a full 350-frame run is 4.8 MB; on a
%                   20-frame clip crf 16/20/30 give 0.41/0.34/0.19 MB. 20 is deliberately conservative
%                   for a QA figure full of thin contour lines and small text, where compression ringing
%                   would be read as data.
%     .movie.quality  scalar in [1 100] (default 75) -- Motion JPEG quality, 'avi' ONLY (ignored by mp4,
%                   whose intermediate is lossless, and by gif). 75/40/15 -> 6.06/4.21/2.81 MB on a
%                   40-frame clip.
%     .movie.stageDir  where the raw frame file is staged, '' (default) = beside the movie itself.
%                   'mp4' ONLY. REFUSED if it resolves under /local (Seb's own constraint): staging is
%                   ~1.7 GB per 350-frame movie times the worker count, and /local here is a shared
%                   7.3 TB volume sitting at 93% full, while /scratch has its own 916 GB. Checked on the
%                   CANONICAL path, so a symlink into /local is caught too.
%     .movie.burnCounter  true|false (default TRUE) -- burn an 'XX/YY' frame counter into the top-left
%                   corner at encode time via ffmpeg's drawtext. 'mp4' ONLY (gif/avi get no ffmpeg pass
%                   and keep updating the MATLAB title per frame instead). XX counts the frames ACTUALLY
%                   RENDERED, so with opts.frameStride>1 it is the position within the rendered subset,
%                   not the raw series index.
%     .movie.fontFile  TrueType path for that counter, '' (default) = probe a short list. If none is
%                   found the counter is skipped WITH A WARNING and everything else still encodes --
%                   losing an overlay is cosmetic, losing the movie is not.
%     .movie.ffmpeg  path to the ffmpeg binary, '' (default) = resolve on the PATH. 'mp4' ONLY. If it
%                   cannot be found, the call ERRORS naming the alternatives rather than silently
%                   writing some other format -- a caller who asked for mp4 because it is the only
%                   format they can VIEW is not helped by quietly getting an avi.
%     .movie.dir    folder for the written movie (default <pwd>/fitVesselTimeSeriesDiagMovies), created
%                   if absent -- same "<pwd>/<functionName>Something" convention fitVesselDiag.m
%                   (fitVesselDiagnostics/) and fitPatchVesselsDiag.m (fitPatchVesselsDiagnostics/) use.
%     .movie.fileName  '' (default, auto: '<sId>_<label>_<fld>_<outFld>_run<R>') | char, extension
%                   supplied by .format and appended if absent.
%     .figPosition  [x y w h] pixels (default [80 80 1800 1000]) -- the figure is captured with getframe,
%                   so THIS is the movie's own resolution. Confirmed to work with an INVISIBLE figure on
%                   a headless session (R2025a/Linux), so this runs in matlab -batch.
%     .stayOpen     true|false (default FALSE) -- opposite of fitPatchVesselsDiag.m's own vessel-mode
%                   default, deliberately: there the figure IS the deliverable, here the movie file is,
%                   and a 36-vessel batch leaving 36 twenty-tile figures open is a memory problem rather
%                   than a convenience. Set true for a single-vessel call you want to poke at
%                   interactively afterwards.
%     .titleStr     override the auto-built per-vessel title (sId/label + fld/outFld + run + frame).
%                   ONLY valid for a single scalar vessel -- errors if set for more than one
%                   (:titleStrMultiVessel), since every vessel needs its own distinct title/filename.
%     .progress     true|false (default TRUE) -- one terminal line per vessel with frame count/timing,
%                   or which kind of cache hit it took.
%     .parallel     true|false (default TRUE) -- render jobs in a parfor. See PARALLELISM below.
%     .nPool        [] (default, the cluster profile's own worker count) | positive integer -- only ever
%                   used when a pool has to be STARTED; ensureParPool.m reuses any existing pool as-is.
%     .dissim.stat  'mean' (default) | 'median' -- how the per-frame 1-r distances are reduced for the
%                   dissimilarity panel. See that panel's own note under LAYOUT.
%     .dissimGap    scalar in [0,0.5) (default 0.012) -- normalized-figure gap between the correlation
%                   matrix and the dissimilarity panel.
%     .corrMat.indicatorColor  RGB (default [1 0 0], RED) -- the moving frame indicator on both the
%                   correlation matrix and the dissimilarity panel. Red because it has to read against
%                   parula, whose middle is the green the indicator used to be.
%     .corrMat.colorbarGap  scalar in [0,0.5) (default 0.045) -- clearance between the matrix and its
%                   own colorbar underneath. Bigger than a side colorbar would need because this panel
%                   KEEPS its x tick labels and they are drawn in exactly that strip.
%     .force        true|false (default FALSE; 0/1 accepted, since doIt's own forceThis is numeric) --
%                   re-render even when the movie file already exists. See CACHING above; in particular
%                   USE THIS after changing any rendering option, because a cache hit does not check the
%                   existing file against the current opts.
%
% OUTPUT
%   vessel : the MUTATED vessel, in the SAME container type and shape it came in (cell stays a cell,
%         struct array stays a struct array) -- this pipeline's "always returns the mutated vessel"
%         convention, the same one drawVessel.m/fitPatchVesselsDiag.m follow. Each DIAGNOSED vessel
%         gains ONE subfield on the very gauss-fit field that was diagnosed, with ONE ENTRY PER
%         RENDERED RUN in run order -- so a struct ARRAY under the default opts.run=[] (6 entries for a
%         6-run vessel) and a scalar for opts.run='cat' or a single run index:
%
%             vessel.<fld>.<outFld>.qaMovie
%               .fName       full path of the written movie ('' when opts.movie.format=='none').
%                            .fName, not .file -- this pipeline's own name for a stored file path
%                            (alignVessel2.m:385, alignVessel3.m:316, promoteRawIm.m:71,
%                            getVesselRoi2.m:185, tsImToBrain.m:64).
%               .format      opts.movie.format actually used.
%               .rendered    true iff THIS call wrote the file; false on either kind of cache hit.
%                            ONE meaning everywhere, including on diagOut. Note that .createdAt and
%                            .opts are NOT reset on a full skip -- they keep describing the file that
%                            is actually on disk, which is what you want when asking "what is this
%                            movie?"; .rendered is the field that answers "did this call do work?".
%               .createdAt   datetime the movie was written.
%               .diagFld     '<fld>.<outFld>' -- what was diagnosed, so the record is self-describing
%                            even when read out of context.
%               .run         the run THIS entry is for -- a run index, or 'cat'.
%               .frames      [1 x nFrame] series indices actually rendered (post .frames/.frameStride).
%               .runOfFrame  [1 x nFrame] which ORIGINAL run each rendered frame came from.
%               .frameInRun  [1 x nFrame] that frame's own index WITHIN its original run.
%               .nFrame      numel(.frames).
%               .corrMat     struct('source','voxels','nFrame','offDiag') -- the matrix SUMMARIZED, not
%                            the matrix: .offDiag has .mean/.median/.p95/.max over its off-diagonal.
%                            .offDiag.mean is the single most useful triage number here (it separated
%                            vessel 13 at +0.298 from the rest at +0.010..+0.18 on the human cache).
%               .opts        the fully resolved opts struct, for provenance.
%         A SKIPPED vessel gains nothing at all -- no empty placeholder -- the same "only processed
%         vessels carry the field" convention fitPatchVessels.m uses for its own opts.outFld.
%         NOTE this reverses an earlier version of this file, which returned only the diagnostics struct
%         and argued in this very block that there was "nothing worth storing back onto the vessel,
%         every number being reconstructable from the stored fit". That was true when the function only
%         drew; it stopped being true once it WRITES A FILE, because a side-effect artifact's location
%         is precisely the thing no amount of recomputation recovers.
%
%   diagOut : struct array, ONE ENTRY PER ELIGIBLE VESSEL (a skipped vessel gets no entry) -- the bulky
%         and the live, deliberately kept OFF the vessel so it never bloats a saved cache:
%     .vesselIdx  original index into the input vessel array/cell (traceability across skips).
%     .idStr      sId+label short identifier.
%     .file       full path of the written movie ('' when .movie.format=='none').
%     .C          the [nFrame x nFrame] correlation matrix ITSELF -- ~1 MB per vessel at 350 frames,
%                 ~35 MB across 36, which is why only its summary goes onto the vessel.
%     .frames / .runOfFrame / .frameInRun  as above.
%     .fitVals    [1 x nFrame] struct array of the per-frame fitted parameter values used (a re-slice
%                 of what is already on the fit).
%     .ax         [1x5] axes handles (tiles 1/6/11/16 then the spanning correlation tile), or gobjects(0)
%                 if the figure was closed (.stayOpen=false).
%     .rendered   as on .qaMovie above.
%     ON A CACHE HIT the heavy fields are absent by construction: a FULL skip computes nothing, so .C
%     and .fitVals are [] and .ax is gobjects(0) (the frame mapping still comes from the stored record);
%     a PARTIAL skip ran pass 1, so .C and .fitVals ARE populated but .ax is still empty (no figure was
%     built). Guard on .rendered rather than assuming .C is there.
%
% NO-ARG CALL -- fitVesselTimeSeriesDiag() (zero input arguments) renders nothing: it prints a reference
% listing of every opts field and its allowed value(s), then returns the default opts struct IN THE
% FIRST OUTPUT -- same convention as fitVesselDiag.m/fitPatchVesselsDiag.m/plotGaussianFitPanels.m's own
% no-arg call (the first output slot carries the opts here and the mutated vessel otherwise, exactly as
% fitPatchVesselsDiag.m's own single output does).
%
%   opts = fitVesselTimeSeriesDiag;
%   opts.frameStride = 10;                            % one preview movie PER RUN (opts.run=[] default)
%   vesselHuman{13} = fitVesselTimeSeriesDiag(vesselHuman{13}, 'tsIm.gaussMotionEstimate', opts);
%   {vesselHuman{13}.tsIm.gaussMotionEstimate.qaMovie.fName}'   % <- one path per run
%
%   opts = fitVesselTimeSeriesDiag;
%   opts.run = 'cat';                                 % every run end-to-end, run boundaries marked
%   opts.corrMat.source = 'patch';                    % data QA instead of fit QA
%   [vesselHuman, diagOut] = fitVesselTimeSeriesDiag(vesselHuman, 'tsIm.gaussMotionEstimate', opts);
%
%   % Triage 36 vessels by how structured their residuals are, without re-rendering anything. NOTE the
%   % max() over .qaMovie: with the default opts.run=[] there is one entry PER RUN, so this asks "what
%   % is this vessel's WORST run" -- a mean across runs would let five clean runs hide one bad one.
%   m = cellfun(@(v) max([v.tsIm.gaussMotionEstimate.qaMovie.corrMat.offDiag.mean]), vesselHuman);
%   [~,ord] = sort(m,'descend');   % worst-looking fits first
%
% See also FITVESSELDIAG, FITPATCHVESSELSDIAG, BUILDGAUSSIANFITDIAG, PLOTGAUSSIANFITPANELS,
% DRAWFRAMECORRMATRIX, FITVESSELTIMESERIES.

    if nargin==0
        printOptsHelp();
        % No-arg call returns the DEFAULT OPTS in the first output -- every sibling's own convention;
        % the output NAME is just the slot, which carries opts here and the mutated vessel otherwise.
        vessel = defaultOpts();
        diagOut = emptyOut();
        return
    end
    assert(nargin>=2, 'fitVesselTimeSeriesDiag:tooFewArgs', ...
        'fld is required -- fitVesselTimeSeriesDiag(vessel, ''<fld>.<outFld>'', opts).');
    if nargin<3 || isempty(opts); opts = struct(); end
    opts = fillOptsDefaults(opts);

    [fldName, outFld] = splitFldPath(fld);

    % READ through a normalized struct array, WRITE back into the caller's own container -- the input
    % shape/type is preserved exactly (cell stays a cell, struct array stays a struct array), because
    % this function now returns the mutated vessel as its first output (see OUTPUT above).
    wasCell = iscell(vessel);
    if wasCell
        vArr = reshape([vessel{:}], size(vessel));
    else
        vArr = vessel;
    end
    assert(isstruct(vArr) && ~isempty(vArr), 'fitVesselTimeSeriesDiag:badVessel', ...
        'vessel must be a vessel struct, or an array/cell of them.');
    nV = numel(vArr);

    assert(isempty(opts.titleStr) || nV==1, 'fitVesselTimeSeriesDiag:titleStrMultiVessel', ...
        ['opts.titleStr is set but %d vessels were passed -- every vessel needs its own distinct ' ...
         'title/filename here (one movie per vessel). Omit it, or call one vessel at a time.'], nV);

    % Eligibility is resolved in a PRE-PASS: the render loops below may be parfor, and a parfor body
    % cannot grow a struct array or issue an ordered per-vessel warning.
    eligible = false(1,nV);
    for i = 1:nV
        vi = vArr(i);
        [ok, why] = isEligible(vi, fldName, outFld);
        if ~ok
            if nV==1
                error('fitVesselTimeSeriesDiag:notEligible', ...
                    'vessel.%s.%s cannot be diagnosed: %s', fldName, outFld, why);
            end
            warning('fitVesselTimeSeriesDiag:vesselSkipped', ...
                'skipping vessel %d (%s): %s', i, idStrOf(vi), why);
            continue
        end
        eligible(i) = true;
    end
    assert(any(eligible), 'fitVesselTimeSeriesDiag:noEligibleVessels', ...
        ['none of the %d vessel(s) passed has a per-frame fit at vessel.%s.%s -- nothing to show ' ...
         '(see the per-vessel warnings above).'], nV, fldName, outFld);

    % --- PARALLEL AXIS (see PARALLELISM in the file header) -----------------
    % Seb's own rule: parfor across RUNS, except for opts.run='cat' where there is only one job per
    % vessel and the axis has to be VESSELS instead. Kept as two explicit loops rather than one flat
    % (vessel,run) job list, which would parallelize both at once but would have to hand every worker
    % its own copy of a vessel struct for EVERY run of that vessel -- ~18 MB per vessel here, so 6x
    % that per 6-run vessel, for no gain over parallelising the runs of one vessel at a time.
    isCat = ischar(opts.run) || isstring(opts.run);
    % Count the jobs FIRST, and only go parallel if there is genuinely more than one. Without this a
    % single vessel + single run still took the parallel path, which cost nothing in time but DID force
    % opts.stayOpen off -- so asking for one interactive figure silently got you none.
    nJobsTotal = 0;
    for i = 1:nV
        if ~eligible(i); continue; end
        if isCat
            nJobsTotal = nJobsTotal + 1;
        else
            imI = vArr(i).(fldName).im; if ~iscell(imI); imI = {imI}; end
            nJobsTotal = nJobsTotal + numel(resolveRunList(opts.run, numel(imI), fldName, outFld));
        end
    end
    useParallel = opts.parallel && nJobsTotal > 1;
    if useParallel; ensureParPool(opts.nPool); end
    % A figure on a worker never reaches the client's desktop, so stayOpen cannot mean anything there
    % -- fitVesselDiag.m documents the same thing about its own save-and-close default. Forced off with
    % a warning rather than silently ignored, since a caller who asked for it would otherwise be left
    % wondering where the figures went.
    if useParallel && opts.stayOpen
        warning('fitVesselTimeSeriesDiag:stayOpenInParallel', ...
            ['opts.stayOpen=true is ignored while rendering in parallel: a figure built on a parfor ' ...
             'worker is never shown on the client. Set opts.parallel=false to keep the figures.']);
        opts.stayOpen = false;
    end

    diagOut = emptyOut();
    if isCat
        % ONE job per vessel -> parallelise over vessels.
        qaCell = cell(1,nV); recCell = cell(1,nV);
        if useParallel
            % opts IS broadcast to every worker (checkcode's PFBNS) and that is the intended trade:
            % opts is tiny next to the ~90 s of rendering each job drives.
            parfor i = 1:nV
                if ~eligible(i); continue; end
                [recCell{i}, qaCell{i}] = doOneJob(vArr(i), i, opts.run, fldName, outFld, opts); %#ok<PFBNS>
            end
        else
            for i = 1:nV
                if ~eligible(i); continue; end
                [recCell{i}, qaCell{i}] = doOneJob(vArr(i), i, opts.run, fldName, outFld, opts);
            end
        end
        for i = 1:nV
            if isempty(qaCell{i}); continue; end
            diagOut(end+1) = recCell{i}; %#ok<AGROW>
            vessel = writeQa(vessel, wasCell, i, fldName, outFld, qaCell{i});
        end
    else
        % ONE job per RUN of each vessel -> parallelise over that vessel's runs, vessel by vessel.
        for i = 1:nV
            if ~eligible(i); continue; end
            vi = vArr(i);
            imI = vi.(fldName).im; if ~iscell(imI); imI = {imI}; end
            runList = resolveRunList(opts.run, numel(imI), fldName, outFld);
            nJ = numel(runList);
            qaJ = cell(1,nJ); recJ = cell(1,nJ);
            % opts and vi ARE broadcast to every worker (checkcode's PFBNS) and that is the intended
            % trade: vi is one vessel (~18 MB) and opts is tiny, against ~94 s of rendering per job.
            if useParallel && nJ>1
                parfor j = 1:nJ
                    [recJ{j}, qaJ{j}] = doOneJob(vi, i, runList{j}, fldName, outFld, opts);
                end
            else
                for j = 1:nJ
                    [recJ{j}, qaJ{j}] = doOneJob(vi, i, runList{j}, fldName, outFld, opts);
                end
            end
            % One .qaMovie ENTRY PER RENDERED RUN, in run order -- a struct array when opts.run=[] gave
            % several. See OUTPUT in the file header.
            qaArr = qaJ{1};
            for j = 2:nJ; qaArr(j) = qaJ{j}; end
            for j = 1:nJ; diagOut(end+1) = recJ{j}; end %#ok<AGROW>
            vessel = writeQa(vessel, wasCell, i, fldName, outFld, qaArr);
        end
    end
end

% ===========================================================================
function vessel = writeQa(vessel, wasCell, i, fldName, outFld, qa)
    % Provenance onto the vessel, as a subfield of the SAME gauss-fit field that was diagnosed --
    % vessel.<fld>.<outFld>.qaMovie (see OUTPUT in the file header). Written per element into the
    % caller's own container, so an INELIGIBLE vessel simply has no .qaMovie rather than an empty
    % placeholder -- the same "only processed vessels carry the field" convention fitPatchVessels.m
    % already establishes for its own opts.outFld.
    if wasCell
        vessel{i}.(fldName).(outFld).qaMovie = qa;
    else
        vessel(i).(fldName).(outFld).qaMovie = qa;
    end
end

% ---------------------------------------------------------------------------
function [rec, qa] = doOneJob(vi, vesselIdx, runSpec, fldName, outFld, opts)
    % ONE (vessel, run-spec) unit of work: cache check, then render or rebuild. Deliberately
    % side-effect-free apart from the movie file itself and its own printing, so it is safe as a parfor
    % body -- everything it needs arrives by argument and everything it produces comes back by return.
    %
    % runSpec overrides opts.run for the whole of this job's downstream: resolveMovieFile,
    % resolveSeries and buildQaProvenance all read opts.run, so setting it once here is what makes
    % "render run 3" work without threading a second run argument through every one of them.
    opts.run = runSpec;

    % --- CACHING (see the CACHING section in the file header) ---------------
    % The movie's path is fully determined by the vessel's identity + fld/outFld + run + format, so it
    % can be resolved BEFORE doing any work. If that file is already there and opts.force is false,
    % this job's expensive half is skipped.
    targetFile = resolveMovieFile(vi, fldName, outFld, opts);
    cacheHit = ~opts.force && ~isempty(targetFile) && exist(targetFile,'file')==2;
    priorQa  = existingQaRecord(vi, fldName, outFld, targetFile);

    if cacheHit && ~isempty(priorQa)
        % FULL skip -- the file is there AND the vessel already carries the record describing it, so
        % there is nothing this call can add. Costs nothing at all, which is the point: this is the
        % steady state when the doIt gets re-run over a populated cache.
        % Passed through VERBATIM except for .rendered. Keeping .createdAt/.opts as stored is the
        % point -- they describe the FILE that is actually on disk, not this call. But .rendered
        % has exactly one meaning everywhere ("THIS call wrote the file"), so it is overridden to
        % false rather than inherited as the true it was written with; otherwise the same field
        % name would answer a different question on qaMovie than it does on diagOut, which a test
        % caught immediately.
        qa  = priorQa;
        qa.rendered = false;
        rec = reusedRec(vi, vesselIdx, priorQa);
        if opts.progress
            fprintf('  fitVesselTimeSeriesDiag: %-24s run %-3s CACHED (record + file)\n', ...
                idStrOf(vi), runSpecStr(runSpec));
        end
    elseif cacheHit
        % PARTIAL skip -- the file is there but the vessel has no record of it (a movie from an
        % earlier session, or a vessel loaded from a cache predating .qaMovie). Re-deriving the
        % record needs only PASS 1, which is the cheap half (~2 ms/frame of buildGaussianFitDiag
        % plus one model evaluation each, ~1-2 s for a 350-frame run against ~94 s to render), so
        % the record is rebuilt rather than left absent. Refusing to populate it here would make
        % the OUTPUT contract depend on whether a cache hit happened, which is exactly the kind of
        % conditional field this pipeline gets bitten by.
        rec = renderOneVessel(vi, vesselIdx, fldName, outFld, opts, false);
        rec.file = targetFile;
        qa  = buildQaProvenance(rec, fldName, outFld, opts, false);
        if opts.progress
            fprintf('  fitVesselTimeSeriesDiag: %-24s run %-3s CACHED (file only, record rebuilt)\n', ...
                idStrOf(vi), runSpecStr(runSpec));
        end
    else
        rec = renderOneVessel(vi, vesselIdx, fldName, outFld, opts, true);
        qa  = buildQaProvenance(rec, fldName, outFld, opts, true);
    end
end

% ---------------------------------------------------------------------------
function runList = resolveRunList(runOpt, nRun, fldName, outFld)
    % WHICH movies one vessel gets, from opts.run (see that option's own note):
    %   []      -> one per run, [1 2 ... nRun]   (the default)
    %   'cat'   -> a single {'cat'} job, every run's frames end to end
    %   scalar  -> that one run
    % Returned as a CELL so the two kinds (numeric run index, the char 'cat') travel in one container.
    if isempty(runOpt)
        runList = num2cell(1:nRun);
    elseif ischar(runOpt) || isstring(runOpt)
        assert(strcmpi(char(runOpt),'cat'), 'fitVesselTimeSeriesDiag:badRun', ...
            'opts.run must be [] (one movie per run), a positive integer, or ''cat'' -- got ''%s''.', char(runOpt));
        runList = {'cat'};
    else
        assert(isnumeric(runOpt) && isscalar(runOpt) && runOpt==round(runOpt) && runOpt>=1 && runOpt<=nRun, ...
            'fitVesselTimeSeriesDiag:badRun', ...
            'opts.run=%s out of range -- vessel.%s.%s was fit over %d run(s) (or use [] / ''cat'').', ...
            mat2str(runOpt), fldName, outFld, nRun);
        runList = {runOpt};
    end
end

% ===========================================================================
function qa = existingQaRecord(v, fldName, outFld, targetFile)
    % The vessel's OWN pre-existing .qaMovie, but only if it actually describes the file we are about to
    % skip -- a record whose .fName points somewhere else (opts.movie.dir changed, a different run, a
    % different format) is stale for this call and is treated as absent rather than trusted. Returns []
    % when there is nothing usable.
    qa = [];
    if isempty(targetFile); return; end
    if ~isfield(v,fldName) || ~isstruct(v.(fldName)); return; end
    if ~isfield(v.(fldName),outFld) || ~isstruct(v.(fldName).(outFld)); return; end
    r = v.(fldName).(outFld);
    if ~isfield(r,'qaMovie') || ~isstruct(r.qaMovie) || isempty(r.qaMovie); return; end
    % SEARCH the array, do not assume a scalar: since opts.run=[] became the default there is one entry
    % PER RUN, and requiring isscalar here meant a multi-run vessel never matched -- every re-run fell
    % through to the "file only, record rebuilt" path and paid pass 1 again instead of skipping outright
    % (caught by the re-run timing being 0.23 s of real work rather than ~0).
    q = r.qaMovie;
    if ~isfield(q,'fName'); return; end
    for k = 1:numel(q)
        if strcmp(char(string(q(k).fName)), targetFile); qa = q(k); return; end
    end
end

% ---------------------------------------------------------------------------
function rec = reusedRec(v, vesselIdx, qa)
    % A diagOut entry for a FULLY skipped vessel. The frame mapping comes from the stored record (it
    % described this same file), while .C/.fitVals/.ax are empty because nothing was computed -- see the
    % file header's own note that diagOut's heavy fields are absent on a cache hit.
    rec = struct('vesselIdx',vesselIdx, 'idStr',idStrOf(v), 'file',qa.fName, 'C',[], ...
                 'frames',qa.frames, 'runOfFrame',qa.runOfFrame, 'frameInRun',qa.frameInRun, ...
                 'fitVals',[], 'ax',gobjects(0), 'rendered',false);
end

% ---------------------------------------------------------------------------
function qa = buildQaProvenance(rec, fldName, outFld, opts, rendered)
    % What gets stored ON the vessel, versus what stays in the second output. The rule: anything a
    % LATER call or a human triaging 36 vessels would want WITHOUT re-running goes here; anything bulky
    % or live stays in diagOut.
    %   IN  -- the movie's own path (the one thing here that is NOT reconstructable from the stored fit:
    %          it is a side-effect artifact location), which frames were rendered and where they came
    %          from, the resolved opts (full provenance, and tiny), and a SUMMARY of the correlation
    %          matrix.
    %   OUT -- the [nFrame x nFrame] matrix itself (350x350 doubles is ~1 MB per vessel, ~35 MB across
    %          36 -- real cache bloat for something fully re-derivable), the per-frame fitVals struct
    %          array (already on the fit, this is just a re-slice of it), and the axes handles (live
    %          graphics have no business in a saved cache).
    % .fName for the path, NOT .file -- that is this pipeline's own field name for a stored file path
    % (alignVessel2.m:385, alignVessel3.m:316, promoteRawIm.m:71, getVesselRoi2.m:185, tsImToBrain.m:64).
    C = rec.C;
    if isempty(C) || size(C,1)<2
        offStats = struct('mean',NaN, 'median',NaN, 'p95',NaN, 'max',NaN);
    else
        off = C(~eye(size(C),'like',true(size(C))));
        off = off(isfinite(off));
        if isempty(off)
            offStats = struct('mean',NaN, 'median',NaN, 'p95',NaN, 'max',NaN);
        else
            offStats = struct('mean',mean(off), 'median',median(off), ...
                              'p95',prctile(off,95), 'max',max(off));
        end
    end

    qa = struct();
    qa.fName       = rec.file;            % '' when opts.movie.format=='none'
    qa.format      = opts.movie.format;
    qa.rendered    = rendered;            % false = the file was already there and this call reused it
    qa.createdAt   = datetime('now');
    qa.diagFld     = sprintf('%s.%s', fldName, outFld);   % what was diagnosed, self-describing
    qa.run         = opts.run;
    qa.frames      = rec.frames;
    qa.runOfFrame  = rec.runOfFrame;
    qa.frameInRun  = rec.frameInRun;
    qa.nFrame      = numel(rec.frames);
    % The correlation matrix, summarized. offDiag.mean is the single most useful triage number this
    % figure produces: on the human cache it separated vessel 13 (+0.298, a genuinely structured
    % in-mask residual) from the rest (+0.010 to +0.18), which is exactly the sort across-vessel
    % comparison you would otherwise have to re-render 36 movies to make.
    qa.corrMat     = struct('source',opts.corrMat.source, 'voxels',opts.corrMat.voxels, ...
                            'nFrame',size(C,1), 'offDiag',offStats);
    qa.opts        = opts;                % full resolved provenance; small
end

% ===========================================================================
function s = emptyOut()
    s = struct('vesselIdx',{}, 'idStr',{}, 'file',{}, 'C',{}, 'frames',{}, 'runOfFrame',{}, ...
               'frameInRun',{}, 'fitVals',{}, 'ax',{}, 'rendered',{});
end

% ---------------------------------------------------------------------------
function [fldName, outFld] = splitFldPath(fld)
    assert((ischar(fld) || (isstring(fld) && isscalar(fld))) && ~isempty(fld), ...
        'fitVesselTimeSeriesDiag:badFld', 'fld must be a nonempty char, ''<fld>.<outFld>''.');
    fld = char(fld);
    parts = strsplit(fld, '.');
    assert(numel(parts)==2 && ~isempty(parts{1}) && ~isempty(parts{2}), ...
        'fitVesselTimeSeriesDiag:badFldPath', ...
        ['fld must be exactly ''<fld>.<outFld>'' (two dot-separated parts) -- got ''%s''. This reads ' ...
         'vessel.<fld>.<outFld>, i.e. the fld fitVesselTimeSeries.m was called with dot-joined with ' ...
         'its own opts.outFld.'], fld);
    fldName = parts{1}; outFld = parts{2};
end

% ---------------------------------------------------------------------------
function s = idStrOf(v)
    sId = ''; lab = '';
    if isfield(v,'sId');   sId = char(string(v.sId));   end
    if isfield(v,'label'); lab = char(string(v.label)); end
    s = strtrim([sId ' ' lab]);
    if isempty(s); s = '(unidentified vessel)'; end
end

% ---------------------------------------------------------------------------
function [ok, why] = isEligible(v, fldName, outFld)
    % Mirrors buildGaussianFitDiag.m's own front-door asserts, but as a PREDICATE -- so a multi-vessel
    % batch can warn-and-skip one heterogeneous member instead of dying (see file header MULTI-VESSEL
    % SKIPPING). Deliberately does NOT try/catch buildGaussianFitDiag itself: that would also swallow
    % genuine bugs inside it as "ineligible vessel".
    ok = false;
    if ~isfield(v,fldName) || ~isstruct(v.(fldName));  why = sprintf('vessel.%s is missing', fldName); return; end
    s = v.(fldName);
    if ~isfield(s,'im') || isempty(s.im);              why = sprintf('vessel.%s has no .im', fldName); return; end
    if ~isfield(s,outFld) || ~isstruct(s.(outFld));    why = sprintf('vessel.%s.%s is missing', fldName, outFld); return; end
    r = s.(outFld);
    if ~isfield(r,'method') || ~any(strcmp(r.method,{'gaussian','parabolic'}))
        why = sprintf('vessel.%s.%s.method is not ''gaussian''/''parabolic''', fldName, outFld); return
    end
    if ~isfield(r,'paramSpec') || ~isfield(r.paramSpec,'mode')
        why = sprintf('vessel.%s.%s has no .paramSpec.mode (not a fitVessel.m-shaped result)', fldName, outFld); return
    end
    if isfield(r,'timeAvg') && r.timeAvg
        why = sprintf('vessel.%s.%s was fit with opts.timeAvg=true (one frame per run) -- nothing to animate', fldName, outFld); return
    end
    if ~any(strcmp(struct2cell(r.paramSpec.mode), 'perFrame'))
        why = sprintf(['vessel.%s.%s has no mode=''perFrame'' parameter (modes: %s) -- every frame ' ...
                       'would be identical, so there is nothing to animate'], fldName, outFld, modeSummary(r.paramSpec.mode));
        return
    end
    ok = true; why = '';
end

% ---------------------------------------------------------------------------
function s = modeSummary(modeStruct)
    fn = fieldnames(modeStruct);
    parts = cell(1,numel(fn));
    for i = 1:numel(fn); parts{i} = sprintf('%s=%s', fn{i}, char(string(modeStruct.(fn{i})))); end
    s = strjoin(parts, ' ');
end

% ===========================================================================
function rec = renderOneVessel(v, vesselIdx, fldName, outFld, opts, renderMovie)
    tStart = tic;
    imS = v.(fldName);
    im  = imS.im; if ~iscell(im); im = {im}; end
    nRun = numel(im);
    Ts   = cellfun(@(x) size(x,4), im);

    % --- which frames -------------------------------------------------------
    [runOfFrame, frameInRun, runEdges] = resolveSeries(opts.run, nRun, Ts, fldName, outFld);
    nSeries = numel(runOfFrame);
    sel = 1:nSeries;
    if ~isempty(opts.frames)
        assert(all(opts.frames>=1 & opts.frames<=nSeries) && all(opts.frames==round(opts.frames)), ...
            'fitVesselTimeSeriesDiag:badFrames', ...
            ['opts.frames must be positive integers within the selected series (1..%d for ' ...
             'opts.run=%s) -- got [%s].'], nSeries, runSpecStr(opts.run), num2str(opts.frames(:).'));
        sel = opts.frames(:).';
    end
    sel = sel(1:opts.frameStride:end);
    assert(~isempty(sel), 'fitVesselTimeSeriesDiag:noFramesSelected', ...
        'opts.frames/.frameStride selected zero frames out of %d.', nSeries);
    nF = numel(sel);
    runSel = runOfFrame(sel); frameSel = frameInRun(sel);

    % --- PASS 1: every per-frame quantity, up front -------------------------
    % Computed for the WHOLE selected set BEFORE any drawing, for two reasons: the fixed color/axis
    % scales (see file header FIXED SCALES) need global extrema, and the correlation matrix needs every
    % frame at once. The frame loop below therefore only DRAWS -- it recomputes nothing (see
    % PERFORMANCE). Memory is trivial at this patch size (an 11x11 patch x 2100 frames x 3 stacks is
    % ~6 MB of doubles), so there is no reason to stream it.
    d1 = buildGaussianFitDiag(v, fldName, outFld, runSel(1), frameSel(1));
    [ny,nx] = size(d1.im);
    voxSz2 = d1.voxSz2; nIso = d1.nIso; validMask = d1.validMask;
    ctrCol = (nx+1)/2; ctrRow = (ny+1)/2;
    [XIso, YIso] = buildIsochromatGrid(nIso, nx, ny, ctrCol, ctrRow, voxSz2);
    toX = @(x) (x-ctrCol)*voxSz2(2);
    toY = @(y) (y-ctrRow)*voxSz2(1);
    [Xpix,Ypix] = meshgrid(1:nx,1:ny);
    Xmm = toX(Xpix); Ymm = toY(Ypix);

    patchStack = zeros(ny,nx,nF);
    modelStack = zeros(ny,nx,nF);
    residStack = zeros(ny,nx,nF);
    otherStack = zeros(ny,nx,nF);   % SECONDARY peaks only (0 for a single-component fit)
    fitValsAll = repmat(d1.fitVals, 1, nF);
    otherFitValsAll = cell(1,nF);
    nComponent = 1;
    for k = 1:nF
        if k==1; dk = d1; else; dk = buildGaussianFitDiag(v, fldName, outFld, runSel(k), frameSel(k)); end
        % dk.modelFull, NOT a local re-evaluation (2026-09-04). This loop used to rebuild the model
        % itself from dk.fitVals, which is the MAIN-PEAK-ONLY view -- so for a multi-component fit
        % every secondary peak landed in the residual panel and the correlation matrix as if the model
        % had failed to explain it. buildGaussianFitDiag.m now returns the full multi-component
        % prediction, and taking it from there also removes the second, independently-maintained model
        % evaluation this file used to carry (the exact duplication that has drifted stale repeatedly
        % elsewhere in this family).
        patchStack(:,:,k) = dk.im;
        modelStack(:,:,k) = dk.modelFull;
        residStack(:,:,k) = dk.im - dk.modelFull;
        otherStack(:,:,k) = dk.otherFull;
        fitValsAll(k)     = dk.fitVals;
        otherFitValsAll{k} = dk.otherFitVals;
        nComponent = max(nComponent, dk.nComponent);
    end
    % Measured signal with the secondary peaks removed -- what the 1-D radial profile is actually
    % about, since that axis is distance from the MAIN peak and a secondary sitting at some unrelated
    % radius otherwise reads as scatter the model failed to follow. Identical in spirit to
    % fitPatchVesselsDiag.m's own "(secondary subtracted)" profile/scatter panels. For a
    % single-component fit otherStack is all zeros, so this is exactly patchStack and nothing changes.
    patchMinusOther = patchStack - otherStack;

    % --- fixed scales (see file header FIXED SCALES) ------------------------
    % ALL of them computed over the FIT-MASK VOXELS ONLY, because that is all any panel now shows (see
    % file header MASKED TO THE FIT REGION). This matters, it is not bookkeeping: driving a color scale
    % or an axis limit off voxels that are not drawn is how vessel 13's unmodelled +950 corner used to
    % set the residual scale for the whole movie while sitting outside every panel's own mask.
    mk  = validMask(:);
    pM  = reshape(patchStack, ny*nx, nF); pM = pM(mk,:);
    mM  = reshape(modelStack, ny*nx, nF); mM = mM(mk,:);
    rM  = reshape(residStack, ny*nx, nF); rM = rM(mk,:);

    cLim = [0, max([max(pM(:)) max(mM(:))])];
    if ~isfinite(cLim(2)) || diff(cLim)==0; cLim(2) = cLim(1)+1; end
    residLim = cLim(2)/opts.scale; if residLim==0 || ~isfinite(residLim); residLim = 1; end
    cLimResid = [-residLim residLim];   % the residual IMAGE stays symmetric: it is a diverging map, and
                                         % a diverging map whose neutral is not at zero is a lie

    % HORIZONTAL EXTENT -- the 1-D panels share the IMAGE panels' own x extent, in mm, so column 1 reads
    % as one stack (Seb's own ask, 2026-09-02: "harmonize the horizontal extent of the 1D plots to that
    % of the 2d images"). flushImLim below is the image's own flush edge, half a voxel outside the
    % outermost voxel CENTER -- the same limit the image panels get, so a point at the edge of the patch
    % sits at the edge of the plot.
    % Checked before adopting rather than assumed, because teff is radius-NORMALIZED and so is not the
    % image's own x coordinate: across vessels 1/5/13/20/30/36 the largest |teff| over MASKED voxels is
    % 1.26-2.00 mm against a 2.20 mm image half-extent, i.e. nothing is clipped. Only the out-of-mask
    % voxels reach 3.1-5.4 mm, and those are no longer drawn.
    flushImLim = [toX(0.5) toX(nx+0.5)];
    xlimShared = flushImLim * opts.xlimScale;

    ylimProf = [0, max([max(pM(:)) max(mM(:))])];
    if ~isfinite(ylimProf(2)) || diff(ylimProf)<=0; ylimProf(2) = ylimProf(1)+1; end
    % ASYMMETRIC residual y-axis (Seb's own ask) -- the residual TRACE gets the data's own range rather
    % than +-max. Deliberately different from the residual IMAGE's symmetric CLim above: a scatter has no
    % diverging colormap to misrepresent, and a real fit's residuals are not symmetric about zero, so
    % forcing symmetry throws away half the vertical resolution. The y=0 reference line drawn in the
    % panel is what keeps the sign readable without the symmetry.
    ylimResid = [min(rM(:)) max(rM(:))];
    if ~all(isfinite(ylimResid)) || diff(ylimResid)<=0; ylimResid = [-1 1]; end
    ylimResid = ylimResid + [-1 1]*0.04*diff(ylimResid);   % a little headroom so points are not clipped
                                                            % by the axes box itself

    % --- correlation matrix -------------------------------------------------
    switch opts.corrMat.source
        case 'residual'; srcStack = residStack;
        case 'patch';    srcStack = patchStack;
        case 'model';    srcStack = modelStack;
        otherwise; error('fitVesselTimeSeriesDiag:badCorrSource', 'unreachable -- validated in fillOptsDefaults.');
    end
    C = frameCorrMatrix(srcStack, validMask, opts.corrMat.voxels);

    % --- EARLY RETURN when the movie is not being (re)written ---------------
    % Everything above this line is PASS 1 plus the correlation matrix -- the cheap half, and the half
    % the .qaMovie record is built from. Everything below builds a figure and encodes frames, which is
    % the ~94 s of a full run. On a cache hit with no stored record the caller asks for exactly the
    % first half (see CACHING in the file header).
    if ~renderMovie
        rec = struct('vesselIdx',vesselIdx, 'idStr',idStrOf(v), 'file','', 'C',C, 'frames',sel, ...
                     'runOfFrame',runSel, 'frameInRun',frameSel, 'fitVals',fitValsAll, ...
                     'ax',gobjects(0), 'rendered',false);
        return
    end

    % --- figure -------------------------------------------------------------
    if isempty(opts.titleStr)
        baseTitle = sprintf('%s (%s.%s, %s, run %s)', idStrOf(v), fldName, outFld, d1.method, runSpecStr(opts.run));
    else
        baseTitle = opts.titleStr;
    end
    [fig, hTitle, ax, hRow, hCol, hDis] = buildQaFigure(baseTitle, opts, C, sel, runEdges, cLim);

    % Frame-INVARIANT bundle, built ONCE. The residual colormap in particular MUST be hoisted: with []
    % (the default) it is a colormap_divergingHue.m call, which searches the sRGB gamut in CIE LCH --
    % rebuilding it per frame cost ~1.7 s/frame, i.e. it dominated the entire render (measured, then
    % fixed; 350 frames would have taken ~10 min).
    % The fitted outline is drawn at ONE level: the contour enclosing opts.gaussContourPercIntegral
    % percent of the Gaussian's own integrated signal (95 by default -- Seb's own ask). The conversion
    % to an amplitude fraction is plotVessels.m:3558's own formula verbatim, u = 1 - XX/100, NOT a second
    % derivation of it: for a 2-D Gaussian the ellipse at Mahalanobis radius r encloses 1-exp(-r^2/2) of
    % the total and sits at amplitude fraction exp(-r^2/2), so those two are the same number and 95%
    % of the signal is the 0.05-of-peak contour (r = 2.448).
    gaussLevelFrac = 1 - opts.gaussContourPercIntegral/100;

    G = struct('voxSz2',voxSz2, 'XIso',XIso, 'YIso',YIso, 'toX',toX, 'toY',toY, 'Xmm',Xmm, 'Ymm',Ymm, ...
               'cLim',cLim, 'cLimResid',cLimResid, 'xlimShared',xlimShared, 'ylimProf',ylimProf, ...
               'ylimResid',ylimResid, 'method',d1.method, 'paramNames',{d1.paramNames}, ...
               'validMask',validMask, 'flushImLim',flushImLim, 'modelFun',d1.modelFun, ...
               'gaussLevelFrac',gaussLevelFrac, 'cbWidth',opts.colorbarWidth, 'dotSize',opts.dotSize, ...
               'tCurve',linspace(xlimShared(1), xlimShared(2), 401), ...
               'residCmap',resolveResidualColormap(opts.residualColormap), ...
               'nComponent',nComponent);
    % .profile is what the 1-D panel plots as "measured": the patch with the SECONDARY peaks removed,
    % identical to .patch for a single-component fit. Kept as its own field rather than subtracted at
    % the point of use so the 2-D panels and the 1-D panel can never disagree about which signal they
    % are showing.
    S = struct('patch',patchStack, 'model',modelStack, 'resid',residStack, ...
               'profile',patchMinusOther, 'otherFitVals',{otherFitValsAll});

    % --- FRAME LOOP: draw only ----------------------------------------------
    % PERFORMANCE -- the four small panels are cla'd and redrawn from the pass-1 stacks (a contour object
    % cannot be updated in place, and at this patch size a redraw is cheaper than the bookkeeping to
    % update every other object individually), while the spanning correlation panel is NEVER touched
    % again: only its crosshair's own constant coordinate is reassigned. That asymmetry is the whole
    % reason drawFrameCorrMatrix.m returns those two handles.
    H = initFramePanels(ax(1:4), S, fitValsAll(1), G);

    % SETTLE BEFORE THE WRITER OPENS -- VideoWriter locks its frame dimensions to the FIRST frame
    % written and then rejects any other size outright (MATLAB:audiovideo:VideoWriter:invalidDimensions,
    % "Frame must be 1800 by 978"). On a freshly created INVISIBLE figure the first getframe can come
    % back short by the title-bar height (1000 -> 978 here) and then correct itself on the next capture,
    % so a writer opened before that settles locks in the wrong size and the whole movie dies partway
    % through -- reproduced deterministically on the 3rd vessel of a batch. Absorbed on a throwaway
    % capture, the same "absorb the first render on a throwaway" fix printFigs.m:41-47 already uses for
    % the analogous headless exportgraphics warm-up, and with the same FULL drawnow (not limitrate) that
    % file uses before each capture.
    refSz = settleFigureForCapture(fig);
    [writer, movieFile] = openMovie(v, fldName, outFld, opts);
    nResized = 0;

    % THE TITLE IS STATIC when the counter is burned in by ffmpeg (mp4). Re-setting the annotation's
    % String dirties the whole figure and cost a MEASURED 32 ms/frame -- 44% of a drawnow+getframe
    % cycle (41.5 ms frozen vs 73.7 ms updated) -- so for mp4 the per-frame number moves to ffmpeg's
    % own drawtext (see encodeMp4) and this text never changes again. gif/avi get no ffmpeg pass, so
    % there the per-frame update is kept: correct output beats a speedup that only applies to a
    % non-default format.
    liveTitle = ~(strcmp(opts.movie.format,'mp4') && opts.movie.burnCounter);
    set(hTitle, 'String', sprintf('%s -- %d frame(s) of %d', baseTitle, nF, nSeries));
    for k = 1:nF
        H = updateFramePanels(H, k, S, fitValsAll(k), G);
        set(hRow, 'YData', [sel(k) sel(k)]);
        set(hCol, 'XData', [sel(k) sel(k)]);
        set(hDis, 'YData', [sel(k) sel(k)]);   % same frame, echoed on the dissimilarity panel
        if liveTitle
            set(hTitle, 'String', sprintf('%s -- frame %d/%d (run %d, frame %d of that run)', baseTitle, ...
                sel(k), nSeries, runSel(k), frameSel(k)));
        end
        drawnow limitrate
        nResized = nResized + writeMovieFrame(writer, fig, opts.movie.format, k==1, refSz);
    end
    closeMovie(writer, opts.movie.format, refSz, nF);
    % Reported, not swallowed: if the warm-up above did not fully settle the render, say how often, so a
    % geometrically-wrong movie is never produced quietly.
    if nResized > 0
        warning('fitVesselTimeSeriesDiag:frameResized', ...
            ['%d of %d captured frame(s) came back at a different pixel size than the first and were ' ...
             'resized to [%d %d] to keep the movie writable -- the figure render did not settle. The ' ...
             'movie is valid but those frames are resampled; re-run, or set opts.figPosition to a size ' ...
             'that fits the display without clamping.'], nResized, nF, refSz(1), refSz(2));
    end

    if opts.progress
        fprintf('  fitVesselTimeSeriesDiag: %-28s %4d frame(s) of %4d, corr=%s/%s, %5.1f s%s\n', ...
            idStrOf(v), nF, nSeries, opts.corrMat.source, opts.corrMat.voxels, toc(tStart), ...
            iff(isempty(movieFile), '', sprintf(' -> %s', movieFile)));
    end

    rec = struct('vesselIdx',vesselIdx, 'idStr',idStrOf(v), 'file',movieFile, 'C',C, 'frames',sel, ...
                 'runOfFrame',runSel, 'frameInRun',frameSel, 'fitVals',fitValsAll, 'ax',ax, ...
                 'rendered',renderMovie);
    if ~opts.stayOpen
        close(fig);
        rec.ax = gobjects(0);
    end
end

% ---------------------------------------------------------------------------
function s = runSpecStr(runOpt)
    if ischar(runOpt) || isstring(runOpt); s = char(runOpt); else; s = num2str(runOpt); end
end

% ---------------------------------------------------------------------------
function v = iff(c, a, b)
    if c; v = a; else; v = b; end
end

% ---------------------------------------------------------------------------
function [runOfFrame, frameInRun, runEdges] = resolveSeries(runOpt, nRun, Ts, fldName, outFld)
    % Maps a SERIES position (what opts.frames/.frameStride and the correlation matrix axis both index)
    % onto (original run, frame within that run) -- which is the addressing buildGaussianFitDiag.m
    % itself takes. For a single run this is the identity; for 'all' it is runs end-to-end in order,
    % which is the SAME order fitVessel.m's own combineRuns uses for opts.runCat (cat(4, imCell{:}),
    % runIdx = repelem(1:nPseudo, Ts)) -- so a series position here means the same thing it would have
    % meant to the fitter.
    if ischar(runOpt) || isstring(runOpt)
        assert(strcmpi(char(runOpt),'cat'), 'fitVesselTimeSeriesDiag:badRun', ...
            'opts.run must be [] , a positive integer, or ''cat'' -- got ''%s''.', char(runOpt));
        runOfFrame = repelem(1:nRun, Ts);
        frameInRun = cell2mat(arrayfun(@(t) 1:t, Ts, 'UniformOutput',false));
        runEdges   = cumsum(Ts(1:end-1)) + 0.5;   % between the last frame of run r and the first of r+1
    else
        assert(isnumeric(runOpt) && isscalar(runOpt) && runOpt==round(runOpt) && runOpt>=1 && runOpt<=nRun, ...
            'fitVesselTimeSeriesDiag:badRun', ...
            'opts.run=%s out of range -- vessel.%s.%s was fit over %d run(s) (or use [] / ''cat'').', ...
            mat2str(runOpt), fldName, outFld, nRun);
        runOfFrame = repmat(runOpt, 1, Ts(runOpt));
        frameInRun = 1:Ts(runOpt);
        runEdges   = [];
    end
end

% ---------------------------------------------------------------------------
function mv = toModelVals(fitVals, paramNames, method)
    % radius/aspectRatio -> literal sx,sy for gaussianModel.m, exactly as buildGaussianFitDiag.m and
    % plotGaussianFitPanels.m both do (this project's fits have no .sx/.sy field; see fitVessel.m's own
    % PARAMETERIZATION note). Kept as an unnamed local, never stored.
    %
    % Takes the three pieces POSITIONALLY rather than a diag struct: the frame loop only carries
    % fitVals per frame (paramNames/method are frame-invariant), and passing a cell field through
    % struct() to fake a diag would silently build a struct ARRAY, one element per paramName.
    vals = cellfun(@(n) fitVals.(n), paramNames, 'UniformOutput',false);
    if strcmpi(method,'gaussian')
        sqrtAR = sqrt(vals{5});
        mv = [vals(1:3), {vals{4}*sqrtAR, vals{4}/sqrtAR}, vals(6:7)];
    else
        mv = vals;
    end
end

% ---------------------------------------------------------------------------
function teffSigned = signedTeff(fitVals, paramNames, method, Xmm, Ymm, voxSz2)
    % SIGNED RADIAL AXIS (see file header) -- magnitude is the elliptical-normalized rho scaled by the
    % EQUIVALENT-AREA radius sqrt(s1*s2); sign is each voxel's own raw image x-position relative to the
    % fitted x0 (negative = left of center). Two-sided and never mirrored, so a left/right asymmetry --
    % the signature of a mis-positioned fit, and the thing this figure is looking for -- survives.
    % ONE home for this formula, called both by the global-window pass and by the frame loop, so the
    % window can never be computed on a different axis from the one actually drawn.
    mv = toModelVals(fitVals, paramNames, method);
    s1 = mv{4}; s2 = mv{5};
    if s1==0 || ~isfinite(s1); s1 = mean(voxSz2); end
    if s2==0 || ~isfinite(s2); s2 = mean(voxSz2); end
    rEq = sqrt(abs(s1*s2)); if rEq==0 || ~isfinite(rEq); rEq = mean(voxSz2); end
    th = fitVals.theta;
    dX = Xmm - fitVals.x0; dY = Ymm - fitVals.y0;
    xr =  dX*cos(th) + dY*sin(th);
    yr = -dX*sin(th) + dY*cos(th);
    rho = sqrt((xr/s1).^2 + (yr/s2).^2);
    sideSign = ones(size(Xmm)); sideSign(dX<0) = -1;
    teffSigned = sideSign .* rho .* rEq;
end

% ---------------------------------------------------------------------------
function C = frameCorrMatrix(stack, validMask, voxelMode)
    % [nFrame x nFrame] Pearson correlation, frames as VARIABLES and voxels as OBSERVATIONS (so
    % corrcoef's own column convention applies directly to an [nVox x nFrame] reshape). See the file
    % header CORRELATION MATRIX for the source/voxel choices and their consequences.
    [ny,nx,nF] = size(stack);
    X = reshape(stack, ny*nx, nF);
    switch voxelMode
        case 'valid'; X = X(validMask(:), :);
        case 'all'    % every voxel -- no restriction
        otherwise; error('fitVesselTimeSeriesDiag:badCorrVoxels', 'unreachable -- validated in fillOptsDefaults.');
    end
    % Drop non-finite voxels ONCE, globally, so every matrix entry is computed over the SAME observation
    % set (a pairwise-complete variant would make entries mutually incomparable -- see file header).
    keep = all(isfinite(X), 2);
    assert(any(keep), 'fitVesselTimeSeriesDiag:noFiniteVoxels', ...
        ['no voxel is finite across all %d selected frames (voxel set ''%s'', %d voxel(s)) -- the ' ...
         'correlation matrix would be entirely undefined.'], nF, voxelMode, size(X,1));
    X = X(keep, :);
    if size(X,1) < 3
        warning('fitVesselTimeSeriesDiag:fewCorrVoxels', ...
            ['only %d voxel(s) contribute to the frame-by-frame correlation matrix (voxel set ''%s'') ' ...
             '-- a correlation over so few observations is extremely noisy; consider ' ...
             'opts.corrMat.voxels=''all'' if this is a patch-source matrix.'], size(X,1), voxelMode);
    end
    if size(X,1) < 2
        C = NaN(nF);                  % corrcoef would return a scalar-ish degenerate result here
        return
    end
    C = corrcoef(X);                  % [nF x nF]
    if isscalar(C) && nF>1; C = NaN(nF); end
end

% ===========================================================================
function dis = frameDissimilarity(C, stat)
    % Per-frame DISTANCE to the rest of the series: the mean (or median) of 1-rho across that frame's
    % own row of the correlation matrix, diagonal excluded.
    %
    % 1-rho is Seb's OWN correlation->distance conversion, lifted from bassReg2/QAdendrogram.m:162
    %     linkage(squareform(1-rho), 'average')   % Convert correlation to distance
    % applied there to the very matrix bassReg2/xCorrQA.m:142 builds (cross-frame corr over masked
    % voxels) -- the same object this figure draws. Not a metric invented here.
    % Reducing it per frame with a MEAN is what matches that call's 'average' linkage; 'median' is
    % offered because a series with a few wild frames has every OTHER frame's mean dragged up by the
    % distances to them.
    % The diagonal is removed with NaN + omitnan rather than by index gymnastics, so a matrix that
    % already carries NaN rows (a frame with no finite voxels -- see frameCorrMatrix) stays NaN here
    % instead of quietly averaging in a zero.
    n = size(C,1);
    if isempty(C) || n<2; dis = nan(1, max(n,0)); return; end
    D = 1 - double(C);
    D(1:n+1:end) = NaN;                 % 1-rho is 0 on the diagonal by construction; drop it
    switch stat
        case 'mean';   dis = mean(D, 2, 'omitnan').';
        case 'median'; dis = median(D, 2, 'omitnan').';
        otherwise; error('fitVesselTimeSeriesDiag:badDissimStat', 'unreachable -- validated in fillOptsDefaults.');
    end
end

% ---------------------------------------------------------------------------
function [fig, hTitle, ax, hRow, hCol, hDis] = buildQaFigure(baseTitle, opts, C, sel, runEdges, cLim)
    % Column 1 holds the four per-frame panels stacked with NO gap; one panel spanning the remaining
    % 4x4 block of the notional 5x4 grid holds the correlation matrix. See file header LAYOUT.
    figArgs = {'Color','k', 'Name',['fitVesselTimeSeriesDiag: ' baseTitle], ...
               'Position',opts.figPosition, 'Visible', iff(opts.stayOpen,'on','off')};
    fig = figure(figArgs{:});
    try theme(fig,'dark'); catch; end
    % MANUAL POSITIONING, NOT tiledlayout -- and this is a reversal worth recording, since the layout
    % originally used tiledlayout(4,5) following plotIRFmat.m's own spanning-tile idiom. tiledlayout
    % cannot reserve the LEFT STRIP this figure needs for tiles 1/6's own colorbars: setting
    % t.OuterPosition (or t.InnerPosition) changes the property but does NOT move the child axes --
    % verified directly, ax(1).Position stayed at x=0.013 after t.OuterPosition=[0.12 0 0.88 1] -- so
    % the colorbars landed against the figure edge with their tick labels and 'meas'/'meas-pred' labels
    % clipped off-figure. Padding 'loose' would make room but by adding margin everywhere, which is the
    % opposite of the ask. Explicit positions also make "minimize the space between all graphics
    % elements" something this file states rather than something it hopes a layout manager infers.
    %
    % GEOMETRY -- every cell is SQUARE and the figure aspect is solved to make that exact, so there is
    % no dead space to trim (see squarePlotBox for why square). With column-1 cell side w, the
    % correlation block must be 4w wide to be square at 4w tall, so
    %     mL + w + gapX + 4w + mR = 1   ->   w = (1 - mL - mR - gapX)/5
    % and the 4 stacked cells must fill the vertical span, which fixes the figure's own aspect:
    %     4 * (w * Wpx/Hpx) = 1 - mB - mT
    % opts.figPosition's default is that solution; change the margins and it stops being exact (the
    % panels stay square via pbaspect, you just get slack).
    mL = opts.leftMargin; mR = opts.rightMargin; mB = opts.bottomMargin; mT = opts.topMargin;
    gapX = opts.colGap; gapD = opts.dissimGap;
    % SIX notional columns now (2026-09-03): column 1's four stacked panels, the 4-wide correlation
    % block, and the 1-wide dissimilarity panel to its right.
    %     mL + w + gapX + 4w + gapD + w + mR = 1   ->   w = (1 - mL - mR - gapX - gapD)/6
    w  = (1 - mL - mR - gapX - gapD)/6;
    figAR = opts.figPosition(3)/opts.figPosition(4);
    h  = w * figAR;
    ax = gobjects(1,6);
    for i = 1:4
        ax(i) = axes('Parent',fig, 'Units','normalized', ...
                     'Position',[mL, mB + (4-i)*h, w, h], 'Color','k');
    end
    xCorrX = mL + w + gapX;
    ax(5) = axes('Parent',fig, 'Units','normalized', ...
                 'Position',[xCorrX, mB, 4*w, 4*h], 'Color','k');
    % ax(6): the dissimilarity panel. SAME vertical extent as the correlation block (mB..mB+4h) so the
    % two share a frame axis exactly -- that alignment is the whole point, since this panel is read
    % across from the matrix row by row.
    ax(6) = axes('Parent',fig, 'Units','normalized', ...
                 'Position',[xCorrX + 4*w + gapD, mB, w, 4*h], 'Color','k');
    % The running frame identifier, as a figure-level annotation (there is no layout title to hang it
    % on any more). Kept despite "remove the title" -- that ask was about the per-PLOT titles, and a
    % movie with no frame counter is materially harder to use; it costs one thin strip at the top.
    hTitle = annotation(fig, 'textbox', [0, 1-mT, 1, mT], 'String','', ...
        'Color','w', 'EdgeColor','none', 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontSize',12, 'Interpreter','none');

    % The correlation panel: drawn ONCE, never redrawn (see the frame loop's own PERFORMANCE note).
    % Its axis is the SELECTED series positions, so with opts.frameStride>1 the ticks are the frames
    % actually shown -- never a denser hidden set.
    cmapCorr = opts.corrMat.colormap;
    if isempty(cmapCorr); cmapCorr = parula(256); end   % Seb's own ask; drawFrameCorrMatrix.m keeps its
                                                         % own diverging default for other callers
    % RED indicator (Seb's own ask) -- reads against parula's blue/green/yellow ramp far better than
    % the green it replaced, which collided with the middle of the map.
    [~, hRow, hCol] = drawFrameCorrMatrix(ax(5), C, sel, opts.corrMat.clim, cmapCorr, ...
        opts.corrMat.indicatorColor, 1.25, runEdges);
    % Everything off except the x ruler (Seb's own ask): the frame axis only needs to be read once, and
    % the matrix is symmetric so a y ruler duplicating it is pure clutter. No labels, no title -- the
    % colorbar label names the quantity.
    set(ax(5), 'FontSize',11, 'YTick',[], 'YTickLabel',[]);
    xlabel(ax(5), ''); ylabel(ax(5), '');
    % drawThinColorbar.m here too, NOT colorbar(ax,'eastoutside') -- a native colorbar SHRINKS the axes
    % it belongs to, which took this panel from its allotted 4w wide down to 0.6017 (measured) and so
    % broke both its squareness and the exact-fit geometry the figure aspect is solved for, leaving a
    % strip of dead space at the right. Same reason tiles 1/6 use it (see buildQaFigure's own note).
    % BOTTOM (Seb's own ask) -- the right-hand side is now occupied by the dissimilarity panel, and a
    % scale bar under the matrix reads naturally against its x ruler anyway.
    % A big gap here, not the default 0.01: this panel KEEPS its x tick labels (Seb's own ask) and they
    % are drawn in exactly the strip below the axes that the bar wants.
    cb = drawThinColorbar(ax(5), cmapCorr, opts.corrMat.clim, opts.colorbarWidth, 'bottom', ...
                          opts.corrMat.colorbarGap);
    cb.Color = 'w'; cb.Label.String = 'frame xcorr'; cb.Label.Color = 'w';

    % --- ax(6): per-frame DISSIMILARITY, axes swapped ------------------------
    % 1 - rho is Seb's OWN established correlation->distance conversion, not a choice made here:
    % bassReg2/QAdendrogram.m:162 builds its frame clustering as
    %     linkage(squareform(1-rho), 'average')   % Convert correlation to distance
    % on the very matrix bassReg2/xCorrQA.m:142 produces (corr over masked voxels, cross-frame), which
    % is the same matrix this figure draws. That old pair is also where this panel's ORIENTATION comes
    % from: QAdendrogram drew its dendrogram with 'Orientation','right', i.e. immediately right of the
    % correlation matrix with frames running vertically -- exactly the arrangement asked for here.
    % What is plotted is the per-frame MEAN of that distance to every other frame (diagonal excluded),
    % i.e. "how unlike the rest of the series is this frame" -- the per-frame reduction of the same
    % quantity the clustering consumed pairwise. A spike is a candidate frame to censor, which is what
    % the old pipeline used the clustering for (it wrote censor files from the main cluster).
    dis = frameDissimilarity(C, opts.dissim.stat);
    hold(ax(6),'on'); set(ax(6), 'Color','k', 'XColor','w', 'YColor','w');
    plot(ax(6), dis, sel, '-', 'Color',[1 1 1], 'LineWidth',1.0);
    set(ax(6), 'YDir','reverse', 'YLim',get(ax(5),'YLim'), 'YTick',[], 'YTickLabel',[], 'FontSize',11);
    % YDir reverse + the matrix's own YLim: frame N sits at the same height in both panels, which is
    % the only reason putting them side by side buys anything.
    dLim = [0, max([max(dis(isfinite(dis))) eps])];
    if ~all(isfinite(dLim)) || diff(dLim)<=0; dLim = [0 1]; end
    set(ax(6), 'XLim', dLim*1.05);
    grid(ax(6),'on');
    % Three ticks only: the panel is one column wide, and the default tick density ran the labels into
    % each other and off the figure's right edge (verified visually).
    set(ax(6), 'XTick', linspace(0, dLim(2), 3));
    xtickformat(ax(6), '%.2f');
    xlabel(ax(6), sprintf('1-r (%s)', opts.dissim.stat), 'Color','w', 'FontSize',10);
    % The moving frame indicator is echoed here in the SAME red, so the eye carries one horizontal line
    % straight across from the matrix into this trace. Seb asked for the red line on the matrix; adding
    % the matching one here is the reason the panels are adjacent at all.
    hDis = plot(ax(6), dLim*1.05, [sel(1) sel(1)], '-', 'Color',opts.corrMat.indicatorColor, 'LineWidth',1.25);
    set(ax(6), 'XLim', dLim*1.05);

    % The patch panel's own CLim is fixed here (not per frame) -- see file header FIXED SCALES.
    set(ax(1), 'CLim', cLim);
end

% ---------------------------------------------------------------------------
function H = initFramePanels(ax, S, fitVals, G)
    % Creates every per-frame graphics object ONCE, at selected-frame 1, THROUGH THE SHARED PRIMITIVES
    % (plotVesselPatch.m/drawModelImage.m/drawRadialProfileMeasuredPredicted.m/
    % drawGaussianPeakContour.m) -- so all panel styling and coordinate math keeps exactly one home,
    % the same one fitVesselDiag.m/plotGaussianFitPanels.m use. updateFramePanels below then only ever
    % reassigns DATA properties on the handles returned here; it creates nothing and styles nothing.
    %
    % WHY the split (measured, not guessed -- see this file's own PR): a cla+full-redraw loop costs
    % ~690 ms/frame, and essentially all of it is the image/axes machinery rather than the data --
    % plotVesselPatch 46 ms, drawModelImage 56 ms, axis('image') 37 ms per call, because each rebuilds
    % an image object, reapplies a 256-entry per-axes colormap and re-runs a tiledlayout aspect
    % negotiation. Reassigning CData on an existing image is ~0.1 ms. Every axis limit, colormap, label
    % and title in this figure is frame-INVARIANT by design (see file header FIXED SCALES), so there was
    % never anything for a redraw to recompute in the first place.
    nx = size(S.patch,2); ny = size(S.patch,1);
    gold  = [1 0.84 0.1];   % SAME main-peak color fitPatchVesselsDiag.m's own panel 1 uses; also the
                             % color of the predicted CURVE on tile 3, so the outline and the model
                             % trace are visibly the same object seen two ways (Seb's own ask)
    red   = [1 0.25 0.25];  % SECONDARY peaks -- again the same convention as fitPatchVesselsDiag.m
    white = [1 1 1];
    mk    = G.validMask;
    H = struct();

    % --- tile 1: masked patch + the 95%-of-signal fitted outline ------------
    % MASKED, not outlined (Seb's own ask, 2026-09-02, replacing the dashed-cyan validMask outline this
    % panel used to carry): every voxel the fit did not use is made fully transparent via AlphaData, so
    % the panel shows exactly the data the fit saw and nothing else. This is strictly better than the
    % outline it replaces -- the outline told you where the boundary was but still drew the out-of-mask
    % signal, which is how vessel 13's unmodelled bright corner read as a gross fit failure when it was
    % simply excluded territory.
    plotVesselPatch(ax(1), S.patch(:,:,1), G.voxSz2);
    set(ax(1), 'CLim', G.cLim);       % re-pin: plotVesselPatch sets its own per-image CLim
    H.patchIm = findobj(ax(1), 'Type','image');
    set(H.patchIm, 'AlphaData', double(mk));
    H.contour = gobjects(0);          % created per frame -- see updateFramePanels
    flushImageAxes(ax(1), G.flushImLim);
    H.cbMeas = drawThinColorbar(ax(1), gray(256), G.cLim, G.cbWidth, 'left');
    H.cbMeas.Label.String = 'meas'; H.cbMeas.Label.Color = 'w';

    % --- tile 2 (row 2): masked residual image ------------------------------
    % Same masking as tile 1. NOTE Seb's own list asked for masking explicitly on tiles 1/3/4 and was
    % silent for this one; masked here too, deliberately, because the alternative is a panel that
    % disagrees with the one directly above it about which voxels exist, and because G.cLimResid is
    % computed over masked voxels only -- an unmasked panel would be drawn on a scale derived from data
    % it is showing outside of. Flagged rather than assumed: one AlphaData line to revert.
    drawModelImage(ax(2), G.toX(1:nx), G.toY(1:ny), S.resid(:,:,1), G.cLimResid, G.residCmap);
    axis(ax(2), 'image');             % drawModelImage does not aspect-lock; tile 1 does, and the two
                                       % must read as the same patch
    H.residIm = findobj(ax(2), 'Type','image');
    set(H.residIm, 'AlphaData', double(mk));
    flushImageAxes(ax(2), G.flushImLim);
    H.cbResid = drawThinColorbar(ax(2), G.residCmap, G.cLimResid, G.cbWidth, 'left');
    H.cbResid.Label.String = 'meas-pred'; H.cbResid.Label.Color = 'w';

    teffSigned = signedTeff(fitVals, G.paramNames, G.method, G.Xmm, G.Ymm, G.voxSz2);
    measured = S.profile(:,:,1); residual = S.resid(:,:,1);   % .profile = secondary-subtracted

    % --- tile 3: measured dots + CONTINUOUS predicted curve -----------------
    % Plain small dots via plot('.') rather than drawRadialProfileMeasuredPredicted.m: that primitive
    % draws size-40 filled circles with a black edge for measured and open circles for predicted, and
    % neither survives Seb's own "small dot markers" / "predicted as a continuous function" spec -- at
    % small sizes its black marker edge is most of the marker. The primitive is still the right thing
    % for fitVesselDiag.m/plotGaussianFitPanels.m, which want that heavier styling; this panel just is
    % not one of its cases.
    hold(ax(3),'on'); set(ax(3), 'Color','k', 'XColor','w', 'YColor','w');
    H.measDot = plot(ax(3), teffSigned(mk), measured(mk), '.', 'Color',white, 'MarkerSize',G.dotSize);
    % The predicted profile is the fit's OWN 1-D-mode formula evaluated on a tight grid across this
    % panel's x range -- gaussianModel.m's 1D mode, the same formula the optimizer used, not a
    % re-derivation. Its own X argument is rho*sMajor, so this panel's radiusEq-scaled axis is converted
    % on the way in; algebraically identical, just a different-but-proportional distance axis (the SAME
    % conversion plotGaussianFitPanels.m's own xe panel documents).
    H.predCurve = plot(ax(3), G.tCurve, zeros(size(G.tCurve)), '-', 'Color',gold, 'LineWidth',1.25);
    xlim(ax(3), G.xlimShared); ylim(ax(3), G.ylimProf); grid(ax(3),'on');
    set(ax(3), 'FontSize',10, 'XTick',[]);   % x-axis removed -- tile 4 below carries the shared one
    % Label says so when the dots are secondary-subtracted -- a panel that silently shows a DIFFERENT
    % measured signal from the patch panel above it would be misread as disagreement between them.
    if G.nComponent > 1
        ylabel(ax(3), 'meas-sec,model', 'Color','w', 'FontSize',10);
    else
        ylabel(ax(3), 'meas,model', 'Color','w', 'FontSize',10);
    end
    squarePlotBox(ax(3));

    % --- tile 4: 1-D residual, masked, small white dots ---------------------
    hold(ax(4),'on'); set(ax(4), 'Color','k', 'XColor','w', 'YColor','w');
    % y=0 FIRST so the dots sit on top of it, and because the y-axis here is asymmetric (see the
    % ylimResid note at the call site) the line is the only thing that shows where zero is.
    plot(ax(4), G.xlimShared, [0 0], '-', 'Color',[0.6 0.6 0.6], 'LineWidth',0.75);
    H.residDot = plot(ax(4), teffSigned(mk), residual(mk), '.', 'Color',white, 'MarkerSize',G.dotSize);
    xlim(ax(4), G.xlimShared); ylim(ax(4), G.ylimResid); grid(ax(4),'on');
    set(ax(4), 'FontSize',10);
    xlabel(ax(4), 'signed distance from fitted center (mm)', 'Color','w', 'FontSize',10);
    ylabel(ax(4), 'meas-pred', 'Color','w', 'FontSize',10);
    squarePlotBox(ax(4));

    H.ax = ax; H.gold = gold; H.red = red; H.ny = ny; H.nx = nx; H.mk = mk;
end

% ---------------------------------------------------------------------------
function flushImageAxes(ax, flushImLim)
    % XLim/YLim flush to the image's own edges (Seb's own ask) -- half a voxel outside the outermost
    % voxel CENTER, which is where imagesc actually draws the outer edge of the outer pixel. axis
    % 'image' alone does NOT do this: it fixes the ASPECT and then leaves whatever limits the data
    % implied, which on a patch this small leaves a visible margin of axes background around the image.
    % Ticks off: these panels carry a colorbar for their scale and their extent is the patch itself.
    set(ax, 'XLim',flushImLim, 'YLim',flushImLim, 'XTick',[], 'YTick',[], ...
            'Box','off', 'FontSize',10);
    pbaspect(ax, [1 1 1]);   % see squarePlotBox below for why this is forced rather than left to
                              % axis('image')
end

% ---------------------------------------------------------------------------
function squarePlotBox(ax)
    % Force a SQUARE plot box. Applied to all four column-1 panels so their boxes are identical in size
    % and x-position, which is what actually delivers "harmonize the horizontal extent of the 1D plots to
    % that of the 2d images" (Seb's own ask) -- matching the XLim alone does not, because the two panel
    % kinds fill their tile differently: axis('image') shrinks an image panel inside its tile to honour
    % the data aspect, while a plot panel fills the tile outright. Measured before this was added: image
    % panels came out 0.1753 wide at x=0.0424 against 0.1819 at x=0.0391 for the 1-D panels -- visibly
    % misaligned edges in a column that is supposed to read as one stack.
    % Square specifically, rather than copying the image panels' own aspect, because the image data IS
    % square here (equal x/y mm extent, see flushImageAxes) so a square box is its undistorted shape,
    % and because 4 square cells in column 1 plus a 4x4 square block for the correlation matrix makes
    % the whole layout a 5x4 grid of squares -- which is what opts.figPosition's own 5:4 aspect is
    % chosen to match, leaving no dead space to trim (Seb's own "minimize the space between all
    % graphics elements, and consequently adjust the aspect ratio").
    pbaspect(ax, [1 1 1]);
end

% ---------------------------------------------------------------------------
function H = updateFramePanels(H, k, S, fitVals, G)
    % Reassigns DATA ONLY onto the handles initFramePanels created (see that function's own note for
    % the measurement that motivated the split). Derived quantities come from the SAME signedTeff
    % helper and the SAME pass-1 stacks the init path used, so the two paths cannot drift in what they
    % display -- verified frame-for-frame against a fresh full redraw (identical captured pixels).
    set(H.patchIm, 'CData', S.patch(:,:,k));
    set(H.residIm, 'CData', S.resid(:,:,k));
    % AlphaData is NOT reassigned -- the fit mask is frame-invariant, so the masking installed by
    % initFramePanels stays correct for every frame.

    % The contour is DELETED AND RECREATED through the primitive rather than mutated: at 5 ms/frame it
    % is not worth optimizing, and going through drawGaussianPeakContour.m keeps the contour LEVEL math
    % (peak.b + peak.a*levelFracs, which moves per frame whenever a/b are themselves perFrame) in that
    % file alone instead of being partially reimplemented here as a LevelList assignment.
    %
    % CAPTURED BY DIFFING THE AXES CHILDREN, not from the primitive's return value -- deliberately, and
    % this is a trap worth spelling out: drawGaussianPeakContour.m draws ONE PLOT LINE PER CONTOUR LEVEL
    % (closed-form ellipses since 2026-08-18, no longer a contour() call) but returns only h(1), "one
    % representative handle ... the innermost/first level", per its own documented OUTPUT. Deleting only
    % that handle therefore leaves every OTHER level's line behind on every frame -- at the default two
    % levels that silently accumulated one stale ellipse per frame (350 by the end of a run), which both
    % smears the panel and slows the loop down as the axes fills. Diffing children captures all of them
    % regardless of how many levels the primitive decides to draw.
    if ~isempty(H.contour); delete(H.contour(isgraphics(H.contour))); end
    H.contour = gobjects(0);
    if strcmpi(G.method,'gaussian')
        kidsBefore = H.ax(1).Children;
        drawGaussianPeakContour(H.ax(1), fitVals, G.XIso, G.YIso, H.gold, 1.25, G.gaussLevelFrac);
        % SECONDARY peaks in RED, main in gold -- the same two-colour convention fitPatchVesselsDiag.m
        % already uses, so a peak means the same thing in both figures. Drawn through the same primitive
        % and captured by the same children-diff, so the stale-line trap above is handled once for all
        % of them. Empty for a single-component fit, in which case this loop simply does not run.
        ofv = S.otherFitVals{k};
        for q = 1:numel(ofv)
            drawGaussianPeakContour(H.ax(1), ofv(q), G.XIso, G.YIso, H.red, 0.9, G.gaussLevelFrac);
        end
        H.contour = setdiff(H.ax(1).Children, kidsBefore);
    end

    teffSigned = signedTeff(fitVals, G.paramNames, G.method, G.Xmm, G.Ymm, G.voxSz2);
    measured = S.profile(:,:,k); residual = S.resid(:,:,k);   % .profile = secondary-subtracted
    set(H.measDot,  'XData', teffSigned(H.mk), 'YData', measured(H.mk));
    set(H.residDot, 'XData', teffSigned(H.mk), 'YData', residual(H.mk));
    % The predicted CURVE moves with the frame too: its x grid is fixed but its values are this frame's
    % own 1-D-mode model, so a per-frame change in a/radius (or, with only x0/y0 free, in nothing at all)
    % is reflected honestly rather than a frame-1 curve being left behind under moving data.
    set(H.predCurve, 'YData', predictedCurve(fitVals, G));
end

% ---------------------------------------------------------------------------
function yc = predictedCurve(fitVals, G)
    % The fit's OWN model formula in 1-D mode (gaussianModel.m / parabolicModel.m 1D MODE), evaluated on
    % G.tCurve. See initFramePanels' own note on the rho*sMajor unit conversion.
    mv = toModelVals(fitVals, G.paramNames, G.method);
    s1 = mv{4}; s2 = mv{5};
    rEq = sqrt(abs(s1*s2)); if rEq==0 || ~isfinite(rEq); rEq = mean(G.voxSz2); end
    sMajor = max(abs([s1 s2]));
    yc = G.modelFun((G.tCurve/rEq)*sMajor, [], mv{:});
end

% ---------------------------------------------------------------------------
function cmap = resolveResidualColormap(spec)
    % SAME convention as plotGaussianFitPanels.m's own opts.residualColormap: [] -> this pipeline's
    % blue/red diverging map (colormap_divergingHue.m with EXPLICIT [] args, never bare -- see that
    % file's own xc note for the no-arg self-demo footgun), char -> resolveColormapFcn.m, numeric ->
    % verbatim.
    if isempty(spec)
        if exist('colormap_divergingHue','file')
            cmap = colormap_divergingHue({[250 300],[340 40]}, [0], [1], [0], [1], 256);
        else
            cmap = parula(256);
        end
        return
    end
    if isnumeric(spec); cmap = spec; return; end
    cmap = resolveColormapFcn(char(spec), '', 256);
end

% ===========================================================================
function movieFile = resolveMovieFile(v, fldName, outFld, opts)
    % The movie's own path, derived from the SAME pieces openMovie uses -- extracted so the CACHING
    % check (see the main loop) can know the target path BEFORE committing to any work. Purely a name
    % computation: creates no folder and touches no file, so it is safe to call speculatively.
    movieFile = '';
    if strcmp(opts.movie.format,'none'); return; end
    if isempty(opts.movie.fileName)
        base = sprintf('%s_%s_%s_run%s', idStrOf(v), fldName, outFld, runSpecStr(opts.run));
        base = regexprep(base, '[^\w\-]+', '_');
    else
        base = opts.movie.fileName;
    end
    ext = ['.' opts.movie.format];
    if ~endsWith(base, ext); base = [base ext]; end
    movieFile = fullfile(opts.movie.dir, base);
end

% ---------------------------------------------------------------------------
function [writer, movieFile] = openMovie(v, fldName, outFld, opts)
    writer = [];
    movieFile = resolveMovieFile(v, fldName, outFld, opts);
    if strcmp(opts.movie.format,'none'); return; end
    if ~exist(opts.movie.dir,'dir'); mkdir(opts.movie.dir); end
    switch opts.movie.format
        case 'mp4'
            ffm = resolveFfmpeg(opts.movie.ffmpeg);
            % Frames are staged as ONE FLAT RAW rgb24 FILE and encoded in a single ffmpeg pass by
            % closeMovie. Replaces per-frame lossless PNG (2026-09-03) purely on measured cost: writing
            % a frame cost 66.8 ms as PNG against 3.7 ms as a raw fwrite -- 18x -- and at ~45% of the
            % whole per-frame budget that PNG compression was the single largest item in the loop.
            % Still LOSSLESS, so nothing is given up on quality; what it costs is transient disk,
            % 4.9 MB/frame (~1.7 GB for a 350-frame run) against PNG's 0.07 MB/frame. See
            % FRAME STAGING in the file header.
            [~, stem] = fileparts(movieFile);
            stageDir = resolveStageDir(opts);
            rawFile  = fullfile(stageDir, ['.' stem '_frames.raw']);
            if exist(rawFile,'file'); delete(rawFile); end   % never append onto a stale run's frames
            fid = fopen(rawFile, 'w');
            assert(fid > 0, 'fitVesselTimeSeriesDiag:rawOpenFailed', ...
                'could not open the frame-staging file for writing: %s', rawFile);
            writer = struct('kind','mp4', 'file',movieFile, 'rawFile',rawFile, 'fid',fid, ...
                            'ffmpeg',ffm, 'frameRate',opts.movie.frameRate, 'crf',opts.movie.crf, ...
                            'size',[NaN NaN], 'nFrame',0, 'fontFile',resolveFontFile(opts.movie.fontFile), ...
                            'burnCounter',opts.movie.burnCounter);
        case 'avi'
            % Motion JPEG AVI -- VideoWriter has no H.264/MPEG-4 profile in this install, which is why
            % the 'mp4' format above goes out through ffmpeg instead of through VideoWriter.
            writer = VideoWriter(movieFile, 'Motion JPEG AVI');
            writer.FrameRate = opts.movie.frameRate;
            writer.Quality   = opts.movie.quality;
            open(writer);
        case 'gif'
            if exist(movieFile,'file'); delete(movieFile); end   % imwrite append would extend a stale file
            writer = struct('kind','gif', 'file',movieFile, 'delay',1/opts.movie.frameRate);
    end
end

% ---------------------------------------------------------------------------
function stageDir = resolveStageDir(opts)
    % Where the raw frame file goes. '' (default) = beside the movie, so a crashed run leaves its
    % staging file next to the output it was for rather than in a tempdir nobody thinks to look in.
    %
    % HARD GUARD: NEVER under /local (Seb's own explicit constraint, 2026-09-03). On this machine
    % /local is a separate 7.3 TB volume sitting at 93% full with ~558 GB free and shared with other
    % work, while /scratch (where the projects live) has its own 916 GB and /tmp is on root -- and raw
    % staging is ~1.7 GB per 350-frame movie, multiplied by however many parfor workers are running.
    % Checked on the CANONICAL path, not the literal string, so a symlink pointing into /local is
    % caught too; readlink -f resolves as far as the path exists, which is enough here.
    stageDir = opts.movie.stageDir;
    if isempty(stageDir); stageDir = opts.movie.dir; end
    if ~exist(stageDir,'dir'); mkdir(stageDir); end
    [st, canon] = system(sprintf('readlink -f %s', escapeShellArg(stageDir)));
    if st==0 && ~isempty(strtrim(canon)); canon = strtrim(canon); else; canon = stageDir; end
    assert(~startsWith(canon, '/local/') && ~strcmp(canon, '/local'), ...
        'fitVesselTimeSeriesDiag:stageDirOnLocal', ...
        ['the frame-staging directory resolves to %s, which is under /local -- refused. Raw staging ' ...
         'is ~1.7 GB per 350-frame movie (times the worker count) and /local is the shared, nearly ' ...
         'full volume on this machine. Point opts.movie.stageDir (or opts.movie.dir) somewhere else, ' ...
         'e.g. under the project on /scratch.'], canon);
end

% ---------------------------------------------------------------------------
function s = escapeShellArg(p)
    % Single-quote a path for /bin/sh, escaping any embedded single quotes. Used for every path handed
    % to system() in this file -- vessel labels reach filenames, so a stray quote or space is a real
    % possibility rather than a hypothetical.
    s = ['''' strrep(char(p), '''', '''\''''') ''''];
end

% ---------------------------------------------------------------------------
function f = resolveFontFile(spec)
    % TrueType font for the burned-in frame counter (ffmpeg's drawtext needs a real font file). ''
    % (default) probes a short list of paths present on this machine; returns '' if none is found, in
    % which case encodeMp4 skips the counter WITH A WARNING rather than failing the whole encode --
    % losing an overlay is a cosmetic degradation, losing the movie is not.
    if ~isempty(spec)
        assert(exist(spec,'file')==2, 'fitVesselTimeSeriesDiag:fontNotFound', ...
            'opts.movie.fontFile = ''%s'' does not exist.', spec);
        f = spec; return
    end
    cand = {'/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf', ...
            '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', ...
            '/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf', ...
            '/usr/share/fonts/truetype/tuffy/tuffy.ttf'};
    f = '';
    for i = 1:numel(cand)
        if exist(cand{i},'file')==2; f = cand{i}; return; end
    end
end

% ---------------------------------------------------------------------------
function ffm = resolveFfmpeg(spec)
    % Locate ffmpeg, or ERROR with something actionable -- never silently downgrade the requested
    % format to one that happens to work (see [[feedback_no_fallbacks]]): a caller who asked for mp4
    % because it is the only format they can actually VIEW is not helped by quietly getting an AVI.
    if ~isempty(spec)
        assert(exist(spec,'file')==2, 'fitVesselTimeSeriesDiag:ffmpegNotFound', ...
            'opts.movie.ffmpeg = ''%s'' does not exist.', spec);
        ffm = spec; return
    end
    [st, out] = system('command -v ffmpeg');
    assert(st==0 && ~isempty(strtrim(out)), 'fitVesselTimeSeriesDiag:ffmpegNotFound', ...
        ['opts.movie.format=''mp4'' needs ffmpeg on the PATH and it was not found. MATLAB''s own ' ...
         'VideoWriter cannot write H.264 in this install (no MPEG-4 profile), so mp4 has to go ' ...
         'through ffmpeg. Either point opts.movie.ffmpeg at the binary, or use ''gif'' (also viewable ' ...
         'in VSCode, larger, not scrubbable) or ''avi'' (small-ish, NOT viewable in VSCode).']);
    ffm = strtrim(out);
end

% ---------------------------------------------------------------------------
function refSz = settleFigureForCapture(fig)
    % Force a FULL render and absorb the first capture on a throwaway, then take the reference frame
    % size from a second one -- see the call site for the VideoWriter dimension-locking bug this exists
    % to prevent. drawnow (not drawnow limitrate) deliberately: limitrate is allowed to skip the render
    % entirely, which is exactly what leaves the first capture unsettled.
    drawnow;
    getframe(fig);        % throwaway -- absorbs the unsettled first render
    drawnow;
    fr = getframe(fig);
    refSz = size(fr.cdata, [1 2]);
end

% ---------------------------------------------------------------------------
function nResized = writeMovieFrame(writer, fig, format, isFirst, refSz)
    % isFirst is passed IN rather than tracked on the writer struct: a struct argument is passed BY
    % VALUE in MATLAB, so a 'first' flag flipped inside this function never reaches the next call --
    % every gif frame then took the LoopCount (create) branch instead of the append branch, and the
    % file ended up containing exactly ONE frame, silently. Caught by asserting imfinfo's own frame
    % count in this file's own test, not by eye (a 1-frame gif still opens and still looks plausible).
    %
    % Returns 1 if the captured frame had to be resized to refSz, 0 otherwise -- counted by the caller
    % and reported at the end of the movie rather than warning per frame.
    nResized = 0;
    if strcmp(format,'none'); return; end
    fr = getframe(fig);   % works on an INVISIBLE figure headlessly -- confirmed R2025a/Linux
    if ~isequal(size(fr.cdata,[1 2]), refSz)
        fr.cdata = imresize(fr.cdata, refSz);
        nResized = 1;
    end
    switch format
        case 'mp4'
            % Lossless RAW rgb24 append; encoded in one pass by closeMovie.
            % permute([3 2 1]) is what makes the byte order come out as rgb24 expects: MATLAB's own
            % column-major linearisation of a [ch x col x row] array walks channel fastest, then
            % column, then row -- i.e. R,G,B per pixel scanning left-to-right then top-to-bottom,
            % exactly ffmpeg's rgb24 raster order. Verified against a synthetic ramp before adopting,
            % not assumed (getting this wrong transposes or channel-swaps the whole movie silently).
            fwrite(writer.fid, permute(fr.cdata, [3 2 1]), 'uint8');
        case 'avi'
            writeVideo(writer, fr);
        case 'gif'
            [A,map] = rgb2ind(fr.cdata, 256);
            if isFirst
                imwrite(A, map, writer.file, 'gif', 'LoopCount',Inf, 'DelayTime',writer.delay);
            else
                imwrite(A, map, writer.file, 'gif', 'WriteMode','append', 'DelayTime',writer.delay);
            end
    end
end

% ---------------------------------------------------------------------------
function closeMovie(writer, format, refSz, nFrame)
    % refSz/nFrame are passed IN rather than accumulated on the writer: a struct argument is passed by
    % value, so nothing writeMovieFrame counted would survive back to here (the same by-value trap that
    % once made every gif exactly one frame long -- see writeMovieFrame's own note). The caller has both
    % numbers already.
    switch format
        case 'avi'
            close(writer);
        case 'mp4'
            fclose(writer.fid);
            encodeMp4(writer, refSz, nFrame);
    end
end

% ---------------------------------------------------------------------------
function encodeMp4(writer, refSz, nFrame)
    % One ffmpeg pass over the staged RAW rgb24 file -> H.264/yuv420p mp4, then the staging file goes.
    %
    % -f rawvideo needs the geometry stated explicitly (a flat byte stream carries none): -s WxH from
    % the SAME refSz the capture loop pinned every frame to, which is why that settle-and-pin step is
    % load-bearing here and not merely defensive -- a single differently-sized frame would silently
    % shear the whole rest of the movie rather than erroring the way VideoWriter did.
    %
    % Every flag here is load-bearing for "plays in VSCode's own built-in preview", which is the whole
    % point of this format (see MOVIE FORMATS in the file header):
    %   -c:v libx264 / format=yuv420p  H.264 8-bit 4:2:0 is what Chromium (and therefore VSCode's
    %                     Electron webview) will decode. A yuv444p or 10-bit stream is a valid mp4 that
    %                     VSCode shows as a BLACK PLAYER -- the failure looks like a broken file, not
    %                     like an unsupported profile, so it is worth pinning explicitly.
    %   pad=ceil(iw/2)*2  yuv420p needs EVEN dimensions and it is opts.figPosition that decides them.
    %                     Verified: an 801x599 figure captures at 801x599 and pads to 802x600 rather
    %                     than failing the encode outright.
    %   -color_range tv + bt709 tags  the raw stage is full-range RGB; without explicit tagging the
    %                     stream ends up marked yuvj420p (full-range, deprecated) and a player that
    %                     honours the tag differently crushes or washes the blacks -- and this figure is
    %                     mostly black background with white markers, i.e. exactly where that shows.
    %   -movflags +faststart  moves the moov atom to the front so the player can start without reading
    %                     the whole file.
    %
    % BURNED-IN FRAME COUNTER (drawtext, 'XX/YY') -- Seb's own ask, and also the reason the per-frame
    % MATLAB title is gone: re-setting that annotation's String dirtied the whole figure and cost a
    % MEASURED 32 ms/frame (a drawnow+getframe cycle went 41.5 ms with the title frozen to 73.7 ms with
    % it updated -- 44% of the cycle). ffmpeg draws it for free at encode time instead. '%{eif:n+1:d}'
    % is ffmpeg's own 0-based frame counter shifted to 1-based; the ':' inside a filter argument has to
    % reach ffmpeg escaped as '\:', which is why the text expression is built with doubled backslashes
    % here. Counts the frames ACTUALLY RENDERED, so with opts.frameStride>1 it is the position within
    % the rendered subset, not the raw series index (the static title names the run; the crosshair on
    % the correlation panel is what locates the frame within the series).
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p';
    if writer.burnCounter
        if isempty(writer.fontFile)
            warning('fitVesselTimeSeriesDiag:noFontForCounter', ...
                ['no TrueType font was found for the burned-in frame counter, so the movie is being ' ...
                 'written WITHOUT it (everything else is unaffected). Set opts.movie.fontFile to a ' ...
                 '.ttf path, or opts.movie.burnCounter=false to stop asking.']);
        else
            drawTxt = sprintf(['drawtext=fontfile=%s:text=''%%{eif\\:n+1\\:d}/%d'':' ...
                               'x=12:y=10:fontsize=22:fontcolor=white:box=1:boxcolor=black@0.55:boxborderw=6'], ...
                              writer.fontFile, nFrame);
            vf = [drawTxt ',' vf];
        end
    end
    cmd = sprintf(['%s -hide_banner -loglevel error -y ' ...
        '-f rawvideo -pix_fmt rgb24 -s %dx%d -framerate %g -i %s ' ...
        '-c:v libx264 -crf %g -preset medium -vf "%s" ' ...
        '-color_range tv -colorspace bt709 -color_primaries bt709 -color_trc bt709 ' ...
        '-movflags +faststart %s'], ...
        escapeShellArg(writer.ffmpeg), refSz(2), refSz(1), writer.frameRate, ...
        escapeShellArg(writer.rawFile), writer.crf, vf, escapeShellArg(writer.file));
    [st, out] = system(cmd);
    % The staged frames are kept ON FAILURE deliberately -- they are the expensive part (one full
    % render each) and re-running ffmpeg over them by hand is trivial, whereas re-rendering is not.
    assert(st==0, 'fitVesselTimeSeriesDiag:ffmpegFailed', ...
        ['ffmpeg failed (status %d) encoding %s\n  command: %s\n  output: %s\n  The staged RAW frames ' ...
         'were KEPT at %s so you can retry the encode without re-rendering.'], ...
        st, writer.file, cmd, strtrim(out), writer.rawFile);
    delete(writer.rawFile);
end

% ===========================================================================
function opts = defaultOpts()
    opts = struct();
    opts.run              = [];
    opts.frames           = [];
    opts.frameStride      = 1;
    opts.corrMat          = struct('source','residual', 'voxels','valid', 'clim',[0 1], ...
                                   'colormap',[], 'indicatorColor',[1 0 0], 'colorbarGap',0.045);
    opts.scale            = 10;
    opts.xlimScale        = 1;
    opts.residualColormap = [];
    opts.gaussContourPercIntegral = 95;
    opts.dotSize          = 6;
    opts.colorbarWidth    = 0.010;
    opts.leftMargin       = 0.105;
    opts.rightMargin      = 0.085;
    opts.bottomMargin     = 0.115;   % roomier than before: the correlation colorbar sits here, BELOW
                                      % the matrix's own x tick labels (see corrMat.colorbarGap)
    opts.topMargin        = 0.035;
    opts.colGap           = 0.012;
    opts.dissimGap        = 0.028;   % wider than colGap: the matrix's last x tick label and the
                                      % dissimilarity panel's first one collided at 0.012 (verified)
    opts.dissim           = struct('stat','mean');
    opts.movie            = struct('format','mp4', 'frameRate',10, 'crf',20, 'quality',75, ...
                                   'ffmpeg','', 'stageDir','', 'burnCounter',true, 'fontFile','', ...
                                   'dir',fullfile(pwd,'fitVesselTimeSeriesDiagMovies'), 'fileName','');
    opts.figPosition      = [80 80 1706 1030];   % re-solved for SIX columns, see buildQaFigure
                                                  % (both dims EVEN -- yuv420p requires it)
    opts.stayOpen         = false;
    opts.titleStr         = '';
    opts.progress         = true;
    opts.force            = false;
    opts.parallel         = true;
    opts.nPool            = [];
end

% ---------------------------------------------------------------------------
function opts = fillOptsDefaults(opts)
    assert(isstruct(opts) && isscalar(opts), 'fitVesselTimeSeriesDiag:badOpts', ...
        'opts must be a scalar struct (start from the no-arg call: opts = fitVesselTimeSeriesDiag;).');
    d = defaultOpts();
    unknown = setdiff(fieldnames(opts), fieldnames(d));
    assert(isempty(unknown), 'fitVesselTimeSeriesDiag:unknownOptsField', ...
        ['unknown opts field(s): %s -- see the no-arg call (opts = fitVesselTimeSeriesDiag;) for every ' ...
         'valid field.'], strjoin(unknown', ', '));
    % Nested groups are filled field-by-field so a caller can set ONE subfield without wiping the rest
    % (fillOptsFromDefaults.m's own convention).
    for grp = {'corrMat','movie','dissim'}
        g = grp{1};
        if ~isfield(opts,g) || isempty(opts.(g)); opts.(g) = d.(g); continue; end
        assert(isstruct(opts.(g)) && isscalar(opts.(g)), 'fitVesselTimeSeriesDiag:badOptsGroup', ...
            'opts.%s must be a scalar struct.', g);
        badSub = setdiff(fieldnames(opts.(g)), fieldnames(d.(g)));
        assert(isempty(badSub), 'fitVesselTimeSeriesDiag:unknownOptsField', ...
            'unknown opts.%s field(s): %s', g, strjoin(badSub', ', '));
        fn = fieldnames(d.(g));
        for i = 1:numel(fn)
            if ~isfield(opts.(g),fn{i}) || (isempty(opts.(g).(fn{i})) && ~strcmp(fn{i},'fileName'))
                opts.(g).(fn{i}) = d.(g).(fn{i});
            end
        end
    end
    flat = setdiff(fieldnames(d), {'corrMat','movie','dissim'});
    for i = 1:numel(flat)
        f = flat{i};
        if ~isfield(opts,f) || (isempty(opts.(f)) && ~any(strcmp(f,{'frames','titleStr','residualColormap'})))
            opts.(f) = d.(f);
        end
    end

    % --- validation (no silent coercion: a wrong value errors with its own identifier) ---
    assert(any(strcmp(opts.corrMat.source, {'residual','patch','model'})), ...
        'fitVesselTimeSeriesDiag:badCorrSource', ...
        'opts.corrMat.source must be ''residual''|''patch''|''model'' -- got ''%s''.', char(string(opts.corrMat.source)));
    assert(any(strcmp(opts.corrMat.voxels, {'valid','all'})), 'fitVesselTimeSeriesDiag:badCorrVoxels', ...
        'opts.corrMat.voxels must be ''valid''|''all'' -- got ''%s''.', char(string(opts.corrMat.voxels)));
    assert(isnumeric(opts.corrMat.clim) && numel(opts.corrMat.clim)==2 && opts.corrMat.clim(2)>opts.corrMat.clim(1), ...
        'fitVesselTimeSeriesDiag:badCorrClim', 'opts.corrMat.clim must be [lo hi] with hi>lo.');
    assert(isnumeric(opts.corrMat.colorbarGap) && isscalar(opts.corrMat.colorbarGap) ...
        && opts.corrMat.colorbarGap>=0 && opts.corrMat.colorbarGap<0.5, ...
        'fitVesselTimeSeriesDiag:badCorrColorbarGap', ...
        'opts.corrMat.colorbarGap must be a scalar in [0,0.5) (normalized figure units).');
    assert(any(strcmp(opts.dissim.stat, {'mean','median'})), 'fitVesselTimeSeriesDiag:badDissimStat', ...
        'opts.dissim.stat must be ''mean'' or ''median'' -- got ''%s''.', char(string(opts.dissim.stat)));
    assert(isnumeric(opts.corrMat.indicatorColor) && numel(opts.corrMat.indicatorColor)==3, ...
        'fitVesselTimeSeriesDiag:badIndicatorColor', 'opts.corrMat.indicatorColor must be an RGB triplet.');
    assert(any(strcmp(opts.movie.format, {'mp4','avi','gif','none'})), 'fitVesselTimeSeriesDiag:badMovieFormat', ...
        'opts.movie.format must be ''mp4''|''avi''|''gif''|''none'' -- got ''%s''.', ...
        char(string(opts.movie.format)));
    assert(isnumeric(opts.movie.frameRate) && isscalar(opts.movie.frameRate) && opts.movie.frameRate>0, ...
        'fitVesselTimeSeriesDiag:badFrameRate', 'opts.movie.frameRate must be a positive scalar.');
    assert(isnumeric(opts.movie.crf) && isscalar(opts.movie.crf) && opts.movie.crf>=0 ...
        && opts.movie.crf<=51, 'fitVesselTimeSeriesDiag:badCrf', ...
        'opts.movie.crf must be a scalar in [0 51] (x264''s own range; LOWER is better quality).');
    assert(isnumeric(opts.movie.quality) && isscalar(opts.movie.quality) && opts.movie.quality>=1 ...
        && opts.movie.quality<=100, 'fitVesselTimeSeriesDiag:badQuality', ...
        'opts.movie.quality must be a scalar in [1 100].');
    assert(ischar(opts.movie.ffmpeg) || (isstring(opts.movie.ffmpeg) && isscalar(opts.movie.ffmpeg)), ...
        'fitVesselTimeSeriesDiag:badFfmpeg', 'opts.movie.ffmpeg must be a char path ('''' = find on PATH).');
    opts.movie.ffmpeg = char(opts.movie.ffmpeg);
    for cf = {'stageDir','fontFile'}
        val = opts.movie.(cf{1});
        assert(ischar(val) || (isstring(val) && isscalar(val)), 'fitVesselTimeSeriesDiag:badMovieField', ...
            'opts.movie.%s must be a char path ('''' = default).', cf{1});
        opts.movie.(cf{1}) = char(val);
    end
    assert(islogical(opts.movie.burnCounter) && isscalar(opts.movie.burnCounter), ...
        'fitVesselTimeSeriesDiag:badBurnCounter', 'opts.movie.burnCounter must be a logical scalar.');
    assert(isnumeric(opts.frameStride) && isscalar(opts.frameStride) && opts.frameStride==round(opts.frameStride) ...
        && opts.frameStride>=1, 'fitVesselTimeSeriesDiag:badFrameStride', ...
        'opts.frameStride must be a positive integer.');
    assert(isnumeric(opts.scale) && isscalar(opts.scale) && opts.scale>0, ...
        'fitVesselTimeSeriesDiag:badScale', 'opts.scale must be a positive scalar.');
    assert(isnumeric(opts.xlimScale) && isscalar(opts.xlimScale) && opts.xlimScale>0, ...
        'fitVesselTimeSeriesDiag:badXlimScale', 'opts.xlimScale must be a positive scalar.');
    assert(isnumeric(opts.gaussContourPercIntegral) && isscalar(opts.gaussContourPercIntegral) ...
        && opts.gaussContourPercIntegral>0 && opts.gaussContourPercIntegral<100, ...
        'fitVesselTimeSeriesDiag:badGaussContourPercIntegral', ...
        ['opts.gaussContourPercIntegral must be a scalar in (0,100) -- got %s. SAME meaning and SAME ' ...
         'open interval as plotVessels.m''s own option of that name.'], mat2str(opts.gaussContourPercIntegral));
    assert(isnumeric(opts.dotSize) && isscalar(opts.dotSize) && opts.dotSize>0, ...
        'fitVesselTimeSeriesDiag:badDotSize', 'opts.dotSize must be a positive scalar.');
    assert(isnumeric(opts.colorbarWidth) && isscalar(opts.colorbarWidth) && opts.colorbarWidth>0 ...
        && opts.colorbarWidth<0.5, 'fitVesselTimeSeriesDiag:badColorbarWidth', ...
        'opts.colorbarWidth must be a scalar in (0,0.5) (normalized figure units).');
    for mf = {'leftMargin','rightMargin','bottomMargin','topMargin','colGap','dissimGap'}
        val = opts.(mf{1});
        assert(isnumeric(val) && isscalar(val) && val>=0 && val<0.5, ...
            'fitVesselTimeSeriesDiag:badMargin', ...
            'opts.%s must be a scalar in [0,0.5) (normalized figure units).', mf{1});
    end
    assert(opts.leftMargin+opts.rightMargin+opts.colGap+opts.dissimGap < 1, ...
        'fitVesselTimeSeriesDiag:badMargin', ...
        'the margins plus colGap plus dissimGap must be < 1 (they leave no width for the panels).');
    assert(opts.bottomMargin+opts.topMargin < 1, 'fitVesselTimeSeriesDiag:badMargin', ...
        'opts.bottomMargin+opts.topMargin must be < 1 (they leave no height for the panels).');
    assert(isnumeric(opts.figPosition) && numel(opts.figPosition)==4 && all(opts.figPosition(3:4)>0), ...
        'fitVesselTimeSeriesDiag:badFigPosition', 'opts.figPosition must be [x y w h] with w,h > 0.');
    assert(islogical(opts.stayOpen) && isscalar(opts.stayOpen), 'fitVesselTimeSeriesDiag:badStayOpen', ...
        'opts.stayOpen must be a logical scalar.');
    assert(islogical(opts.progress) && isscalar(opts.progress), 'fitVesselTimeSeriesDiag:badProgress', ...
        'opts.progress must be a logical scalar.');
    assert((islogical(opts.force) || (isnumeric(opts.force) && any(opts.force==[0 1]))) ...
        && isscalar(opts.force), 'fitVesselTimeSeriesDiag:badForce', ...
        'opts.force must be a logical scalar (0/1 accepted -- doIt''s own forceThis is numeric).');
    opts.force = logical(opts.force);
    assert(islogical(opts.parallel) && isscalar(opts.parallel), 'fitVesselTimeSeriesDiag:badParallel', ...
        'opts.parallel must be a logical scalar.');
    assert(isempty(opts.nPool) || (isnumeric(opts.nPool) && isscalar(opts.nPool) && opts.nPool>=1), ...
        'fitVesselTimeSeriesDiag:badNPool', 'opts.nPool must be [] or a positive integer.');
    assert(isempty(opts.run) || (isnumeric(opts.run) && isscalar(opts.run)) ...
        || ((ischar(opts.run)||isstring(opts.run)) && strcmpi(char(opts.run),'cat')), ...
        'fitVesselTimeSeriesDiag:badRun', ...
        'opts.run must be [] (one movie per run), a positive integer, or ''cat''.');
    assert(isempty(opts.frames) || (isnumeric(opts.frames) && isvector(opts.frames)), ...
        'fitVesselTimeSeriesDiag:badFrames', 'opts.frames must be [] or a numeric vector.');
end

% ---------------------------------------------------------------------------
function printOptsHelp()
    fprintf('\nFITVESSELTIMESERIESDIAG -- QA movie for a PER-FRAME vessel fit (4x5: 4 per-frame panels + a spanning frame-by-frame correlation matrix)\n');
    fprintf('  fitVesselTimeSeriesDiag(vessel, ''<fld>.<outFld>'', opts)   e.g. ''tsIm.gaussMotionEstimate''\n\n');
    fprintf('  opts.run                   : [] = ONE MOVIE PER RUN | positive integer | ''cat''         ([])\n');
    fprintf('  opts.parallel              : true | false -- parfor over runs, or vessels for ''cat''  (true)\n');
    fprintf('  opts.nPool                 : [] = profile default | positive integer worker count      ([])\n');
    fprintf('  opts.dissim.stat           : ''mean'' | ''median'' -- per-frame 1-r reduction        (''mean'')\n');
    fprintf('  opts.frames                : [] = every frame of the selected run(s) | index vector     ([])\n');
    fprintf('  opts.frameStride           : positive integer, keep every Nth frame                     (1)\n');
    fprintf('  opts.corrMat.source        : ''residual'' (fit QA) | ''patch'' (data QA) | ''model''    (''residual'')\n');
    fprintf('  opts.corrMat.voxels        : ''valid'' (the fit''s own mask) | ''all''                    (''valid'')\n');
    fprintf('  opts.corrMat.clim          : [lo hi]                                                    ([-1 1])\n');
    fprintf('  opts.corrMat.colormap      : [] = blue/red diverging | [N x 3]                           ([])\n');
    fprintf('  opts.scale                 : residual image CLim = +-patchCLim(2)/scale                   (10)\n');
    fprintf('  opts.xlimScale             : 1D XLim = image mm extent * xlimScale (1 = flush match)       (1)\n');
    fprintf('  opts.residualColormap      : [] = manuscript diverging map | ''fname(args)'' | [N x 3]    ([])\n');
    fprintf('  opts.gaussContourPercIntegral : (0,100) -- outline at this %% of integrated signal        (95)\n');
    fprintf('  opts.dotSize               : MarkerSize of the 1D dot markers                              (6)\n');
    fprintf('  opts.colorbarWidth         : thin-colorbar thickness, normalized figure units          (0.010)\n');
    fprintf('  opts.leftMargin            : strip reserved for the left colorbars + ylabels           (0.105)\n');
    fprintf('  opts.rightMargin           : strip reserved for the xcorr colorbar                     (0.085)\n');
    fprintf('  opts.bottomMargin          : strip reserved for the shared x ruler                     (0.055)\n');
    fprintf('  opts.topMargin             : strip reserved for the frame counter                       (0.035)\n');
    fprintf('  opts.colGap                : gap between column 1 and the xcorr panel                  (0.012)\n');
    fprintf('  opts.movie.format          : ''mp4'' (H.264, VIEWABLE IN VSCODE) | ''gif'' | ''avi'' | ''none''  (''mp4'')\n');
    fprintf('                               mp4 goes via ffmpeg (VideoWriter has no H.264 here); avi is NOT viewable in VSCode\n');
    fprintf('  opts.movie.frameRate       : positive scalar                                            (10)\n');
    fprintf('  opts.movie.crf             : [0 51] x264 quality, LOWER=better; mp4 only               (20)\n');
    fprintf('  opts.movie.quality         : [1 100] avi JPEG quality; avi only                         (75)\n');
    fprintf('  opts.movie.ffmpeg          : path to ffmpeg, '''' = find on PATH; mp4 only               ('''')\n');
    fprintf('  opts.movie.dir             : output folder                        (<pwd>/fitVesselTimeSeriesDiagMovies)\n');
    fprintf('  opts.movie.fileName        : '''' = auto from sId/label/fld/run                          ('''')\n');
    fprintf('  opts.figPosition           : [x y w h] pixels -- this IS the movie resolution   ([80 80 1800 1000])\n');
    fprintf('  opts.stayOpen              : true | false (the movie file is the deliverable)        (false)\n');
    fprintf('  opts.titleStr              : '''' = auto (single-vessel calls only)                      ('''')\n');
    fprintf('  opts.progress              : true | false                                             (true)\n');
    fprintf('  opts.force                 : true | false -- re-render even when the movie file exists  (false)\n\n');
end
