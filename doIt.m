% doIt.m -- minimal end-to-end example of sub-voxel vessel tracking with fitVesselPatchTimeSeries.m
%
% Pipeline: load a preprocessed full-image fMRI timeseries -> click the center pixel of one or more
% vessels -> crop one small patch per vessel -> fit a 2D Gaussian (+ flat background) to every frame of
% the patch, with each parameter's granularity chosen independently (perVessel / perRun / perFrame /
% perRunPoly) -> visualize the fits.
%
% Run it top to bottom (F5), or section by section (Ctrl+Enter). Step 3 opens an interactive figure the
% FIRST time only: the clicked centers are cached under drawVesselCenterCache/ (git-tracked), so later
% runs are non-interactive unless forceThis is raised (see step 3).
clearvars
close all
restoredefaultpath



%%%%%%%%%%%%%%%%%%%%%
%% Set up environment
%%%%%%%%%%%%%%%%%%%%%
workDir = fileparts(mfilename('fullpath'));
addpath(fullfile(workDir, 'vesselFit'));     % fitVesselPatchTimeSeries.m + its dependencies and the
                                             % two visualizers, ported from huMoMain2 (see README.md)
addpath(fullfile(workDir, 'vesselPatch'));   % loadNiftiTs.m / drawVesselCenter.m / makeVessel.m
figDir = fullfile(workDir, 'figures');       % exported figures (gitignored)
if ~exist(figDir, 'dir'); mkdir(figDir); end
%% %%%%%%%%%%%%%%%%%%
forceThis = 0;   % manual-work reuse ladder (drawVesselCenter.m's opts.force): 0 = reuse cached clicks,
                 % 1 = re-open the click figure seeded from the cache, 2 = re-open it blank



%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%
% One preprocessed (motion-corrected, dummy-scans dropped) single-slice 4D NIfTI per run. Several runs
% are simply several entries of fList; everything below is per-run aware (mode='perRun'/'perRunPoly'
% give one value / one polynomial per entry of fList). The example ships ONE run.
fList = {fullfile(workDir, 'exampleData', 'mc_run1_mag.nii.gz')};
tsIm  = loadNiftiTs(fList);
MRIread(fList)
fprintf('image %dx%d, %d run(s), %d frame(s), voxel %.2fx%.2f mm, TR %.3f s\n', ...
    size(tsIm.im{1},1), size(tsIm.im{1},2), numel(tsIm.im), size(tsIm.im{1},4), tsIm.vSize(1), tsIm.vSize(2), tsIm.dt(1));
%% %%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select vessel centers (manual)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Click the center pixel of one or more vessels on the time-averaged image; Enter to accept. The
% accepted list is cached per opts.id, so this only prompts the first time (or when forceThis >= 1).
opts = drawVesselCenter;
opts.id       = 'example';
opts.cacheDir = fullfile(workDir, 'drawVesselCenterCache');
opts.force    = forceThis;
opts.zoomXlim = [];   % initial view only -- zoom/pan freely from the figure toolbar
opts.zoomYlim = [];
center = drawVesselCenter(tsIm, opts);
disp(center)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the vessel variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One (cropSz+1)-pixel square patch per clicked center, in the vessel struct shape
% fitVesselPatchTimeSeries.m reads (vessel{v}.tsIm.im / .vSize, plus sId/label for figure titles).
opts = makeVessel;
opts.cropSz = 10;            % 11 x 11 pixels around the clicked center
opts.sId    = 'example';
vessel = makeVessel(tsIm, center, opts);
clear tsIm                   % the full image is no longer needed
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%


imagesc(mean(vessel{1}.tsIm.im{:},4))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit -- time-averaged anatomy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One Gaussian + background on the TIME-AVERAGED patch (opts.timeAvg=true): the anatomical starting
% point every fit below is seeded from. The default seed ('heuristic') needs no prior fit -- a start
% is derived from the patch itself, centered on the patch center, i.e. the clicked pixel. Bounds are
% opened up (Inf) so the anatomy can settle wherever the data says.
opts = fitVesselPatchTimeSeries('gaussian + background');
opts.outFld  = 'gaussAnat';
opts.timeAvg = true;
opts.model.gaussian.seedBounds.x0y0        = inf;
opts.model.gaussian.seedBounds.theta       = inf;
opts.model.gaussian.seedBounds.a           = inf;
opts.model.gaussian.seedBounds.radius      = inf;
opts.model.gaussian.seedBounds.aspectRatio = inf;
opts.model.background.seedBounds.b         = inf;
vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);
vessel{1}.tsIm.gaussAnat.fits   % the fitPatchVessels.m-shaped view of a time-averaged fit

%%% Visualize the anatomy fit
% (a) Every vessel in one grid: time-averaged patch + the fitted main-peak contour (showVessel.m).
opts = showVessel;
showVessel(vessel, {'tsIm.gaussAnat.fits'}, opts);
exportgraphics(gcf, fullfile(figDir, 'gaussAnat_contour.png'), 'Resolution', 150);
% (b) One six-panel diagnostic per vessel (plotGaussianFitPanels.m): measured patch, fitted model at
% isochromat resolution, residual, two-sided radial profiles of the data and of the model, and a
% measured-vs-predicted scatter. buildGaussianFitDiag.m rebuilds the model image from the stored fit
% (run 1, frame 1 -- a time-averaged fit has exactly one frame per run).
opts = plotGaussianFitPanels;
opts.xlimScale = 6;           % radial-profile x-range, in fitted equivalent-area radii
% Numeric residual colormap (blue-white-red): the function's own default calls util's
% colormap_divergingHue, which this standalone tool does not carry.
opts.residualColormap = interp1([0 0.5 1], [0.1 0.4 1; 1 1 1; 1 0.3 0.1], linspace(0, 1, 256).');
for v = 1:numel(vessel)
    [diagS, titleStr] = buildGaussianFitDiag(vessel{v}, 'tsIm', 'gaussAnat', 1, 1);
    plotGaussianFitPanels(diagS, titleStr, opts);
    exportgraphics(gcf, fullfile(figDir, sprintf('%s_%s_gaussAnat_panels.png', vessel{v}.sId, vessel{v}.label)), 'Resolution', 150);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit -- perVessel / perRun / perFrame / perRunPoly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each fit is seeded from the anatomy fit above via the 2-arg call (opts derived from the existing
% result; its .mode/.seedBounds copied), then only the granularity is changed. Every parameter has
% its own mode -- mix them freely (e.g. position perFrame, shape perVessel).
%   'perVessel'          one value for the whole vessel, across all runs (joint fit of every frame)
%   'perRun'             one value per run
%   'perFrame'           one value per frame
%   'perRunPoly([0 1 2])' a Legendre polynomial in within-run time, its own coefficients per run


% %%% perVessel -- the anatomy refit jointly on every raw frame (no time averaging)
% opts = fitVesselPatchTimeSeries(vessel, 'tsIm.gaussAnat');
% opts.outFld  = 'gaussPerVessel';
% opts.timeAvg = false;
% vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);

% %%% perRun -- every parameter free per run (with one run this equals perVessel; add runs to fList)
% opts = fitVesselPatchTimeSeries(vessel, 'tsIm.gaussAnat');
% opts.outFld  = 'gaussPerRun';
% opts.timeAvg = false;
% for p = gaussParam; opts.model.gaussian.mode.(p{1}) = 'perRun'; end
% opts.model.background.mode.b = 'perRun';
% vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);

vessel{1}.tsIm
%%% perFrame -- vessel POSITION tracked frame by frame, shape and background shared across frames
opts = fitVesselPatchTimeSeries(vessel, 'tsIm.gaussAnat');
opts.outFld  = 'gaussPerFrame';
opts.model.gaussian.mode.x0     = 'perFrame';
opts.model.gaussian.mode.y0     = 'perFrame';
opts.model.gaussian.mode.a      = 'perFrame';
opts.model.gaussian.mode.radius = 'perFrame';
vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);

figure
params = {'x0','y0','a','radius'};
for p = 1:numel(params)
    subplot(2,2,p); hold on
    plot(0:2:(size(vessel{1}.tsIm.gaussPerFrame.(params{p}),2)*2-1),vessel{1}.tsIm.gaussPerFrame.(params{p}));
    yline(vessel{1}.tsIm.gaussAnat.(params{p}));
    ylabel(params{p})
    grid on
end





%%% perRunPoly -- position as a smooth quadratic drift within each run
opts = fitVesselPatchTimeSeries(vessel, 'tsIm.gaussAnat');
opts.outFld  = 'gaussPerRunPoly';
opts.timeAvg = false;
opts.model.gaussian.mode.x0     = 'perRunPoly([0 1 2])';
opts.model.gaussian.mode.y0     = 'perRunPoly([0 1 2])';
opts.model.gaussian.mode.a      = 'perFrame';
opts.model.gaussian.mode.radius = 'perFrame';
vessel = fitVesselPatchTimeSeries(vessel, 'tsIm', opts);

figure
params = {'x0','y0','a','radius'};
for p = 1:numel(params)
    subplot(2,2,p); hold on
    plot(0:2:(size(vessel{1}.tsIm.gaussPerRunPoly.(params{p}),2)*2-1),vessel{1}.tsIm.gaussPerRunPoly.(params{p}));
    yline(vessel{1}.tsIm.gaussAnat.(params{p}));
    ylabel(params{p})
    grid on
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%
%% Visualize the fit
%%%%%%%%%%%%%%%%%%%%
%%% Fitted position over time, all granularities on one axis per vessel (mm, patch-centered)
for v = 1:numel(vessel)
    fig = figure('Name', sprintf('%s %s -- fitted position', vessel{v}.sId, vessel{v}.label), 'NumberTitle', 'off');
    for iParam = 1:2
        param = {'x0','y0'}; param = param{iParam};
        ax = subplot(2, 1, iParam, 'Parent', fig); hold(ax, 'on'); grid(ax, 'on');
        t0 = 0;
        for r = 1:numel(vessel{v}.tsIm.im)
            T  = size(vessel{v}.tsIm.im{r}, 4);
            t  = t0 + (0:T-1) * vessel{v}.dt(r);
            yF = perRunVal(vessel{v}.tsIm.gaussPerFrame.(param),   r);   % [1 x T]
            yP = perRunVal(vessel{v}.tsIm.gaussPerRunPoly.(param), r);   % [1 x T] (polynomial expanded per frame)
            yR = perRunVal(vessel{v}.tsIm.gaussPerRun.(param),     r);   % scalar
            yV = perRunVal(vessel{v}.tsIm.gaussPerVessel.(param),  r);   % scalar
            hF = plot(ax, t, yF, '.', 'Color', [0.6 0.6 0.6]);
            hP = plot(ax, t, yP, '-', 'Color', [0.9 0.5 0.1], 'LineWidth', 2);
            hR = plot(ax, t([1 end]), [yR yR], '-', 'Color', [0.1 0.6 0.9], 'LineWidth', 1.5);
            hV = plot(ax, t([1 end]), [yV yV], '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
            t0 = t(end) + vessel{v}.dt(r);
        end
        ylabel(ax, [param ' (mm)']);
        if iParam == 1
            title(ax, sprintf('%s %s -- fitted position (mm, relative to patch center)', vessel{v}.sId, vessel{v}.label), 'Interpreter', 'none');
            legend(ax, [hF hP hR hV], {'perFrame', 'perRunPoly([0 1 2])', 'perRun', 'perVessel'}, 'Location', 'northeast');
        else
            xlabel(ax, 'time (s)');
        end
    end
    exportgraphics(fig, fullfile(figDir, sprintf('%s_%s_position.png', vessel{v}.sId, vessel{v}.label)), 'Resolution', 150);
end

%%% Per-frame QA movie of the perFrame fit (fitVesselTimeSeriesDiag.m): patch + fitted outline,
%%% residual, 1D profiles, and a frame-by-frame residual correlation matrix, animated over frames.
%%% 'gif' needs no external binary; 'mp4' (the default, scrubbable) needs ffmpeg on the PATH.
opts = fitVesselTimeSeriesDiag;
opts.movie.format = 'gif';
opts.movie.dir    = fullfile(workDir, 'qaMovies');
opts.frameStride  = 5;       % every 5th frame keeps the example quick; 1 = every frame
opts.parallel     = false;   % one vessel at a time (no Parallel Computing Toolbox needed)
opts.force        = 1;
vessel = fitVesselTimeSeriesDiag(vessel, 'tsIm.gaussPerFrame', opts);
{vessel{1}.tsIm.gaussPerFrame.qaMovie.fName}'
%% %%%%%%%%%%%%%%%%%%%
