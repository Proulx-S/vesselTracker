function [h, Z, levelFracs] = drawGaussianPeakContour(ax, peak, XIso, YIso, color, lineWidth, levelFracs)
% DRAWGAUSSIANPEAKCONTOUR  Draws one rotated-2D-Gaussian peak's own iso-contour (gaussianModel.m) onto
% ax, at an isochromat-resolution grid (buildIsochromatGrid.m) the CALLER builds once and passes in --
% same "caller owns the axes" convention as drawMaskOutline.m, and same "grid built once, reused across
% every peak on one patch" convention fitPatchVesselsDiag.m/fitVessel.m already established for
% buildIsochromatGrid.m itself (this file does not call buildIsochromatGrid.m internally, so multiple
% peaks on the SAME patch -- e.g. a main + several secondary peaks -- share one grid construction
% rather than each rebuilding it).
%
% Promoted out of fitPatchVesselsDiag.m's own drawThreeMetricPanels, which had two near-identical inline
% gaussianModel+contour call sites (main peak gold, each secondary peak red/flagged-gray) -- both now
% call this instead, keeping their own surrounding bookkeeping (peak-index text labels, .excluded
% flagging, ZfitOther/otherContribVox accumulation) local to that file, since none of that generalizes.
% showVessel.m's own fit-type overlay support (see its file header PATH SPEC note) is the second caller.
%
% INPUT
%   ax      : target axes.
%   peak    : struct with .a/.x0/.y0/.radius/.aspectRatio/.theta/.b (fitVessel.m's/fitPatchVessels.m's
%             own 'gaussian' output convention -- radius=sqrt(sx*sy), aspectRatio=sx/sy, see fitVessel.md
%             -- this project's fits have no literal .sx/.sy field at all; converted to gaussianModel.m's
%             own native sx,sy right here, as plain local variables, never stored anywhere).
%   XIso, YIso : isochromat-resolution mm grid (buildIsochromatGrid.m), SAME size, caller-built.
%   color   : RGB triplet for the contour line.
%   lineWidth : contour LineWidth. Default 0.5 (this pipeline's own established thinnest-practical
%             value -- see fitPatchVesselsDiag.m's own contour history).
%   levelFracs : contour level(s), as fraction(s) of the peak's own amplitude (peak.b + peak.a*levelFracs
%             is what's actually contoured) -- fractions, not absolute levels, so a caller can report
%             what the drawn contour actually MEANS (e.g. in a legend/title, see contourLevelFracStr.m)
%             without re-deriving it from peak.a/.b itself. Default [exp(-0.5) 0.9] (2026-08-17,
%             narrowed 2026-08-17 -- was [0.1 exp(-0.5) 0.9], the inner 10% level dropped as visual
%             clutter once the radius level below made it redundant). exp(-0.5)=~61% is the peak's own
%             1-SIGMA ellipse: the isocontour where (x'/sx)^2+(y'/sy)^2=1 in the peak's own rotated frame,
%             which by construction has semi-axes EXACTLY sx,sy -- so its own equivalent-area radius
%             (sqrt of the two semi-axes' product) is EXACTLY peak.radius (=sqrt(sx*sy)), not an
%             approximation. Included by DEFAULT (not an opt-in extra) so any default display of the
%             main peak's own gold contour -- this file's own two callers, fitPatchVesselsDiag.m's panel
%             1 and showVessel.m's own FIT-type overlay, both omit this argument and rely on the default
%             -- visually shows the contour that actually corresponds to the fitted radius value, not
%             just an arbitrary outer bracket (Seb's own report, 2026-08-16: a screenshot showing
%             radius's own bounds-table value with no way to see where that radius actually sits on the
%             patch itself).
%
% OUTPUT
%   h : one representative line handle for the drawn contour (the innermost/first level -- callers only
%       use this for identity/legend purposes, e.g. showVessel.m's own legHandles, same "one handle
%       stands for the whole multi-level contour" behavior a single contour() object already had).
%   Z : the evaluated model surface (gaussianModel.m on XIso/YIso) -- returned so a caller that also
%       needs it for something else (e.g. fitPatchVesselsDiag.m's own composite panel) doesn't have to
%       evaluate gaussianModel.m a second time.
%   levelFracs : echoes the INPUT (resolved to its own default if the caller left it empty/omitted) --
%       so a caller can label what was actually drawn (e.g. contourLevelFracStr(levelFracs)) even when
%       it didn't pass an explicit value itself.
%
%   [XIso,YIso] = buildIsochromatGrid(7, nx, ny, (nx+1)/2, (ny+1)/2, voxSz2);
%   [h,~,levelFracs] = drawGaussianPeakContour(ax, peak, XIso, YIso, [1 0.85 0.1], 0.5);

    if nargin<7 || isempty(levelFracs); levelFracs = [exp(-0.5) 0.9]; end
    if nargin<6 || isempty(lineWidth); lineWidth = 0.5; end
    sqrtAR = sqrt(peak.aspectRatio); sx = peak.radius*sqrtAR; sy = peak.radius/sqrtAR;
    Z = gaussianModel(XIso, YIso, peak.a, peak.x0, peak.y0, sx, sy, peak.theta, peak.b);

    % Drawn as CLOSED-FORM ellipses, NOT via contour() sampled on XIso/YIso (2026-08-18, replacing the
    % former [~,h]=contour(ax,XIso,YIso,Z,levels,...) call here). A 2D anisotropic rotated Gaussian's
    % own iso-density level set IS exactly an ellipse -- gaussianModel's formula reduces at
    % v=b+a*levelFrac to (xr/(sx*rho))^2+(yr/(sy*rho))^2=1 with rho=sqrt(-2*log(levelFrac)), (xr,yr) the
    % (X,Y) frame rotated by theta around (x0,y0) -- so sampling that level on a FINITE grid is not just
    % unnecessary, it actively fails when the peak's own fitted center lies outside the grid's own
    % sampled extent: contour() still returns a valid handle (per this function's own former OUTPUT
    % note) but its ContourMatrix comes back EMPTY -- confirmed against a real fitted vessel (a
    % bound-clipped edge peak, x0 converged at its own lower bound, ~2x the patch's own half-width
    % outside XIso/YIso's domain) -- meaning NOTHING was ever actually drawn, regardless of axes
    % XLim/YLim/zoom, since there was no sampled data to contour in the first place. The analytic
    % ellipse below has no such limitation: it draws correctly no matter where (x0,y0) sits relative to
    % XIso/YIso, which remains passed in/used only for this function's own Z return value (still needed
    % by callers that accumulate/composite Z further, e.g. fitPatchVesselsDiag.m's own composite panel).
    t = linspace(0, 2*pi, 100);
    h = gobjects(1, numel(levelFracs));
    for i = 1:numel(levelFracs)
        rho = sqrt(-2*log(levelFracs(i)));
        xr = rho*sx*cos(t); yr = rho*sy*sin(t);
        xw = peak.x0 + xr*cos(peak.theta) - yr*sin(peak.theta);
        yw = peak.y0 + xr*sin(peak.theta) + yr*cos(peak.theta);
        h(i) = plot(ax, xw, yw, 'Color',color, 'LineWidth',lineWidth);
    end
    h = h(1);
end
