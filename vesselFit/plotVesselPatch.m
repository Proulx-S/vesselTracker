function [toX, toY, cLim] = plotVesselPatch(ax, im, voxSz2)
% PLOTVESSELPATCH  Draws a vessel patch (this pipeline's own term -- a square roi image of a central
% vessel cross-section) into a GIVEN axes, in physical mm, center-voxel-relative -- the one shared
% low-level primitive for the pixel-drawing mechanics every vessel-patch-panel caller in this pipeline
% needs (imagesc + axis 'image' + grayscale colormap + a degenerate-image-safe CLim + no ticks + hold
% on). Promoted out of showVessel.m's own showOneVessel subfunction -- fitVesselDiag.m needed the exact
% same handful of lines and had grown its own independent copy instead (see this project's own
% multi-vessel/multi-metric terminology note on this drift).
%
% Deliberately UNOPINIONATED about everything else a caller might want layered on top -- roi outlines
% (drawMaskOutline.m), model contours, centroid markers, titles -- draw those into the SAME axes
% afterward using the returned toX/toY, so any further mm-space plotting stays in the SAME coordinate
% convention this panel used. Same "caller owns the axes/figure lifecycle" convention as
% drawMaskOutline.m -- this function never creates a figure or axes of its own.
%
% XLim/YLim LOCKED TIGHT (2026-09-07, Seb's own report -- an roi outline overlay pushed xlim/ylim
% beyond the patch): imagesc+axis('image') sets XLim/YLim to the tight patch extent but leaves
% XLimMode/YLimMode at 'auto' (SAME root cause fitPatchVesselsDiag.m's own axPatch/axIsoOverlay panels
% already worked around locally, see that file's own "Lock the DEFAULT view" note) -- any overlay drawn
% afterward whose own vertices sit outside that extent (a mask outline's own pixel-EDGE vertices can
% land fractionally past the image's pixel-CENTER extent; likewise a fit-peak contour or footprint box)
% then auto-EXPANDS the axes to fit it. Explicitly re-applying XLim/YLim and switching both to 'manual'
% here, once, centralizes the fix in this shared primitive instead of leaving every caller (this
% pipeline's roi/fit-peak/footprint overlay drawers included) to rediscover and re-fix it individually.
% A caller that genuinely wants a wider default view (e.g. plotVessels.m's own vesselPatchFrame
% zoomWindow override) still can -- an explicit xlim(ax,...)/ylim(ax,...) call afterward works exactly
% the same in 'manual' mode as it did in 'auto'.
%
% INPUT
%   ax     : target axes (an existing axes handle).
%   im     : [ny x nx] image to display -- ALREADY reduced to one 2D frame (time-/run-averaged,
%            censored, whichever THIS caller's own convention is) -- this function does no further
%            reduction of its own.
%   voxSz2 : [row col] mm/voxel, in-plane -- same convention as fitVessel.m/showVessel.m/drawVessel.m.
%
% OUTPUT
%   toX, toY : @(x) / @(y) pixel-index -> physical-mm conversion functions (center-voxel-relative), for
%              the caller's own further mm-space overlays to stay in the SAME coordinate space this
%              panel was drawn in.
%   cLim     : [min max] (or [c c-1 c+1] broadened, if im is uniform) -- the CLim actually applied to ax,
%              for a caller that wants a SECOND, related axes (e.g. a model-image panel) to share
%              exactly the same color scale rather than recomputing this same degenerate-image guard.
%
%   [toX, toY] = plotVesselPatch(ax, imAvg, [1 1]);
%   drawMaskOutline(ax, roiMask, [1 1], [1 1 1], '-', 1.5);          % further overlay, same coordinate space
%   [toX, toY, cLim] = plotVesselPatch(ax1, im, voxSz2);
%   set(ax2, 'CLim',cLim);                                          % a second panel sharing the same scale

    im = double(im);
    [ny,nx] = size(im);
    ctrCol = (nx+1)/2; ctrRow = (ny+1)/2;   % center-voxel-relative mm -- see fitVessel.m's own POSITION CONVENTION
    toX = @(x) (x-ctrCol)*voxSz2(2);
    toY = @(y) (y-ctrRow)*voxSz2(1);
    cLim = [min(im(:)) max(im(:))]; if diff(cLim)==0; cLim = cLim + [-1 1]; end

    cla(ax);
    imagesc(ax, toX(1:nx), toY(1:ny), im); axis(ax,'image'); colormap(ax,gray);
    set(ax, 'CLim',cLim, 'XTick',[], 'YTick',[]); hold(ax,'on');
    % Lock the tight extent 'axis image' just computed -- see XLim/YLim LOCKED TIGHT above.
    set(ax, 'XLim',xlim(ax), 'YLim',ylim(ax), 'XLimMode','manual', 'YLimMode','manual');
end
