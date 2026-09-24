function h = drawMaskOutline(ax, mask, voxSz2, color, lineStyle, lineWidth, faceAlpha)
% DRAWMASKOUTLINE  Voxel-following mask outline (same convention as every other ROI polygon in this
% pipeline) rather than a smooth contour() interpolation through pixel centers, mapped to physical
% (mm) coordinates via the SAME center-voxel-origin convention as every caller's own toX/toY
% (fitVesselProfile.m's UNITS/ORIGIN note: mm = (pixelIndex-(N+1)/2)*voxSz, N = size(mask) along that
% dimension) before being plotted onto ax -- polyOutlineOrEmpty.m returns vertices in raw crop-local
% pixel-index coordinates (pixel i at coordinate i), which is NOT yet physical mm on its own; a bare
% scale() by voxSz2 without first re-centering left every outline offset by roughly half the crop
% width/height, landing outlines off in a corner instead of aligned on the image underneath them.
% GHOST HANDLE if the mask is null/all-false (2026-08-17, Seb's own ask -- a null mask is a real,
% expected case, not just tissueMask, see maskFromThreshold.m/gaussEllipseMask.m's OUTPUT docs):
% rather than drawing nothing AND returning [], plots an EMPTY polyshape with the SAME
% color/lineStyle/lineWidth/faceAlpha -- a valid graphics handle that renders nothing on the axes
% (zero vertices) but still carries the right styling for a legend icon, so a caller building a
% legend from this handle (showVessel.m's own context overlays, drawVessel.m's own cache-hit
% confirmation and voxel-mode accept step) can label an empty roi ("this roi exists, nothing is
% currently excluded") instead of silently omitting it -- confirmed empirically that
% plot(ax,polyshape(),...) produces a real, legend-able Polygon object with a correctly colored/
% styled icon despite having nothing to actually render.
%
% Promoted out of defineLumenAndTissueRoi.m's old local plotSegmentVesselLumen subfunction --
% plotRoiFitFirstPass.m/plotRoiFitTissue.m/plotRoiFitFinal.m all need it now.
%
% INPUT
%   ax        : target axes.
%   mask      : logical mask, or [] / all-false (ghost handle, see GHOST HANDLE above).
%   voxSz2    : [row col] mm/voxel (in-plane), used to scale the outline to physical coordinates.
%   color     : RGB triplet for the outline (and the fill, when faceAlpha>0 -- always the SAME color
%               as the edge, no separate fill-color input).
%   lineStyle : default '-'.
%   lineWidth : default 1.5 (this file's original, unchanged default -- pass a smaller value for a
%               deliberately thin outline, e.g. plotRoiFitFinal.m's panel 3a).
%   faceAlpha : fill opacity, 0-1 (default 0 -- no fill, i.e. this file's original behavior,
%               unchanged for every existing caller that doesn't pass it).
%
% OUTPUT
%   h : the plotted Polygon graphics object handle -- ALWAYS valid (2026-08-17; previously [] for a
%       null mask, see GHOST HANDLE above), never [] itself. Callers that need a specific z-order
%       (e.g. a dashed reference outline drawn on top of a solid one, even after later plot calls on
%       the same axes might otherwise reorder things) can pass this to uistack(h,'top').
%
%   drawMaskOutline(ax1, mask1, voxSz2, [1 0.2 0.2]);
%   h = drawMaskOutline(ax1, manMask, voxSz2, [1 1 0.3], '-.'); uistack(h,'top');
%   h = drawMaskOutline(ax1, gaussMask, voxSz2, [1 0.85 0.1], '-', 1);        % thin
%   h = drawMaskOutline(ax1, tissueMask, voxSz2, [0.2 1 0.4], '-', 1.5, 0.2); % thin + translucent fill

    if nargin<7 || isempty(faceAlpha); faceAlpha = 0;   end
    if nargin<6 || isempty(lineWidth); lineWidth = 1.5; end
    if nargin<5 || isempty(lineStyle); lineStyle = '-'; end
    P = polyOutlineOrEmpty(mask);
    if ~isempty(P.Vertices)
        % Re-center to the SAME origin toX/toY use (pixel-index (N+1)/2 -> mm 0) before scaling to mm
        % -- see the header note above.
        [ny,nx] = size(mask);
        P = translate(P, [-(nx+1)/2, -(ny+1)/2]);
        P = scale(P, [voxSz2(2) voxSz2(1)]);
    end
    h = plot(ax, P, 'FaceColor',color, 'FaceAlpha',faceAlpha, 'EdgeColor',color, ...
        'LineWidth',lineWidth, 'LineStyle',lineStyle);
end
