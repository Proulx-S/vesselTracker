function [offRow, offCol, nR, xSmall, ySmall] = finalCropWindow(com, cropXlim, cropYlim, cropSz)
% FINALCROPWINDOW  The final-crop small window (nR = cropSz+1 voxels) within a vessel's CURRENT
% (still large) image, centered on COM (com doesn't depend on crop size, only the window placed
% around it does).
%
% Generic geometry helper: this file has no vessel-struct field-name assumptions -- callers extract
% whatever their own vessel-struct convention calls the center point / current-window bounds and
% pass them positionally. This is the promoted-to-vfMRItools sibling of olsFaa2's own
% finalCropWindow.m (identical math, `vessel.com`/`.cropXlim`/`.cropYlim` field access left at the
% call site instead of baked in here, so the shared math has no dependency on that project's
% specific vessel-struct shape).
%
% INPUT
%   com             : [x y] center point, in the CURRENT (still-large) window's full-image pixel
%                     coordinates.
%   cropXlim,cropYlim : the CURRENT, still-large window's full-image bounds (2-element [min max]).
%   cropSz          : final patch half-window (e.g. 10 -> 11x11 voxels).
%
% OUTPUT
%   offRow, offCol : row/col offset of the small window's [1,1] corner within the CURRENT (large)
%                    image array (0-based, i.e. add 1 for a 1-based index).
%   nR             : small window side length (cropSz+1).
%   xSmall, ySmall : small window bounds in FULL-IMAGE pixel coordinates (same convention as
%                    cropXlim/cropYlim).
%
%   [offRow, offCol, nR, xSmall, ySmall] = finalCropWindow(com, cropXlim, cropYlim, cropSz)

    nR = cropSz + 1;
    xSmall = round(com(1) + [-1 1]*cropSz/2);
    ySmall = round(com(2) + [-1 1]*cropSz/2);
    offCol = xSmall(1) - cropXlim(1);
    offRow = ySmall(1) - cropYlim(1);
end
