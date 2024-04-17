function N = getWallNormalVector(surfaceofimpact)
% GETWALLNORMALVECTOR Get the normal vector of a surface
switch surfaceofimpact
    case 1 %"SmallWall"
        N = [1 0 0];
    case 2 %"OpSmallWall"
        N = [-1 0 0];
    case 3 %"LargeWall"
        N = [0 1 0];
    case 4 %"OpLargeWall"
        N = [0 -1 0];
    case 5 %"Floor"
        N = [0 0 1];
    case 6 %"Ceiling"
        N = [0 0 -1];
end

end