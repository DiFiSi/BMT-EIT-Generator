function [x,y,z] = getCoords(fmdl)
    % 4 for 2D, 3 for 3D
    dims = size(fmdl.nodes,2);
    interp_no = 6 - dims;
    try 
        interp_no = fmdl.elem_select.interp_no;
    end

    pts = interp_mesh(fmdl, interp_no);
    x = squeeze(pts(:,1,:));
    y = squeeze(pts(:,2,:));
    if dims == 2
        z = 0 * x;
    else
        z = squeeze(pts(:,3,:));
    end
end