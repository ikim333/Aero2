%% Local coordinates from the panel

function ctrl_pan_local=coord_change_ex2(ctrl_points, X, pan_sin,pan_cos)

    ctrl_pan_local=[((ctrl_points(1)-X(1))*pan_cos-(ctrl_points(2)-X(2))*pan_sin),((ctrl_points(1)-X(1))*pan_sin+(ctrl_points(2)-X(2))*pan_cos)];
    
end
