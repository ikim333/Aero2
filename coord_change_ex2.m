function coord_punt_control=coord_change_ex2(control, x, s_panell,c_panell)

    coord_punt_control=[((control(1)-x(1))*c_panell-(control(2)-x(2))*s_panell),((control(1)-x(1))*s_panell+(control(2)-x(2))*c_panell)];
    
end
