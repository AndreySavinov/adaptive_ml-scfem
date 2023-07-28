function a = stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input)
    M = input(1);
    a = 1.1*ones(size(x1)); % a_0
    gap_size = 0.1;
    cookie_size = (1 - 4*gap_size)/3; 
    for m = 1:M
        % Set a_m
        if m == 1
            a_m = 1.0*(x1<=cookie_size+gap_size).*(x1>=gap_size).*(x2<=cookie_size+gap_size).*(x2>=gap_size);
        elseif m == 2
            a_m = 0.8*(x1<=2*(cookie_size+gap_size)).*(x1>=cookie_size+2*gap_size).*(x2<=cookie_size+gap_size).*(x2>=gap_size);
        elseif m == 3
            a_m = 0.4*(x1<=3*(cookie_size+gap_size)).*(x1>=2*cookie_size+3*gap_size).*(x2<=cookie_size+gap_size).*(x2>=gap_size);
        elseif m == 4
            a_m = 0.2*(x1<=cookie_size+gap_size).*(x1>=gap_size).*(x2<=2*(cookie_size+gap_size)).*(x2>=cookie_size+2*gap_size);
        elseif m == 5
            a_m = 0.1*(x1<=3*(cookie_size+gap_size)).*(x1>=2*cookie_size+3*gap_size).*(x2<=2*(cookie_size+gap_size)).*(x2>=cookie_size+2*gap_size);
        elseif m == 6
            a_m = 0.05*(x1<=cookie_size+gap_size).*(x1>=gap_size).*(x2<=3*(cookie_size+gap_size)).*(x2>=2*cookie_size+3*gap_size);
        elseif m == 7
            a_m = 0.02*(x1<=2*(cookie_size+gap_size)).*(x1>=cookie_size+2*gap_size).*(x2<=3*(cookie_size+gap_size)).*(x2>=2*cookie_size+3*gap_size);
        elseif m == 8
            a_m = 0.01*(x1<=3*(cookie_size+gap_size)).*(x1>=2*cookie_size+3*gap_size).*(x2<=3*(cookie_size+gap_size)).*(x2>=2*cookie_size+3*gap_size);
        else
            a_m = 0;
        end
        a = a + a_m * yy(m);
    end
end