function KL_gradcoeff = stochcol_diffusion_grad_x1_spatial_expansion(x1, x2, yy, input)
M = input(1);
alpha_bar = input(2);
sigma_tilde = input(3);
KL_gradcoeff = zeros(size(x1)); % a_0_dx
for m = 1:M
    km = floor(-0.5e0+sqrt(0.25e0+2*m));
    beta_x1 = m -km*(km + 1)/2;
    beta_x2 = km -beta_x1;
    if m == 1
        KL_gradcoeff = KL_gradcoeff - 2*pi*beta_x1*0.498/(m^(sigma_tilde))*sin(2*pi*beta_x1*x1).*cos(2*pi*beta_x2*x2)*yy(m);
    else
        KL_gradcoeff = KL_gradcoeff - 2*pi*beta_x1*alpha_bar/(m^(sigma_tilde))*sin(2*pi*beta_x1*x1).*cos(2*pi*beta_x2*x2)*yy(m);
    end
end
end