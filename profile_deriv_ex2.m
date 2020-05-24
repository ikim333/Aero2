function dz=profile_deriv_ex2(x,p,m)
if x<p
    dz=(2*m/p^2)*(p-x);
else
    dz=(2*m/(1-p)^2)*(p-x);
end
end


