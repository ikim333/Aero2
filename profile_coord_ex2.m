%% Profile coordinates

function z=profile_coord_ex2(x,p,m)
    if x<p
        z=(m/p^2)*(2*p*x-x^2);
    else
        z=(m/(1-p)^2)*(1-2*p+2*x*p-x^2);
    end
end
