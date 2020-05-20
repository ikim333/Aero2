function [X, Xp, c, twist, dy] = coordinates (N, b, twist_t, c_t, c_r)

inc_theta = pi/N;

X = zeros(N+1,3);
Xp = zeros(N,3);

X(1,2) = -b/2; %first coordinate defined
twist_X(1,1) = twist_t; %twist at tip defined

for i=1:N
    
   X(i+1,2) = -(b/2)*cos(i*inc_theta); %cosine destribution
   
   dy(i,1) = X(i+1,2) - X(i,2); %section thickness
   
   Xp(i,2) = X(i+1,2) - dy(i,1)/2; %central point in each section
   
   L = b/(2-(2*c_t/c_r)); %semi_span if chord at tip is 0 (geometry)
   
   c(i,1) = c_r*(1-norm(Xp(i,2))/L); %chord of each section (geometry)
   
   twist_X(i+1,1) = twist_t*norm(cos(i*inc_theta)); %cosine distribution
   
   delta_twist(i,1) = twist_X(i+1,1) - twist_X(i,1); %inc of twist at each section
   
   twist(i,1) = twist_X(i+1,1) - delta_twist(i,1)/2; %twist at Xp
    
end

end