%% Point extrapolation (Lagrangian representation)
function y0 = approximate(x,y,x0)
n = length(x);
y0 = zeros(1,size(y,2));
for k = 1:n
   w = 1;
   for j = [1:k-1 k+1:n]
      w = (x0-x(j))./(x(k)-x(j)).*w;
   end
   y0 = y0 + w*y(k,:);
end
