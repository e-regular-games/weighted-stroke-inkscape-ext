function [D] = bezier_derivative(C)
  n = columns(C);
  B = factorial(n-1) ./ (factorial([1:n-1]) .* factorial(n - [1:n-1]));
  D = B .* n .* (C(:,2:n) - C(:,1:n-1));
end