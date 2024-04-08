function [M] = bezier_creator_matrix(n)
  i = ones(n) .* [0:n-1]';
  j = i';
  na = i-j < 0;
  i(na) = 0;
  j(na) = 0;
  M = factorial(n-1) ./ (factorial(j) .* factorial(n - j - 1)) .* factorial(n-j - 1) ./ (factorial(i - j) .* factorial(n-i-1))  .* power(-1, i-j);
  M(na) = 0;
  M = M';
end
 
