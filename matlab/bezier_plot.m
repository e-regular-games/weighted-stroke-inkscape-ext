function bezier_plot(curves)
  c = 1;
  if columns(size(curves)) == 3
    c = size(curves)(1);
  else
    tmp = zeros(1, rows(curves), columns(curves));
    tmp(1,:,:) = curves;
    curves = tmp;
  endif
  
  hold on;
  for i=1:c
    C = squeeze(curves(i,:,:));
    n = columns(C);
    M = bezier_creator_matrix(n);
    
    dist = 0;
    for j=2:n
      dist += norm(C(:,j) - C(:,j-1));
    endfor
    
    resolution = 1 / (dist * 100);
    T = power([0:resolution:1], [0:n-1]');
    
    P = C * M * T;
    
    plot(P(1,:), P(2,:), "-")    
  endfor
  hold off;
end
  