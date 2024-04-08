function curves = bezier_cubic_divide(C)
  n = columns(C);
  M = bezier_creator_matrix(n);
  curves = zeros(n-1, 2, 4);
  for i=1:n-1
    t = 1.0 / (n-i);
    [Q, R] = casteljau(C, t);
    C = R;
    
    s0 = Q(:,2) - Q(:,1);
    s2 = Q(:,n) - Q(:,n-1);
    s = [s0, s2];
    q = [bezier_eval(Q, 1/3), bezier_eval(Q, 2/3)];

    Q3 = [Q(:,1), [0 0; 0 0], Q(:,n)];
  
    K = inv(s) * (27 * q - [Q3(:,1), Q3(:,4)] * [20 7; 7 20]) * inv([12 6; -6 -12]);
    a = K(1,1);
    b = K(2,2);

    Q3(:,2) = Q3(:,1) + a * s0;
    Q3(:,3) = Q3(:,4) - b * s2;

    curves(i,:,:) = Q3;
  endfor
end