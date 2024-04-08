function [Q, R] = casteljau(C, t)
  Ci = C;
  n = columns(Ci);  
  Q = [];
  R = [];
  while n != 1
    Q = [Q, Ci(:,1)];
    R = [Ci(:,n), R];
    Ci = Ci(:,1:n-1) * (1-t) + Ci(:,2:n) * t;
    n = columns(Ci);
  end
  Q = [Q, Ci];
  R = [Ci, R];
end