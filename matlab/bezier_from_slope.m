
# 3 curves as represented in inkscape. (note the derivative is not continuous)
c = 2;
P = zeros(c,2, 4);
P(1,:,:) = [113.71428457 115.67558525 119.12438391 122.58252597; 42.85372511 40.73576398  37.88677577  36.48778526];
P(2,:,:) = [122.58252597 122.21719616 120.18474665 118.50788518; 36.48778526 40.30105504 44.31520358  46.67225794];

C = zeros(c, 2, 4);
S = zeros(c, 2, 2);

for i=1:c
  C(i,:,:) = bezier_from_points(squeeze(P(i,:,:)));
  S(i,:,1) = C(i,:,2) - C(i,:,1);
  S(i,:,2) = C(i,:,4) - C(i,:,3);
endfor

bezier_plot(C);

for i=1:c-1
  S(i,:,2) = (S(i,:,2) + S(i+1,:,1)) / 2;
  S(i+1,:,1) = S(i,:,2); 
endfor

S

Cp = zeros(c, 2, 4);
for i=1:c
  Si = squeeze(S(i,:,:))
  Pi = squeeze(P(i,:,:));
  Cp(i,:,:) = bezier3_match_slopes(Si(:,1), Si(:,2), Pi(:,1), Pi(:,2), Pi(:,3), Pi(:,4));
endfor

bezier_plot(Cp);