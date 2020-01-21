# pca analysis
pkg load statistics
data = csvread ('unified.csv');

d1 = data(2:size(data)(1),2:size(data)(2));
d2 = d1 > 0;

d3 = d1 (sum (d2,2) > 6, :);

[a, b] = princomp (transpose(d3), 'econ');

scatter (b(:,1), b(:,2))
