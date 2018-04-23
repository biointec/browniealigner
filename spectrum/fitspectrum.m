load spectrum.txt;
hold off;
plot(spectrum(:,1),spectrum(:,2));
axis([0, 80, 0, 4000000]);
x = 1:80;

clear components;
components(1,:) = [1.073 1.078 1e9];
components(2,:) = [10.89 12.89 8.16e8];
components(3,:) = [26.95 55.7 4.83e7];
components(4,:) = [51.91 79.68 7.353e7];

components = mixtureModel(spectrum(1:80,:), components, 150)

hold on;

y = zeros(1,length(x));
for i = 1:size(components,1)
  y += components(i,3)*negbinomial(x, components(i,1), components(i,2));
  plot(x, components(i,3)*negbinomial(x, components(i,1), components(i,2)), 'linestyle', ':');
endfor

plot(x, y, 'linestyle', '-.', 'color', 'r');
