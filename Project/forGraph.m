clc;
clear all;
steps = 100; % quantity of steps for grid
step = 1 / steps;  % step of grid
checkTime = zeros(1, 2);
tic;

Const = (3 * cosh(sqrt((pi.^2 + 1) / 3)) / (pi.^2 + 1) - 3 / (pi.^2 + 1)) / sinh(sqrt((pi.^2 + 1) / 3));
U = zeros(steps + 1, steps + 1);
for ix = 1 : steps + 1
    for iy = 1 : steps + 1
%         U(ix, iy) = sin(pi * (ix - 1) * step) * (-3 * cosh((iy - 1) * step * sqrt((pi.^2 + 1) / 3)) / (pi.^2 + 1) + 3 / (pi.^2 + 1) + Const * sinh((iy - 1) * step * sqrt((pi.^2 + 1) / 3)));
        U(ix, iy) = sin(pi * (ix - 1) * step) * exp((ix - 1) * step) * sin(pi * (iy - 1) * step);
    end
end
[X, Y] = meshgrid(0 : step : 1, 0 : step : 1);
figure
surf(X, Y, U);
title('Exact solution');
ylabel('Y'), xlabel('X'), zlabel('U(X, Y)');
colorbar;
savefig('ES');

fSin = zeros((steps - 1), (steps - 1));
for iy = 1 : steps - 1
    for ix = 1 : steps - 1
        fSin(ix, iy) = step.^2 * 2 * pi * exp(ix * step) * sin(pi * iy * step) * (2 * pi * sin(pi * ix * step) - cos(pi * ix * step));
%         fSin(ix, iy) = 3 * step.^2 * sin(pi * ix * step);
    end
end
checkTime(1,1) = toc;
% multiply by W^(-1)
fSin = fSin';
uSin = idst(fSin);
uSin = uSin';
uSin = idst(uSin) / step.^2;

lambda = zeros((steps - 1), (steps - 1));
for m = 1 : steps - 1
    for k = 1 : steps - 1
        lambda(k, m) = step.^2 / (4 * (sin(pi * k * step * 0.5)).^2 + 12 * ...
        (sin(pi * m * step * 0.5)).^2 + step.^2);
    end
end
uSin = lambda .* uSin;

% multiply by W
uSin = uSin';
uSin = dst(uSin);
uSin = uSin';
uSin = dst(uSin);
checkTime(1,2) = toc;
disp("Time for Double fast sine conversion: " ), disp(checkTime(1,2) - checkTime(1,1));
USin = zeros(steps + 1, steps + 1);
for iy = 1 : steps - 1
    for ix = 1 : steps - 1
        USin(ix + 1, iy + 1) = uSin(ix, iy);
    end
end

% figure
% surf(X, Y, USin);
% title('Double fast sine conversion');
% ylabel('Y'), xlabel('X'), zlabel('USin(X, Y)');
% colorbar;
% savefig('DFSC');

dif2 = U - USin;
figure
surf(X, Y, dif2);
title('Difference between U & USin');
ylabel('Y'), xlabel('X'), zlabel('dif2(X, Y)');
colorbar;
savefig('Dif2');

%% constructing right part of equation
fJ = zeros(steps + 1, steps + 1);
for iy = 1 : steps - 1
    for ix = 1 : steps - 1
        fJ(ix + 1, iy + 1) =  step.^2 * 2 * pi * exp(ix * step) * sin(pi * iy * step) * (2 * pi * sin(pi * ix * step) - cos(pi * ix * step));
    end
end

%% Jacobi method
UJ = zeros(steps + 1, steps +1);
iterations = 0;
UJPrev = UJ;

while (max(max(U - UJ)) > step^2)
   for iy = 2 : steps
       for ix = 2 : steps
           UJ(ix, iy) = (fJ(ix, iy) + UJPrev(ix - 1, iy) + UJPrev(ix + 1, iy) + 3 * UJPrev(ix, iy - 1) + 3 * UJPrev(ix, iy + 1)) / (8 + step.^2);
       end
   end
   UJPrev = UJ;
   iterations = iterations + 1;
end

figure
difJ = U - UJ;
surf(X, Y, difJ);
title('Difference between U & UJ');
ylabel('Y'), xlabel('X'), zlabel('difJ(X, Y)');
colorbar
savefig('Dif3');

