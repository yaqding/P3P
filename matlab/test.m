% mex solver_p3p.cpp %(Require Eigen3) comment this line after completed
clc;
clear;
X = rand(3,3);
R = orth(rand(3,3));
T = rand(3,1);
Y = R*X+T;
x2d = normc(Y);
data = [X(:,1); x2d(:,1); X(:,2); x2d(:,2); X(:,3); x2d(:,3)];
sols = solver_p3p(data);
for i = 1:size(sols,2)/4
    Re = sols(:,4*i-3:4*i-1);
    Te = sols(:,4*i);
    err_r(i) = 2*asind(norm(R-Re)/(2*sqrt(2)));
    err_t(i) = norm(T-Te);
end