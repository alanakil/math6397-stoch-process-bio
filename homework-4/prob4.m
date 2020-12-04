%% Hwk 4 Prob 4b

syms eps
assume(eps,'real')

% eps=0.00001;

%% TFT-TFT
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps) eps.^2;...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps);...
    eps*(1-eps), (1-eps).^2, eps.^2, eps*(1-eps);...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_TFT_TFT = [ 1-sum(y), y(1), y(2), y(3) ] 

%% TFT-GRIM
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps) eps.^2;...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2;...
    eps*(1-eps), (1-eps).^2, eps.^2, eps*(1-eps);...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_TFT_GRIM = [ 1-sum(y), y(1), y(2), y(3) ] 


%% TFT-ALLC
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps) eps.^2;...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps);...
    (1-eps).^2, eps*(1-eps), eps*(1-eps), eps.^2;...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps)];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_TFT_ALLC = [ 1-sum(y), y(1), y(2), y(3) ] 


%% GRIM-GRIM
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps) eps.^2;...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2;...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2;...
    eps.^2, eps*(1-eps), eps*(1-eps), (1-eps).^2];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_GRIM_GRIM = [ 1-sum(y), y(1), y(2), y(3) ] 


%% GRIM-ALLC
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps) eps.^2;...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps);...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps);...
    eps*(1-eps), eps.^2, (1-eps).^2, eps*(1-eps)];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_GRIM_ALLC = [ 1-sum(y), y(1), y(2), y(3) ] 

%% ALLC-ALLC
M = [ (1-eps).^2, eps*(1-eps), eps*(1-eps), eps.^2;...
    (1-eps).^2, eps*(1-eps), eps*(1-eps), eps.^2;...
    (1-eps).^2, eps*(1-eps), eps*(1-eps), eps.^2;...
    (1-eps).^2, eps*(1-eps), eps*(1-eps), eps.^2];

% We can simplify and solve the system using A and b.
A = [ 1+M(1,2)-M(2,2), M(1,2)-M(3,2), M(1,2)-M(4,2);...
    M(1,3)-M(2,3), 1+M(1,3)-M(3,3), M(1,3)-M(4,3);...
    M(1,4)-M(2,4), M(1,4)-M(3,4), 1+M(2,4)-M(4,4)];

b = [ M(1,2); M(1,3); M(1,4)];

y = A\b;

% Stationary distribution.
pi_ALLC_ALLC = [ 1-sum(y), y(1), y(2), y(3) ] 

