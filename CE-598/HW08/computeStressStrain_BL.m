function [stress strain] = computeStressStrain ...
            (xCur, Xref, L, LatRotAngle, bondList, E, nu, planeStressFlag)
numPtcls = size(xCur, 1);
numBonds = size(bondList, 2);
% oppBond = [2 1 4 3 6 5];

LatRotAngle = LatRotAngle*pi/180.;  %convert to radians

latRotZ = [cos(LatRotAngle) sin(LatRotAngle) 0;
           -sin(LatRotAngle)  cos(LatRotAngle) 0
           0 0 1];

%for plane stress only
D = (E/(1-nu^2))* ...
      [1    nu  0;
       nu   1   0;
       0    0   1];

% D = (E/((1+nu)*(1-2*nu)))* ...
%   [(1-nu)  nu   0          nu     0            0;
%    nu    (1-nu) 0          nu     0            0;
%    0      0     (1-2*nu)/2 0      0            0;
%    nu     nu    0          (1-nu) 0            0;
%    0      0     0          0      (1-2*nu)/2   0;
%    0      0     0          0      0           (1-2*nu)/2];
   
stress = zeros(numPtcls, 6);
strain = zeros(numPtcls, 6);
stretch_3 = zeros(numBonds/2,1);
stretch_6 = zeros(numBonds,1);
N = [   1           0               0;
        (0.5)^2     (sqrt(3)/2)^2   0.5*sqrt(3)/2;
        (-0.5)^2    (sqrt(3)/2)^2   -0.5*sqrt(3)/2];
    
N_inv = inv(N);

for iPtcl = 1:numPtcls
    for jBond = 1:numBonds
        if bondList(iPtcl,jBond)~= 0
            x_ptcl = xCur(iPtcl,1);
            y_ptcl = xCur(iPtcl,2);
            x_adj = xCur(bondList(iPtcl,jBond),1);
            y_adj = xCur(bondList(iPtcl,jBond),2);
            L_star = sqrt((x_ptcl - x_adj)^2 + (y_ptcl - y_adj)^2);
            stretch_6(jBond) = (L_star - L)/L;
        end
    end
    for jBond = 1:numBonds/2
        Bond = [1 3 5];
        coBond = [2 4 6];
        if stretch_6(Bond(jBond)) && stretch_6(coBond(jBond)) ~= 0
            stretch_3(jBond) = (stretch_6(Bond(jBond)) + ...
                stretch_6(coBond(jBond)))/2;
        elseif stretch_6(Bond(jBond)) || stretch_6(coBond(jBond)) == 0
            stretch_3(jBond) = stretch_6(Bond(jBond)) + ...
                stretch_6(coBond(jBond));
        end
    end

    strain(iPtcl,1:3) = latRotZ * N_inv * stretch_3;
    strain(iPtcl,4) = (strain(iPtcl,1) + strain(iPtcl,2))/2; % out of plane
    stress(iPtcl,1:3) = latRotZ * D * transpose(strain(iPtcl,1:3));
end

return
end