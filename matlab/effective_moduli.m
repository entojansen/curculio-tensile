%Author: Andrew Jansen
%Date: January 2019
%Title: Effective moduli of Curculio cuticle

%--------------------------------------------------------------------------
%Definition of elastic constants by cuticle region
%--------------------------------------------------------------------------
%Elastic constants for exocuticle, denoted by 'x'
%E denotes Young's modulus along material axis denoted by number
Ex1 = 6.4;
Ex2 = 6.4;

%nu denotes Poisson ratio along material axes denoted by numbers
nux12 = 0.36;
nux21 = 0.36;

%G denotes shear modulus along material axes denoted by numbers
Gx12 = 2.5;

%Matrix of constants for exocuticle
exo = [Ex1 Ex2 nux12 nux21 Gx12];

%--------------------------------------------------------------------------
%For macrofiber (a.k.a. 'balken') layer in endocuticle, denoted by 'b'
Eb1 = 8.5;
Eb2 = 0.52;
nub12 = 0.49;
nub21 = 0.03;
Gb12 = 0.17;

endo = [Eb1 Eb2 nub12 nub21 Gb12];

cuticle = [exo; endo];
%--------------------------------------------------------------------------
%Cuticle lay-up: layer thicknesses, orientation angles & stacking sequence  
%--------------------------------------------------------------------------
%sequence for basal portion of head
base_z = (1 / 480) * [-24 -12:-1 1:12 24]; %48 layers, 50 microns => mm
base_theta = 30 * [0 3:8 -(-5:0); 0:5 -(-8:-3) 0]; %correct angle sequence
base_seq = [1, repelem(2, 24) 1]; %1 for exocuticle, 2 for endocuticle
base_layup = [base_z; base_theta(1:26); base_seq];

%--------------------------------------------------------------------------
%sequence for apical portion of head
apex_z = (1 / 800) * [-40 -36:4:-4 -3:-1 1:3 4:4:36 40];
apex_theta = [0 repelem(-45, 6) repelem(45, 6);
    repelem(45, 6) repelem(-45, 6) 0];
apex_seq = [1, repelem(2, 24) 1]; %same sequence as base
apex_layup = [apex_z; apex_theta(1:26); apex_seq];

%--------------------------------------------------------------------------
%permutations: exo thickness, endo angle; exo angle, endo thickness
mod_theta_layup = [base_layup(1, :); apex_layup(2, :); base_layup(3, :)];

mod_z_layup = [apex_layup(1, :); base_layup(2, :); base_layup(3, :)];

%--------------------------------------------------------------------------
%Run main function
%--------------------------------------------------------------------------
main(base_layup, mod_z_layup, mod_theta_layup, apex_layup, cuticle);

%--------------------------------------------------------------------------
%Function declarations
%--------------------------------------------------------------------------
function Q = reduced(constants)
%function to calculate 2D reduced stiffness matrix
%Extract elastic constants
E1 = constants(1);
E2 = constants(2);
nu12 = constants(3);
nu21 = constants(4);
G12 = constants(5);

%p denotes the denominator used to find elements of stiffness matrix
p = 1 - nu12 * nu21;

%E*nu is averaged to correct for rounding error in the published values
%because Eb1 * nub21 = Eb2 * nub12 is actually true for exo and endocuticle
Enu_avg = (E1 * nu21 + E2 * nu12) / 2;

Q11 = E1 / p;
Q12 = Enu_avg / p;
Q16 = 0;
Q21 = Enu_avg / p;
Q22 = E2 / p;
Q26 = 0;
Q61 = 0;
Q62 = 0;
Q66 = G12;
Q = [Q11 Q12 Q16; Q21 Q22 Q26; Q61 Q62 Q66];
end

%--------------------------------------------------------------------------
function T = transform(theta)
%function to generate transformation matrix (rotation about angle theta)
%this uses DEGREES, not radians
c = cosd(theta);
s = sind(theta);
T = [(c ^ 2) (s ^ 2) (2 * c * s);
    (s ^ 2) (c ^ 2) (-2 * c * s);
    (-c * s) (c * s) ((c ^ 2) - (s ^ 2))];
end

%--------------------------------------------------------------------------
function [A, D] = stiffness(materials, profile)
%funtion to calculate A, B, D matrices from material and layup
%materials(i,:) for the ith row, all columns (ie vector of elements in row)
A = zeros(3);
%Use this to check that B is actually zero, as needed
%B = zeros(3);
D = zeros(3);

for i = 1:size(profile, 2)
    lamina = profile(:, i);
    
    if lamina(1) < 0
        z_k1 = lamina(1);
        z_k = profile(1, (i + 1));
    else
        z_k = lamina(1);
        z_k1 = profile(1, (i - 1));
    end
    
    if (z_k == -z_k1) && (lamina(1) > 0)
        continue
    else
        theta = lamina(2);
        seq = lamina(3);
        Q = reduced(materials(seq, :));
        T = transform(theta);
        Q_k = T \ Q / T';
        A = A + Q_k * (z_k - z_k1);
        %Use this to check that B is actually zero, as needed
        %B = B + (1 / 2) * Q_k * ((z_k ^ 2) - (z_k1 ^ 2));
        D = D + (1 / 3) * Q_k * ((z_k ^ 3) - (z_k1 ^ 3));
    end
end

%Use this to check that B is actually zero, as needed
%disp(B)

end

%--------------------------------------------------------------------------
function elasticity = moduli(A, D, profile)
%funtion to calculate membrane and flexural elastic moduli along x-axis
A_star = inv(A);
D_star = inv(D);
z = profile(1, end) - profile(1, 1);

E_mx = 1 / (A_star(1) * z);
E_fx = 12 / (D_star(1) * (z ^ 3));
elasticity = [E_mx, E_fx];
end

%--------------------------------------------------------------------------
function main(base, layers, angles, apex, laminae)
%function accepts layups and produces elastic moduli
%call stiffness for each permutation
[A_base, D_base] = stiffness(laminae, base);
[A_layers, D_layers] = stiffness(laminae, layers);
[A_angles, D_angles] = stiffness(laminae, angles);
[A_apex, D_apex] = stiffness(laminae, apex);

base_moduli = moduli(A_base, D_base, base);
layer_moduli = moduli(A_layers, D_layers, layers);
angle_moduli = moduli(A_angles, D_angles, angles);
apex_moduli = moduli(A_apex, D_apex, apex);

fprintf('\nType\tMembrane\tFlexural\n')
fprintf('Base\t%6.4f GPa\t%6.4f GPa\n', base_moduli(1), base_moduli(2))
fprintf('Layers\t%6.4f GPa\t%6.4f GPa\n', layer_moduli(1), layer_moduli(2))
fprintf('Angles\t%6.4f GPa\t%6.4f GPa\n', angle_moduli(1), angle_moduli(2))
fprintf('Apex\t%6.4f GPa\t%6.4f GPa\n', apex_moduli(1), apex_moduli(2))
end














