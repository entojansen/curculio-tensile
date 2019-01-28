%Author: Andrew Jansen
%Date: January 2019
%Title: Effective moduli of Curculio cuticle

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
%p denotes the denominator used to find elements of stiffness matrix
px = 1 - nux12 * nux21;

%--------------------------------------------------------------------------
%For macrofiber (a.k.a. 'balken') layer in endocuticle, denoted by 'b'
Eb1 = 8.5;
Eb2 = 0.52;
nub12 = 0.49;
nub21 = 0.03;
Gb12 = 0.17;

function Q = lamina(S)
%Elements of 2D reduced stiffness matrix, denoted Q
pb = 1 - nub12 * nub21;
%E*nu is averaged to correct for rounding error in the published values
%because Eb1 * nub21 = Eb2 * nub12 is actually true for this material.
Enub_avg = (Eb1 * nub21 + Eb2 * nub12) / 2;

Qb11 = Eb1 / pb;
Qb12 = Enub_avg / pb;
Qb16 = 0;
Qb21 = Enub_avg / pb;
Qb22 = Eb2 / pb;
Qb26 = 0;
Qb61 = 0;
Qb62 = 0;
Qb66 = Gb12;
Q = [Qb11 Qb12 Qb16;
      Qb21 Qb22 Qb26;
      Qb61 Qb62 Qb66];
end