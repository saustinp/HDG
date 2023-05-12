// Element sizes
h_min = 0.05;
h_max = 0.5;

//////////////////////////////////////////////////////////////////////////

Include "domain.geo";

Field[1] = Threshold;
Field[1].IField = 1;
Field[1].LcMin = h_min;
Field[1].LcMax = h_max;
Field[1].DistMin = 0.05;
Field[1].DistMax = 0.25;

Background Field = 1;

//Recombine Surface{1000};
Recombine Surface "*";

//Mesh.Algorithm = 8


