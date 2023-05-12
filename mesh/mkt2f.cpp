#ifndef __MKT2F
#define __MKT2F

// Written by: C. Nguyen & P. Fernandez

vector<Int> mkface(Int dim, Int elemtype)
{    
    vector <vector<Int> > face(nfe, vector<Int>());
    
    switch (dim) {
        case 1: // line element
            nfe = 2;                                    
            face[0].resize(1);
            face[1].resize(1);
            face[0][0] = 0;
            face[1][0] = 1;
            break;
        case 2:
            if (elemtype==0) // triangular
                nfe = dim+1;                            
                face[0].resize(2);
                face[0][0] = 1;
                face[0][1] = 2;
                
                face[1].resize(2);
                face[1][0] = 2;
                face[1][1] = 0;
                
                face[2].resize(2);
                face[2][0] = 0;
                face[2][1] = 1;                            
            else if (elemtype==1) // quadrilateral
                //face = {0,1,2,3,1,2,3,0}; // [[1,2];[2,3];[3,4];[4,1]] - 1;
                nfe = 2*dim;
                face[0].resize(2);
                face[0][0] = 0;
                face[0][1] = 1;
                
                face[1].resize(2);
                face[1][0] = 1;
                face[1][1] = 2;
                
                face[2].resize(2);
                face[2][0] = 2;
                face[2][1] = 3;                            
                
                face[3].resize(2);
                face[3][0] = 3;
                face[3][1] = 0;                            
            break;
        case 3:
            if (elemtype==0) { // tetrahedral
                // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
                nfe = dim+1;
                face[0].resize(3);
                face[0][0] = 1;
                face[0][1] = 2;
                face[0][2] = 3;
                
                face[1].resize(3);
                face[1][0] = 0;
                face[1][1] = 3;
                face[1][2] = 2;
                
                face[2].resize(3);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][1] = 3;
                
                face[3].resize(3);
                face[3][0] = 0;
                face[3][1] = 2;
                face[3][2] = 1;                
            }
            else if (elemtype==1) { // hexes
                // [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;
                nfe = 2*dim;                
                nfe = dim+1;
                face[0].resize(4);
                face[0][0] = 0;
                face[0][1] = 3;
                face[0][2] = 2;
                face[0][1] = 1;
                
                face[1].resize(4);
                face[1][0] = 4;
                face[1][1] = 5;
                face[1][2] = 6;
                face[1][3] = 7;
                
                face[2].resize(4);
                face[2][0] = 0;
                face[2][1] = 1;
                face[2][1] = 3;
                
                face[3].resize(4);
                face[3][0] = 0;
                face[3][1] = 2;
                face[3][2] = 1;                
            }
            break;
        default:
            error("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    return face;
}

void mkt2t(vector<Int> *t2t, vector<Int> *t, Int ne, Int nfe)
{
    
    
}

void mkt2f(vector<Int> *f_p, vector<Int> *t2f_p, vector<Int> *t2t_p, vector<Int> *t_p, Int dim, Int elemtype, Int ne)
{
    // Note:
    //  - The face count starts at 0 (not at 1)
    //  - All entries in t2f are non-negative, even for boundary faces
    
    Int i, j, k, l, nfe, nvf, nve;
    Int reorderFaces = 1;       // 0: No reorder. 1: Boundary first, then interior.
    
    Int *t2t = &t2t_p[0][0];
    Int *t = &t_p[0][0];
    
    Int face[];
    vector<Int> f_tmp;
    vector<Int> t2f_tmp;
    vector<Int> faceMapping;
    vector<Int> bouFaces;
    vector<Int> notBouFaces;
    vector<Int> a;
    vector<Int> b;
    vector<Int> b_sorted;
    vector<Int> c;
    
    switch (dim) {
        case 1:
            nfe = 2;
            nvf = 1;
            face = {0,1};
            break;
        case 2:
            if (elemtype==0)
                nfe = dim+1;
                nvf = dim;
                face = {1,2,0,2,0,1};      // [[2,3];[3,1];[1,2]] - 1;
            else if (elemtype==1)
                nfe = 2*dim;
                nvf = 2*(dim-1);
                face = {0,1,2,3,1,2,3,0}; // [[1,2];[2,3];[3,4];[4,1]] - 1;
            break;
        case 3:
            if (elemtype==0) {
                nfe = dim+1;
                nvf = dim;
                face = {1,0,0,0,2,3,1,2,3,2,3,1};        // [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] - 1
            }
            else if (elemtype==1) {
                nfe = 2*dim;
                nvf = 2*(dim-1);
                face = {0,4,0,2,1,3,3,5,1,3,2,0,2,6,5,7,6,4,1,7,4,6,5,7};  // [[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]] - 1;
            }
            break;
        default:
            error("Only can handle dim=1, dim=2 or dim=3\n");
    }
    
    // Compute number of faces:
    Int nf = 0;
    for (i=0; i < ne; i++)
        for j=0; j < nfe; j++)
            if (t2t[j*ne+i] > i || t2t[j*ne+i] < 0) {
                nf++;
            }
    
    // Allocate memory:
    t2f_p[0].resize(ne*nfe);
    Int *t2f = &t2f_p[0][0];
    f_p[0].resize(nf*(nvf+2));
    Int *f = &f_p[0][0];
    for (i = 0; i < f.size(); i++)
        f[i] = 0;
    a.resize(nvf*nfe);
    b.resize(nvf*nfe);
    b_sorted.resize(nvf);
    c.resize(nfe);
    
    // Compute f and t2f:
    Int ie, jf = 0;
    for (i=0; i < ne; i++) {
        for j=0; j < nfe; j++) {
            if (t2t[j*ne+i] > i || t2t[j*ne+i] < 0) {
                ie = t2t[j*ne+i];
                for (Int k = 0; k < nvf; k++)
                    f[k*nf+jf] = t[face[k*nfe+j]*ne+i];
                f[nvf*nf+jf] = i;
                f[(nvf+1)*nf+jf] = ie;
                t2f[j*ne+i] = jf;

                if (ie >= 0) {
                    for (k = 0; k < nfe; k++)
                        for (l = 0; l < nvf; l++)
                            a[k*nvf+l] = t[face[l*nfe+k]*ne+ie];
                    for (k = 0; k < nvf; k++)
                        b[k] = f[k*nf+jf];
                    sort(&b_sorted[0], &b[0], nvf);
                    for (k = 0; k < nfe; k++)
                        for (l = 0; l < nvf; l++)
                            b[k*nvf+l] = b_sorted[l];
                    for (k = 0; k < nfe; k++) {
                        c[k] = 0.0;
                        for (l = 0; l < nvf; l++)
                            c[k] += abs(a[k*nvf+l] - b[k*nvf+l]);
                        if (c[k] == 0)
                            break;
                    }
                    t2f[k*ne+ie] = jf;
                }
                jf++;
            }
        }
    }
    if (jf != nf)
        error("Error 5G7JH9M in mkt2f.cpp\n");
    
    if (reorderFaces == 1) {        // Reorder faces - First interior then boundary
        Int numBouFaces = 0, numNotBouFaces = 0;
        bouFaces.resize(nf);
        notBouFaces.resize(nf);
        for (i = 0; i < nf; i++) {
            if (f[(nvf+1)*nf+i] == 0) {
                bouFaces[numBouFaces] = i;
                numBouFaces++;
            }
            else {
                notBouFaces[numNotBouFaces] = i;
                numNotBouFaces++;
            }
        }

        f_tmp.resize(nf*(nvf+2));
        faceMapping.resize(nf);
        for (i = 0; i < numBouFaces; i++)
            for (j = 0; j < nvf+2; j++) {
                f_tmp[j*nf+i] = f[j*nf+bouFaces[i]];
                faceMapping[bouFaces[i]] = i;
            }
        for (i = numBouFaces; i < nf; i++)
            for (j = 0; j < nvf+2; j++) {
                f_tmp[j*nf+i] = f[j*nf+notBouFaces[i]];
                faceMapping[notBouFaces[i]] = i;
            }
        for (i = numBouFaces; i < nf*(nvf+2); i++)
            f[i] = f_tmp[i];

        // Modify t2f accordingly based on new face ordering:
        t2f_tmp.resize(ne*nfe);
        for (i = 0; i < ne; i++)
            for (j = 0; j < nfe; j++)
                t2f_tmp[j*ne+i] = faceMapping[t2f[j*ne+i]];
        for (i = 0; i < ne*nfe; i++)
            t2f[i] = t2f_tmp[i];
    }
}

#endif
