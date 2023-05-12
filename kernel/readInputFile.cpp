#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

struct solstruct {
    vector<double> UDG;
    vector<double> UH;
};

struct meshstruct {
    vector<double> dgnodes;
    vector<int> elcon;
    vector<int> bf;
    vector<int> f;
    vector<int> t2f;
    vector<int> f2f;
    vector<int> perm;
    vector<int> permgeom;
};
 
struct masterstruct {
    vector<double> plocvl;
    vector<double> gpvl;
    vector<double> gwvl;
    vector<double> plocfc;
    vector<double> gpfc;
    vector<double> gwfc;
    vector<double> shapvt;
    vector<double> shapvg;
    vector<double> shapvgdotshapvl;
    vector<double> shapft;
    vector<double> shapfg;
    vector<double> shapfgdotshapfc;
    vector<double> shapmv;
    vector<double> shapmf;    
};

struct appstruct {
    vector<int> bcm;
    vector<double> bcs;
    vector<int> bcd;
    vector<double> bcv;
    vector<double> dt;
    vector<double> param;
    vector<int> flag;
    vector<double> factor;
};

struct intstruct {
    int ndomains;    
};

void readdarray(ifstream &in, vector<double> &a, int N)
{
    a.resize(N);
  
    in.read( reinterpret_cast<char*>( &a[0] ), sizeof(double)*N );        
}

void readiarray(ifstream &in, vector<int> &a, int N)
{
    a.resize(N);
    in.read( reinterpret_cast<char*>( &a[0] ), sizeof(int) * N );       
}

void readiarrayfromdouble(ifstream &in, vector<int> &a, int N)
{
    a.resize(N);
  
    double read;
    for (unsigned i = 0; i < N; i++) {
        in.read( reinterpret_cast<char*>( &read ), sizeof read );
        a[i] = (int) round(read);            
    }
}

void readInput(char* filename, vector<int> &ndims, meshstruct &mesh, masterstruct &master, 
        appstruct &app, intstruct &intf, solstruct &sol) 
{    
    // open this file to read
    ifstream in(filename, ios::in | ios::binary);
    
    //If file is not valid
    if (!in) {
        cout <<"Unable to open file" << filename << endl;
    }
    
    // If this is a valid file
    if (in) {
        int N = 100;                
        readiarrayfromdouble(in, ndims, N);       
        
        int nd = ndims[0];
        int ncd = ndims[1];
        int nfe = ndims[2];
        int nve = ndims[3];
        int nvf = ndims[4];
        int ne = ndims[5];
        int nf = ndims[6];
        int nv = ndims[7];
        int ndh = ndims[8]; 
        int npe = ndims[9];
        int npf = ndims[10];
        int nme = ndims[11];
        int nmf = ndims[12];
        int nge = ndims[13]; 
        int ngf = ndims[14];
        int porder = ndims[15]; 
        int morder = ndims[16]; 
        int torder = ndims[17]; 
        int nstage = ndims[18]; 
        int nc = ndims[19];
        int ncu = ndims[20];
        int ncq = ndims[21];
        int ncp = ndims[22];
        int nch = ndims[23];
        int ns  = ndims[24];
        int nb  = ndims[25];
        int ndt = ndims[26];
        int nparam = ndims[27];
        int nflag = ndims[28];
        int nfactor = ndims[29];
        
        N = npe*nc*ne;          
        readdarray(in, sol.UDG, N);
        N = nch*ndh;
        readdarray(in, sol.UH, N);
        
        N = nme*ncd*ne;
        readdarray(in, mesh.dgnodes, N);
        N = npf*nfe*ne;
        readiarrayfromdouble(in, mesh.elcon, N);                        
        N = nfe*ne;
        readiarrayfromdouble(in, mesh.bf, N);                                
        N = nf*(nvf+2);
        readiarrayfromdouble(in, mesh.f, N);                                
        N = ne*nfe;
        readiarrayfromdouble(in, mesh.t2f, N);                                
        N = nf*(2*nfe-1);
        readiarrayfromdouble(in, mesh.f2f, N);                                
                
        N = nfe*npf;
        readiarrayfromdouble(in, mesh.perm, N);                                
        N = nmf*nfe;
        readiarrayfromdouble(in, mesh.permgeom, N);                                      
        N = npe*nd;
        readdarray(in, master.plocvl, N);        
        N = nge*nd;
        readdarray(in, master.gpvl, N);                
        N = nge;
        readdarray(in, master.gwvl, N);                
        N = npf*(nd-1);
        readdarray(in, master.plocfc, N);        
        N = ngf*(nd-1);
        readdarray(in, master.gpfc, N);                        
        N = ngf;
        readdarray(in, master.gwfc, N);                                       
        N = nge*npe*(nd+1);
        readdarray(in, master.shapvt, N);        
        N = npe*nge*(nd+1);
        readdarray(in, master.shapvg, N);                
        N = npe*npe*nge*(nd+1);
        readdarray(in, master.shapvgdotshapvl, N);                
        N = ngf*npf*nd;
        readdarray(in, master.shapft, N);        
        N = npf*ngf*nd;
        readdarray(in, master.shapfg, N);                
        N = npf*npf*ngf*nd;
        readdarray(in, master.shapfgdotshapfc, N);                
        N = nge*nme*(nd+1);
        readdarray(in, master.shapmv, N);        
        N = ngf*nmf*nd;
        readdarray(in, master.shapmf, N);                        
        
        N = nb;
        readiarrayfromdouble(in, app.bcm, N);                                     
        N = nb*nch;
        readdarray(in, app.bcs, N);                        
        N = nb;
        readiarrayfromdouble(in, app.bcd, N);                                     
        N = nb*nch;
        readdarray(in, app.bcv, N);                        
        N = ndt;
        readdarray(in, app.dt, N);                        
        N = nparam;
        readdarray(in, app.param, N);                        
        N = nflag;
        readiarrayfromdouble(in, app.flag, N);                                     
        N = nfactor;
        readdarray(in, app.factor, N);     
        
        intf.ndomains = ns;        
    }
    
    in.close();
}

void writedarray(ofstream &out, vector<double> &a, int N)
{    
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );       
}

void writeiarray(ofstream &out, vector<int> &a, int N)
{    
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(int) * N );       
}

void writeiarraytodouble(ofstream &out, vector<int> &a, int N)
{    
    double b;
    for (unsigned i = 0; i < N; i++) {
        b = (double) a[i];  
        out.write( reinterpret_cast<char*>( &b ), sizeof(double) );                  
    }
}

void writeOutput(char* filename, vector<int> &ndims, meshstruct &mesh, masterstruct &master, 
        appstruct &app, intstruct &intf, solstruct &sol) 
{    
    // open this file to read
    ofstream out(filename, ios::out | ios::binary);
    
    //If file is not valid
    if (!out) {
        cout <<"Unable to open file" << filename << endl;
    }
    
    // If this is a valid file
    if (out) {
        int N = 100;                
        writeiarraytodouble(out, ndims, N);                
        
        int nd = ndims[0];
        int ncd = ndims[1];
        int nfe = ndims[2];
        int nve = ndims[3];
        int nvf = ndims[4];
        int ne = ndims[5];
        int nf = ndims[6];
        int nv = ndims[7];
        int ndh = ndims[8]; 
        int npe = ndims[9];
        int npf = ndims[10];
        int nme = ndims[11];
        int nmf = ndims[12];
        int nge = ndims[13]; 
        int ngf = ndims[14];
        int porder = ndims[15]; 
        int morder = ndims[16]; 
        int torder = ndims[17]; 
        int nstage = ndims[18]; 
        int nc = ndims[19];
        int ncu = ndims[20];
        int ncq = ndims[21];
        int ncp = ndims[22];
        int nch = ndims[23];
        int ns  = ndims[24];
        int nb  = ndims[25];
        int ndt = ndims[26];
        int nparam = ndims[27];
        int nflag = ndims[28];
        int nfactor = ndims[29];
        
        N = npe*nc*ne;
        writedarray(out, sol.UDG, N);
        N = nch*ndh;
        writedarray(out, sol.UH, N);
        
        N = nme*ncd*ne;
        writedarray(out, mesh.dgnodes, N);               
        N = npf*nfe*ne;
        writeiarraytodouble(out, mesh.elcon, N);                        
        N = nfe*ne;
        writeiarraytodouble(out, mesh.bf, N);                                
        N = nf*(nvf+2);
        writeiarraytodouble(out, mesh.f, N);                                
        N = nfe*ne;
        writeiarraytodouble(out, mesh.t2f, N);                                
        N = nf*(2*nfe-1);
        writeiarraytodouble(out, mesh.f2f, N);                                
        N = nfe*npf;
        writeiarraytodouble(out, mesh.perm, N);                                        
        N = nmf*nfe;
        writeiarraytodouble(out, mesh.permgeom, N);                                     
                
        N = npe*nd;
        writedarray(out, master.plocvl, N); 
        N = nge*nd;
        writedarray(out, master.gpvl, N);      
        N = nge;
        writedarray(out, master.gwvl, N);    
        N = npf*(nd-1);
        writedarray(out, master.plocfc, N);  
        N = ngf*(nd-1);
        writedarray(out, master.gpfc, N);   
        N = ngf;
        writedarray(out, master.gwfc, N);     
        N = nge*npe*(nd+1);
        writedarray(out, master.shapvt, N);        
        N = npe*nge*(nd+1);
        writedarray(out, master.shapvg, N);     
        N = npe*npe*nge*(nd+1);
        writedarray(out, master.shapvgdotshapvl, N);   
        N = ngf*npf*nd;
        writedarray(out, master.shapft, N);    
        N = npf*ngf*nd;
        writedarray(out, master.shapfg, N);       
        N = npf*npf*ngf*nd;
        writedarray(out, master.shapfgdotshapfc, N); 
        N = nge*nme*(nd+1);
        writedarray(out, master.shapmv, N);   
        N = ngf*nmf*nd;
        writedarray(out, master.shapmf, N);                        
        
        N = nb;
        writeiarraytodouble(out, app.bcm, N);      
        N = nb*nch;
        writedarray(out, app.bcs, N);                                
        N = nb;
        writeiarraytodouble(out, app.bcd, N);    
        N = nb*nch;
        writedarray(out, app.bcv, N);  
        N = ndt;
        writedarray(out, app.dt, N); 
        N = nparam;
        writedarray(out, app.param, N);                                
        N = nflag;
        writeiarraytodouble(out, app.flag, N); 
        N = nfactor;
        writedarray(out, app.factor, N);             
    }
    
    out.close();
}

int main(int argc, char** argv) {
   
    if( argc == 3 ) {  
      cout << argv[1] << endl;
      cout << argv[2] << endl;
    }
    else {
      cout << "Usage: ./cppfile InputFile OutputFile\n";
      return 1;
    }
    
    vector<int> ndims; 
    meshstruct mesh;
    masterstruct master; 
    appstruct app;
    intstruct intf;
    solstruct sol;    
            
    readInput(argv[1], ndims, mesh, master, app, intf, sol);
    writeOutput(argv[2], ndims, mesh, master, app, intf, sol);
    
    return 0;
     
}
