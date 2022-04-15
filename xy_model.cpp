#include<iostream>
#include<math.h>
#include<cstdlib>
#include<stdio.h>
#include<ctime>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<string>
#include<ctype.h>
#include<unistd.h>
#include<sstream>

using namespace std;

typedef unsigned long  int ULONG;
typedef unsigned short int USHORT;
typedef signed long int    LONG;
typedef signed short int   SHORT;
typedef double REAL;

#include "random.hpp"



//file to output data into
ofstream DATA("DATA_4.dat", ios::out);
ofstream CI("C(i)_4.dat", ios::out);
ofstream SF("Rho_sf.dat", ios::out);
fstream G_cur;
fstream C_cur;
fstream EI;
fstream MI;
fstream RHO;
fstream bond_cur;
fstream Stat;

int d; //dimensions
int L; //lattice size
double T; //temperature
long long int transient;
long long int sweeps;
long long int Z1_max;
long long int Z_max; 
int step_z; //number of Z between saving statistics
int step_p; //number of measurements between printing
int step_w; //number of measurements between writing to files  2^20
const double delta=0.0001; // the acceptance of error difference to stop the binning loop

void read_params(){
  ifstream inFile("params_XY");
  string s;
  inFile >> d; getline(inFile, s);
  inFile >> L; getline(inFile, s);
  inFile >> transient; getline(inFile,s);
  inFile >> sweeps; getline(inFile,s);
  inFile >> Z1_max; getline(inFile,s);
  inFile >> step_z; getline(inFile,s);
  inFile >> step_p; getline(inFile,s);
  inFile >> step_w; getline(inFile,s);
}

//associations
int ass(int &id, int &site){
int i3=pow(L, abs(id));
int i4=pow(L, (abs(id)-1));
int i5=(id/abs(id));
int ass = site + i5*i4; // inside associations
int i6 = (site/i3);
int i7 = ((site-i4)/i3);
int i8 = ((site+i4)/i3);
if(id>0){
if((site%i3)==0){
ass = ass - i3;	
}
// just for 3D and upper dimensions
if(((site+i4)%i3)!=0&&i6!=i8){
ass = ass - i3;	
}
}
if(id<0){
if(((site-i4)%i3)==0){
ass = ass + i3;	
}
// just for 3D and upper dimensions
if((site%i3)!=0&&i6!=i7){
ass = ass + i3;	
}
}
if(ass > pow(L, d)){
ass = ass - i3;
}
if(ass < 1){
ass = ass + i3;
}
return ass;
}

// Associates an integer index with a vector
int NUM_VEC(int vec[]){
int num = 1;
for(int j6=1;j6<=d;j6++){
int j7=pow(L, (j6-1))*vec[j6];
num=num+j7;
}
return num;
}

// integrate
double integrate(double &p){ 
double Q=1/(2*p+2);
return Q;
}


//main program
int main(int argc, char *argv[])
{

const int N = 60;
double temp[N];
for (int h = 0; h<N; h++)

{

for (int c=0; c<100; c++)
{

    temp[h]= 0.4+h*0.028;
    printf("%f\n",temp[h]);



Random ran;
ran("random");
read_params();
T=temp[h];
double K=(1/(2*T)); //inverse temperature (Js^2)/(2T)
int n=pow(L, d); //number of spin points on lattice
transient=transient*2*d*32768;
sweeps=sweeps*2*d*32768;
Z1_max=Z1_max*2*d*32768;
step_w=step_w*2*d*32768;
step_z=step_z*2*d*16384;


long long int it=0;
Z_max=Z1_max+transient+sweeps;
long long int length=((Z_max-transient-sweeps)/step_z); //length of the E, Msq array
//initialization
int ira=1;
int masha=ira;
long double num_bonds=0;
int** bond;
bond = new int* [(2*d+1)];
for(int x0=0;x0<(2*d+1);x0++){
bond[x0] = new int[(n+1)]; 
}
int vec[(d+1)]; // vector to measure distance
for(int i=1;i<=(2*d);i++){
for(int j=1;j<=n;j++){
bond[i][j]=0;	
}	
}
for(int k=1;k<=d;k++){
vec[k]=0;
}

//initializing statistics
long double Z=0; 
long double av_bonds=0; // av_bonds = average number of bonds times Z
long double sq_bonds=0; // sum of (num_bonds^2-num_bonds)
long double* G = new long double[(n+1)]; // spin correlator
long double* Msq = new long double[(length+1)]; // Magnetization squared per spin
long double* E = new long double[(length+1)]; // Energy per site
long double Msq_av=0;
long double E_av=0; // average of energy per site
long double Msq_sigma=0;
long double E_sigma=0; // error of energy per site
long double* C = new long double[(length+1)]; // specific heat
long double C_av=0;
long double C_sigma=0;

for(int i0=1;i0<=n;i0++){
G[i0]=0;	
}
for(int i1=1;i1<=length;i1++){
Msq[i1]=0;
E[i1]=0;
}
for(int h0=1;h0<=length;h0++){
C[h0]=0;
}

// Winding number & rho_sf
long double* num_winding = new long double[(d+1)];
long double av_winding=0;
long double* rho_sf = new long double[(length+1)];
long double rho_sf_av=0;
long double rho_sf_sigma=0;
for(int w1=1;w1<=d;w1++){
num_winding[w1]=0;
}
for(int w4=1;w4<=length;w4++){
rho_sf[w4]=0;
}

// read config & stat
bond_cur.open("bond_cur");

if(bond_cur){
bond_cur.close();
Stat.open("Stat", ios::binary | ios::in);
Stat.read(reinterpret_cast<char*>(&Z), sizeof(Z));
Stat.read(reinterpret_cast<char*>(&ira), sizeof(ira));
Stat.read(reinterpret_cast<char*>(&masha), sizeof(masha));
Stat.read(reinterpret_cast<char*>(&num_bonds), sizeof(num_bonds)); 
Stat.read(reinterpret_cast<char*>(&av_bonds), sizeof(av_bonds));
Stat.read(reinterpret_cast<char*>(&sq_bonds), sizeof(sq_bonds));
Stat.read(reinterpret_cast<char*>(&it), sizeof(it)); 
for(int b1=1;b1<=d;b1++){
Stat.read(reinterpret_cast<char*>(&vec[b1]), sizeof(vec[b1])); 
}
for(int b7=1;b7<=d;b7++){
Stat.read(reinterpret_cast<char*>(&num_winding[b7]), sizeof(num_winding[b7])); 
}
Stat.read(reinterpret_cast<char*>(&av_winding), sizeof(av_winding));
Stat.close();
bond_cur.open("bond_cur", ios::binary | ios::in);
int d1;
int d2;
int d3;
for(int b2=1;b2<=(2*d*n);b2++){
d1=(b2/n);
d2=d1+1;
d3=(b2%n);
if(d3!=0){
bond_cur.read(reinterpret_cast<char*>(&bond[d2][d3]), sizeof(bond[d2][d3])); 
}
if(d3==0){
bond_cur.read(reinterpret_cast<char*>(&bond[d1][n]), sizeof(bond[d1][n]));
}
}
bond_cur.close();
G_cur.open("G_cur", ios::binary | ios::in); 
for(int b3=1;b3<=n;b3++){
G_cur.read(reinterpret_cast<char*>(&G[b3]), sizeof(G[b3])); 
}
G_cur.close();
C_cur.open("C_cur", ios::binary | ios::in); 
for(int b10=1;b10<=length;b10++){
C_cur.read(reinterpret_cast<char*>(&C[b10]), sizeof(C[b10])); 
}
C_cur.close();
EI.open("EI", ios::binary | ios::in); 
for(int b4=1;b4<=length;b4++){
EI.read(reinterpret_cast<char*>(&E[b4]), sizeof(E[b4])); 
}
EI.close();
MI.open("MI", ios::binary | ios::in);
for(int b5=1;b5<=length;b5++){
MI.read(reinterpret_cast<char*>(&Msq[b5]), sizeof(Msq[b5])); 
}
MI.close();
RHO.open("RHO", ios::binary | ios::in);
for(int b8=1;b8<=length;b8++){
RHO.read(reinterpret_cast<char*>(&rho_sf[b8]), sizeof(rho_sf[b8])); 
}
RHO.close();
}
else{
bond_cur.close();
}

// The main MC loop
while(Z<Z_max){

//printf("%Lf = %lld\n", Z, Z_max);



it=it+1;
//SHIFT(ira, masha, num_bonds, vec, bond);
if(ira==masha){
ira=(int)ceil(ran()*n);
while(ira==0){
ira=(int)ceil(ran()*n);
}
masha=ira;
}
// bond[id][site]: an arrow from this site to the next one following "positive" direction id; e.g. bond[1][1]: 1 to 2
// bond[id+d][site]: an arrow from this site to the next one following "negative" direction id; e.g. 2D: bond[1+2][2]: 2 to 1
// prepare for update 
int id1;
int site;
int nb; // outgoing current
int mb; // incoming current
int charge1; // |nb-mb| before update
int charge2; // |nb-mb| after update
int v;

// choose direction
double id0=((2*ran()-1)*d);
while(id0==0){
id0=((2*ran()-1)*d);
}

if(id0<0){ // negative direction
id1=(int)floor(id0);
site=ass(id1, ira); // the site which ira will shift if the update is accepted
mb=bond[abs(id1)][site];
nb=bond[(abs(id1)+d)][ira];
v=vec[abs(id1)]-1;
if(v==-1){
v=L-1;	
}
}
if(id0>0){ // positive direction
id1=(int)ceil(id0);
site=ass(id1, ira); // the site which ira will shift if the update is accepted
nb=bond[abs(id1)][ira];
mb=bond[(abs(id1)+d)][site];
v=vec[abs(id1)]+1;
if(v==L){
v=0;	
}
}


double rand=ran();
while(rand==0.5){
rand=ran();
}

// Draw
if(rand<0.5&&ran()<(K/(nb+1))){
nb=nb+1;
if(id1>0){
bond[abs(id1)][ira]=nb;
num_winding[abs(id1)]=num_winding[abs(id1)]+1;
}
if(id1<0){
bond[(abs(id1)+d)][ira]=nb;
num_winding[abs(id1)]=num_winding[abs(id1)]-1;
}
vec[abs(id1)]=v;
num_bonds=num_bonds+1;
ira=site;	
}
// Erase
if(rand>0.5&&ran()<(mb/K)){
mb=mb-1;
if(id1>0){
bond[(abs(id1)+d)][site]=mb;
num_winding[abs(id1)]=num_winding[abs(id1)]+1;
}
if(id1<0){
bond[abs(id1)][site]=mb;
num_winding[abs(id1)]=num_winding[abs(id1)]-1;
}
vec[abs(id1)]=v;
num_bonds=num_bonds-1;
ira=site;
}


//MEASURE(ira, masha, vec, G, Z, num_bonds, av_bonds, E, Msq);
// update the spin correlator
if(Z>=transient){
int m=NUM_VEC(vec);
G[m]=G[m]+1;
}
// update Z
if(ira==masha){
Z=Z+1;

if(Z>transient){
av_bonds=av_bonds+num_bonds;
sq_bonds=sq_bonds+(num_bonds*num_bonds-num_bonds);
for(int w5=1;w5<=d;w5++){
av_winding=av_winding+(num_winding[w5]/L)*(num_winding[w5]/L);
}
/*long double Winding[(d+1)];
for(int w9=1;w9<=d;w9++){
Winding[w9]=0;
}
int site0;
for(int w2=1;w2<=d;w2++){
for(int w3=1;w3<=n;w3++){
site0=ass(w2, w3);
Winding[w2]=Winding[w2]+bond[w2][w3]-bond[(w2+d)][site0];
}
Winding[w2]=Winding[w2]/L;
av_winding=av_winding+Winding[w2]*Winding[w2];
}*/
}

// store statistics after number of Z
long long int Z1=Z;

if(Z1>(transient+sweeps)&&((Z1-transient-sweeps)%step_z)==0){
long long int j1=((Z1-transient-sweeps)/step_z);
long long int Z2=Z1-transient;
// Energy
E[j1]=-T*(av_bonds/Z2);
// Energy per site
E[j1]=E[j1]/n;
// Specific heat
C[j1]=(sq_bonds/Z2)-(av_bonds/Z2)*(av_bonds/Z2);
// Specific heat per site
C[j1]=C[j1]/n;
//printf("%Lf\n", C[j1]);
// rho(s)=(2mT/d)*L^(2-d)*<M^2>
rho_sf[j1]=(T/d)*pow(L, (2-d))*(av_winding/Z2);
// Sum all spin correlators
for(int j2=1;j2<=n;j2++){
Msq[j1]=Msq[j1]+G[j2];	
}
// Magnetization squared per spin
Msq[j1]=Msq[j1]/(Z2*n);
if(j1>0&&(j1%step_p)==0){
printf("it = %Ld\n", it);
printf("Z = %Lf\n", Z);
printf("E[%Ld] = %Lf\n", j1, E[j1]);
printf("Msq[%Ld] = %Lf\n", j1, Msq[j1]);
printf("C[%Ld] = %Lf\n\n", j1, C[j1]);
}
}

// write config to files
if((Z1%step_w)==0){
Stat.open("Stat", ios::binary | ios::out | ios::trunc); 
Stat.write(reinterpret_cast<char*>(&Z), sizeof(Z));
Stat.write(reinterpret_cast<char*>(&ira), sizeof(ira));
Stat.write(reinterpret_cast<char*>(&masha), sizeof(masha));
Stat.write(reinterpret_cast<char*>(&num_bonds), sizeof(num_bonds)); 
Stat.write(reinterpret_cast<char*>(&av_bonds), sizeof(av_bonds));
Stat.write(reinterpret_cast<char*>(&sq_bonds), sizeof(sq_bonds));
Stat.write(reinterpret_cast<char*>(&it), sizeof(it));
for(int a1=1;a1<=d;a1++){
Stat.write(reinterpret_cast<char*>(&vec[a1]), sizeof(vec[a1])); 
}
for(int a9=1;a9<=d;a9++){
Stat.write(reinterpret_cast<char*>(&num_winding[a9]), sizeof(num_winding[a9])); 
}
Stat.write(reinterpret_cast<char*>(&av_winding), sizeof(av_winding));
Stat.close();
bond_cur.open("bond_cur", ios::binary | ios::out | ios::trunc); 
for(int a2=1;a2<=(2*d);a2++){
for(int a3=1;a3<=n;a3++){
bond_cur.write(reinterpret_cast<char*>(&bond[a2][a3]), sizeof(bond[a2][a3])); 
}
}
bond_cur.close();
G_cur.open("G_cur", ios::binary | ios::out | ios::trunc); 
for(int a4=1;a4<=n;a4++){
G_cur.write(reinterpret_cast<char*>(&G[a4]), sizeof(G[a4])); 
}
G_cur.close();
EI.open("EI", ios::binary | ios::out | ios::trunc); 
MI.open("MI", ios::binary | ios::out | ios::trunc);
RHO.open("RHO", ios::binary | ios::out | ios::trunc);
C_cur.open("C_cur", ios::binary | ios::out | ios::trunc);
for(int a5=1;a5<=length;a5++){
EI.write(reinterpret_cast<char*>(&E[a5]), sizeof(E[a5])); 
MI.write(reinterpret_cast<char*>(&Msq[a5]), sizeof(Msq[a5]));
RHO.write(reinterpret_cast<char*>(&rho_sf[a5]), sizeof(rho_sf[a5]));
C_cur.write(reinterpret_cast<char*>(&C[a5]), sizeof(C[a5]));
}
EI.close();
MI.close();
RHO.close();
C_cur.close();
ran.flush("random.iseed");
}

}
}

//AVERAGES(E, Msq, E_av, Msq_av, E_sigma, Msq_sigma, C);
for(int j3=1;j3<=length;j3++){
E_av=E_av+E[j3]; // sum of energy per site
Msq_av=Msq_av+Msq[j3]; // sum of magnetization squared per spin
rho_sf_av=rho_sf_av+rho_sf[j3];
C_av=C_av+C[j3];
}
E_av=(E_av/length);
Msq_av=(Msq_av/length);
rho_sf_av=(rho_sf_av/length);
C_av=(C_av/length);
// Calculate 1st sigma
long double tempE1=0;
long double tempE2=0;
long double tempMsq1=0;
long double tempMsq2=0;
long double tempW1=0;
long double tempW2=0;
long double tempC1=0;
long double tempC2=0;

for(int j4=1;j4<=length;j4++){
tempE1=tempE1+(E[j4]-E_av)*(E[j4]-E_av); // variance of E
tempMsq1=tempMsq1+(Msq[j4]-Msq_av)*(Msq[j4]-Msq_av); // variance of Msq
tempW1=tempW1+(rho_sf[j4]-rho_sf_av)*(rho_sf[j4]-rho_sf_av);
tempC1=tempC1+(C[j4]-C_av)*(C[j4]-C_av);
}
E_sigma=sqrt(tempE1)/length;
Msq_sigma=sqrt(tempMsq1)/length;
rho_sf_sigma=sqrt(tempW1)/length;
C_sigma=sqrt(tempC1)/length;
// Binning method
long double C_sigma_cur=C_sigma+2*delta; // current error of C to compare with previous error
long double* E_cur = new long double[(length+1)]; // E array for binning
long double* Msq_cur = new long double[(length+1)]; // Msq array for binning
long double* rho_sf_cur = new long double[(length+1)];
long double* C_temp = new long double[(length+1)];
for(int j5=1;j5<=length;j5++){
E_cur[j5]=E[j5];
Msq_cur[j5]=Msq[j5];
rho_sf_cur[j5]=rho_sf[j5];
C_temp[j5]=C[j5];
}
int k0=0; // group time
int k1;
int m0; // number of data points
while((C_sigma_cur-C_sigma)>delta||(C_sigma_cur-C_sigma)<(-delta)){
C_sigma=C_sigma_cur;
k0=k0+1;
k1=(pow(2, k0));
m0=(length/k1);
//if((m0%2)==1){
	//printf("False, length should be a larger power of 2.\n\n");
	//break;
//}
tempE2=0;
tempMsq2=0;
tempW2=0;
tempC2=0;
for(int k1=1;k1<=m0;k1++){
E_cur[k1]=(E_cur[(2*k1-1)]+E_cur[(2*k1)])/2;
Msq_cur[k1]=(Msq_cur[(2*k1-1)]+Msq_cur[(2*k1)])/2;
rho_sf_cur[k1]=(rho_sf_cur[(2*k1-1)]+rho_sf_cur[(2*k1)])/2;
C_temp[k1]=(C_temp[(2*k1-1)]+C_temp[(2*k1)])/2;
tempE2=tempE2+(E_cur[k1]-E_av)*(E_cur[k1]-E_av);
tempMsq2=tempMsq2+(Msq_cur[k1]-Msq_av)*(Msq_cur[k1]-Msq_av);
tempW2=tempW2+(rho_sf_cur[k1]-rho_sf_av)*(rho_sf_cur[k1]-rho_sf_av);
tempC2=tempC2+(C_temp[k1]-C_av)*(C_temp[k1]-C_av);	
}
C_sigma_cur=sqrt(tempC2)/m0;
Msq_sigma=sqrt(tempMsq2)/m0;
rho_sf_sigma=sqrt(tempW2)/m0;
E_sigma=sqrt(tempE2)/m0;
}

//output data to file
for(int i2=1;i2<=length;i2++){
DATA<<E[i2]<< //energy per site
"\t      "<<Msq[i2]<<endl; // Msq per spin array
}

double X=(Msq_av/T); //susceptibility per spin
double X_sigma=(Msq_sigma/T);
// output averages & errors & C & susceptibility
DATA<<E_av<< 
"\t    "<<Msq_av<<
"\t    "<<E_sigma<<
"\t    "<<Msq_sigma<<
"\t    "<<X<<         //susceptibility per spin
"\t      "<<X_sigma<<endl;


for(int h1=1;h1<=length;h1++){
CI<<C[h1]<<endl; // specific heat
}
CI<<C_av<< 
"\t      "<<C_sigma<<endl;

for(int a8=1;a8<=length;a8++){
SF<<rho_sf[a8]<<endl;
}
SF<<rho_sf_av<< 
"\t      "<<rho_sf_sigma<<endl;

printf(" d = %d\n L = %d\n T = %f\n\n", d, L, T);
printf("E per site = %Lf +/- %Lf\n", E_av, E_sigma);
printf("Msq_av = %Lf +/- %Lf\n", Msq_av, Msq_sigma);
printf("C_av = %Lf +/- %Lf\n", C_av, C_sigma);
printf("rho_sf = %Lf +/- %Lf\n", rho_sf_av, rho_sf_sigma);
printf("Susceptibility per spin = %f +/- %f\n\n", X, X_sigma);
printf("G[1] = %Lf\n", G[1]);
printf("Z = %Lf\n", Z);
printf("it = %Ld\n", it);

FILE *out_file = fopen("output.dat","a");

for(int i=1;i<=(2*d);i++){
for(int j=1;j<=n;j++){

fprintf(out_file, "%d", bond[i][j]);
  
printf("%d", bond[i][j]);
}	
}
fprintf(out_file, "\n");
fclose(out_file);



std::remove("bond_cur");



}
    return EXIT_SUCCESS;

}
}