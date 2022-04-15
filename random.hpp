/* random number generator defined as class */
class Random
{
 public:
        // constructors
        Random(void);
        ~Random(void);

        // overloaded operators
        REAL operator () (void);                  // return uniform deviate
        void operator () (string title);  // set seed
     
        // general accessors
        void flush (string title); 
        REAL gauss (const REAL sigma);          // return gaussian deviate


 private:
        ULONG i,j,k,m,n,iseed,ii,jj;
        ULONG ij,kl;
        USHORT init,seed;
        REAL s,t,c,cd,cm;
        REAL getran();
		double u[98];
		string runid;
};

/*   Methods for Random class */
Random::Random()
{
 init=0;
 ij=1802;
 kl=9373;
 i=1;
 while ( i <= 97)
 {
  u[i]=0;
  i++;
 }
}

Random::~Random()
{
}

void Random ::operator () (string title)
{
  fstream file,logfile;
  runid = title + ".iseed";
  file.open(runid.c_str(),ios::binary|ios::in);
  runid = title + ".rnd";
  logfile.open(runid.c_str(),ios::out);
  if(!file)
  {
	  seed = 0; // there is no restart file; start from scratch
	  if(init)
	  { 
		  logfile << " seed has already been set : ignoring " << endl;
	  }
	  init=1;
  } 
  else
  {
	  seed = 1; // there is a restart file; continue from where it ended
	  for(ii=1; ii<=97; ii++) 
	  {
		  file.read(reinterpret_cast<char*>(&u[ii]), sizeof(u[ii]));
	  }
		file.read(reinterpret_cast<char*>(&c), sizeof(c));
 		file.read(reinterpret_cast<char*>(&cd), sizeof(cd));
 		file.read(reinterpret_cast<char*>(&cm), sizeof(cm));
 		file.read(reinterpret_cast<char*>(&i), sizeof(i));
 		file.read(reinterpret_cast<char*>(&j), sizeof(j));
	  file.close();
	  logfile << "valid restart file found for random number generator " 
		  << endl;
  }
/*  starting from scratch, create file */
  if(!seed)
  {
	  iseed=0;
	  i=((ij/177)%177)+2;
	  j=(      ij%177)+2;
	  k=((kl/169)%178)+1;
	  m=(      kl%169)  ;
	  for(ii=1; ii<=97; ii++)
	  {
		  s=0.;
		  t=0.5;
		  for(jj=1; jj<=24; jj++)
		  {
			  n=((((j*i)%179)*k)%179);
			  i=j; j=k; k=n;
			  m=((53*m+1)%169);
			  if(((m*n)%64) >= 32) s+=t;
			  t*=0.5;
		  }
		  u[ii]=s;
	  }
	  c=362436./16777216.;
	  cd= 7654321./16777216.;
	  cm=16777213./16777216.;
	  i=97;
	  j=33;
	  logfile << "no restart file found for random number generator: "
		  <<  " starting from scratch " << endl;
  }
  logfile.close();
}

REAL Random::operator () (void)
{
 return getran();
}

#include "random/methods"