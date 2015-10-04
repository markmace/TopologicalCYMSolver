#include <iostream>
#include <fstream>

#include "cfile.h"

using std::cout;
using std::cerr;
using std::endl;

Konfig::Konfig()
{
  first = new Konfitem();
}

Konfig::Konfig(const char* nev ,const char* ertek)
{
  first = new Konfitem(nev,ertek);
  last=first;
}

Konfig::Konfig(char* nev ,double ertek)
{
  first = new Konfitem(nev,ertek);
  last=first;
}

Konfig::Konfig(char* nev ,int ertek)
{
  first = new Konfitem(nev,ertek);
  last=first;
}

Konfig::Konfig(char *filenev)
{
  first=new Konfitem("Config ","file");
  last=first;
  read_file(filenev);
}

Konfig::Konfig(int argc,char * argv[])
{
  char buffer[K_BHOSSZ];
  first =new Konfitem("__Command","line");
  last=first;

  if (argc==1) { first =new Konfitem("Config ","file"); last=first; return; }

  int i,j;

  for(i=1;i<argc;i++)
    if (argv[i][0]=='-') 
      { for(j=0; (j<K_BHOSSZ) && (argv[i][j+1]!=0 );j++) 
        	buffer[j]=argv[i][j+1];
        buffer[j]=0;
	if ((i+1)<argc) Set(buffer,argv[i+1]);      
        	else Set(buffer,buffer);
      }
}

void Konfig::read_file(char *filenev)
{
  ifstream fin(filenev);
  char buffer[K_BHOSSZ],nev[K_BHOSSZ],ertek[K_BHOSSZ];

  cerr << "Reading configfile: " << filenev << endl;
  
   int error=0;
  int sor=0;
  int i,j;

 
     do { fin.getline(buffer,K_BHOSSZ); sor++;
     i=0;
     while  ((buffer[i]==' ' || buffer[i]=='\t') && i<K_BHOSSZ) i++;
     if (buffer[i]=='#' || buffer[i]==0) continue; 


 for(i=0;i<K_BHOSSZ && buffer[i]!= '=' ;i++) ;
 if (i == K_BHOSSZ) { 
   error=1; cerr << "Error in params file line " << sor 
               << ": No equality sign"<< endl; break; }

  for(j=0;j<i;j++) nev[j]=buffer[j]; nev[i]= '\0';

  for(j=i+1;buffer[j]!='\0';j++) ertek[j-i-1]=buffer[j];
  ertek[j-i-1]=0;   

  //  cerr << nev << " " << ertek << endl;

  Set(nev,ertek);
 
    }   while ( !fin.eof() && error == 0 );



  fin.close();
}

void Konfig::save_file(char* filenev)
{
  ofstream fout(filenev);
  first->save(fout);
  fout.close();
}

void Konfitem::save(ofstream& fout)
{
  fout << nevem << "=" << ertek << endl;
  if (next!=NULL) next->save(fout);
}

Konfig::~Konfig()
{
  delete first;
}

Konfitem* Konfig::Getlast()
{
  return first->Getlast();
}

Konfitem* Konfig::Getlast_slow()
{
  return first->Getlast();
}

void Konfig::Newitem(char* nev,char* ert)
{ 
 last->Setnext(new Konfitem(nev,ert));
 last=last->Getnext();
}

void Konfig::Newitem(char* nev,double ert)
{
 last->Setnext(new Konfitem(nev,ert));
 last=last->Getnext();
}

void Konfig::Newitem(char* nev,int ert)
{ 
 last->Setnext(new Konfitem(nev,ert));
 last=last->Getnext();
}

void Konfig::printall()
{
  first->print();
}

void Konfig::Set(char* nev,char* ertek)
{ Konfitem * mutato=exists(nev); 
 if (mutato==NULL) 
    {  
    strncpy(tempb,nev,K_BHOSSZ);
    strip(tempb);
    Newitem(tempb,ertek); return; 
    }
 mutato->Setval(ertek);
}

void Konfig::Set(char* nev,double ertek)
{ Konfitem * mutato=exists(nev); 
 if (mutato==NULL) 
 {  
    strncpy(tempb,nev,K_BHOSSZ);
    strip(tempb);
    Newitem(tempb,ertek); return; 
    }
 mutato->Setval(ertek);
}

void Konfig::Set(char* nev,int ertek)
{ Konfitem * mutato=exists(nev); 
 if (mutato==NULL) 
   {  
    strncpy(tempb,nev,K_BHOSSZ);
    strip(tempb);
    Newitem(tempb,ertek); return; 
    }
 mutato->Setval(ertek);
}

int Konfig::Exists(const char* nev)
{ Konfitem* ezaz=exists(nev);
 if (ezaz==NULL) return 0;
 return 1;
}


char* Konfig::Getval(const char* nev) 
{ Konfitem* ezaz=exists(nev);
 if (ezaz==NULL) return NULL;
 return ezaz->Getval();
}

int Konfig::Getval(const char *nev,char* ide) const
{ Konfitem *ezaz=exists(nev);
 if (ezaz==NULL) return 1;
 strncpy(ide,ezaz->Getval(),K_BHOSSZ);
 return 0;
}

int  Konfig::Getval(const char *nev,double& ide) const
{ Konfitem *ezaz=exists(nev);
    if (ezaz==NULL) return 1;
    ide=atof(ezaz->Getval());
    return 0;
}

int  Konfig::Getval(const char *nev,long double& ide) const
{ Konfitem *ezaz=exists(nev);
 if (ezaz==NULL) return 1;
    ide=atof(ezaz->Getval());
 return 0;
}

int  Konfig::Getval(const char *nev,float& ide) const
{ Konfitem *ezaz=exists(nev);
    if (ezaz==NULL) return 1;
    ide=atof(ezaz->Getval());
    return 0;
}

int  Konfig::Getval(const char *nev,int& ide) const
{ Konfitem *ezaz=exists(nev);
 if (ezaz==NULL) return 1;
 ide=atoi(ezaz->Getval());
 return 0;
}

Konfitem* Konfig::exists(const char* ez) const
{ char nev[K_BHOSSZ];
 
 strncpy(nev,ez,K_BHOSSZ);
 strip(nev);
 return first->Whereis(nev);
}


Konfitem* Konfitem::Whereis(const char *nev)
{ int i=0;

 while ( nev[i]==nevem[i] && nev[i]!=0 && nevem[i]!=0 ) i++;
 if (nev[i]==0 && nevem[i]==0 && i!=0) return this;
 if (next!=NULL) return next->Whereis(nev);
 return NULL;
}



void Konfitem::print()
{
  cerr << "nev:" << nevem << "<" << endl;
  cerr << "ertekem:" << ertek << "<" << endl;
  if (next!=NULL) next->print();
}

Konfitem::Konfitem()
{

  nevem[0]=0;
  ertek[0]=0;
  next=NULL;
}

Konfitem::Konfitem(const char* nev,const char* ert)
{  
  Setname(nev);
  Setval(ert);
  next=NULL;
}

Konfitem::Konfitem(char* nev,double ert)
{  
  Setname(nev);
  Setval(ert);
  next=NULL;
}

Konfitem::Konfitem(char* nev,int ert)
{  
  Setname(nev);
  Setval(ert);
  next=NULL;
}

Konfitem::~Konfitem()
{
  delete next;
}

void Konfitem::Setname(const char* erre)
{ 
  strncpy(nevem,erre,K_BHOSSZ);
  nevem[K_BHOSSZ-1]=0;
}

void Konfitem::Setval(const char* erre)
{ 
  strncpy(ertek,erre,K_BHOSSZ);
  ertek[K_BHOSSZ-1]=0;
}

void Konfitem::Setval(int szam)
{ 
  snprintf(ertek,K_BHOSSZ,"%d",szam);
}

void Konfitem::Setval(double szam)
{ 
  snprintf(ertek,K_BHOSSZ,"%g",szam);
}

Konfitem* Konfitem::Getlast()
{
  if (next==NULL) return this;
  return next->Getlast();
}

// kiszedi az elso spaceket, es a szo utani elsotol mindent elhagy
char* strip(char * ezt,int max)
{ int i=0,k=0;
 while ((ezt[i]==' ' || ezt[i]=='\t') && i< max) i++;
 for( ;i<max && ezt[i]!=0 && ezt[i]!=' ' ;i++,k++) 
          { ezt[k]=ezt[i]; } 
 ezt[k]=0;
 return ezt;
}

// kiszedi a sztringbol a spacekat
char* spacekill(char *ezt,int max)
{ int i=0,k=0;
 
 do {
 while (ezt[i]==' ') i++;
 for( ;i<max && ezt[i]!=0 && ezt[i]!=' ' ;i++,k++)  
    { ezt[k]=ezt[i]; }
 } while (i<max && ezt[i]!=0);
 ezt[k]=0;
 return ezt;
}
