#ifndef __CFILE_H
#define __CFILE_H
#include <fstream>

/* Bufferek hossza legalabb ennyi legyen */

using std::ofstream;
using std::ifstream;

#define K_BHOSSZ 80

char* strip(char * ezt,int =K_BHOSSZ);
char* spacekill(char *ezt,int =K_BHOSSZ);

#define getint(valt,def) int valt=def; kapcs.Getval(#valt,valt)
#define getdbl(valt,def) double valt=def; kapcs.Getval(#valt,valt)

#define Kgetint(valt,def) int valt=def; kapcs.Getval(#valt,valt)
#define Kgetdbl(valt,def) double valt=def; kapcs.Getval(#valt,valt)


class Konfitem
{
 public:
  Konfitem();
  Konfitem(const char*,const char*);
  Konfitem(char*,double);
  Konfitem(char*,int);
  ~Konfitem();
   
  Konfitem* Getnext() const {return next;}
  Konfitem* Getlast();
  void Setnext(Konfitem *kov) { next=kov; }
  void Setname(const char*);
  char* Getname() {return nevem;}
  void Setval(const char*);
  void Setval(int);
  void Setval(double);
  char* Getval() {return ertek;}
  void print();
  Konfitem* Whereis(const char *); 
  void save(ofstream& fout);

 private:
  Konfitem* next;
  char nevem[K_BHOSSZ];
  char ertek[K_BHOSSZ];
};

class Konfig
{
 public:
  Konfig();
  Konfig(const char*,const char*);
  Konfig(char*,double);
  Konfig(char*,int);
  Konfig(char *);
  Konfig(int argc,char* argv[]);  
  //  Konfig(char*);
  ~Konfig();

  Konfitem* Getfirst() const {return first;}
  void Set(char*,char*);
  void Set(char*,double);
  void Set(char*,int);
  char * Getval(const char*);
  int Getval(const char * const ,char*) const;
  int Getval(const char * const ,double&) const;
  int Getval(const char * const ,long double&) const;
  int Getval(const char * const ,float&) const;
  int Getval(const char * const ,int&) const;
  void read_file(char*);
  Konfitem* Getlast();
  Konfitem* Getlast_slow();
  void printall();
  Konfitem* exists(const char*) const;
  int Exists(const char *);
  void save_file(char *);


 private:
 Konfitem* first;
 Konfitem* last;
 void Newitem(char*,char*); 
 void Newitem(char*,double);
 void Newitem(char*,int);
 char tempb[K_BHOSSZ];

};








#endif
