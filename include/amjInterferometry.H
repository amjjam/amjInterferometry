#ifndef _AMJINTERFEROMETRY_H_
#define _AMJINTERFEROMETRY_H_

#include <cstring>
#include <complex>
#include <vector>

#include <amjPacket.H>
#include <amjPacketRW.H>

class Phasor: public std::complex<float>{
public:
  Phasor(){};
  virtual ~Phasor(){};
  int memsize(){return 2*sizeof(float);}
  virtual void write(amjPacket &) const;
  virtual void read(amjPacket &);
  // int write(char *);
  //int read(char *);
private:
};


class Phasors: public amjPacketRW{
public:
  Phasors(std::vector<Phasor>::size_type n=0){_p.resize(n);};
  virtual ~Phasors(){};
  const std::string &name(){return _name;};
  void rename(std::string name){_name=name;};
  std::vector<Phasor>::size_type size(){return _p.size();};
  void resize(std::vector<Phasor>::size_type n){_p.resize(n);};
  Phasor &operator[](int i){return _p[i];};
  virtual void write(amjPacket &p) const;
  virtual void read(amjPacket &p);
private:
  std::string _name;
  std::vector<Phasor> _p;
};


class PhasorSets: public amjPacketRW{
public:
  PhasorSets(std::vector<Phasors>::size_type n=0){_s.resize(n);};
  virtual ~PhasorSets(){};
  std::vector<Phasors>::size_type size(){return _s.size();};
  void resize(std::vector<Phasors>::size_type n){_s.resize(n);};
  Phasors &operator[](int i){return _s[i];};
  int memsize();
  virtual void write(amjPacket &p) const;
  virtual void read(amjPacket &p);
private:
  std::vector<Phasors> _s;
};

#endif /* _AMJINTERFEROMETRY */
