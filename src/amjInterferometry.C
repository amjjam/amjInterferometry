#include "../include/amjInterferometry.H"

#include <iostream>

void Phasor::write(amjPacket &p) const{
  float _real=real();
  float _imag=imag();
  p << _real << _imag;
}

void Phasor::read(amjPacket &p){
  float _real;
  float _imag;
  p >> _real >> _imag;
  real(_real);
  imag(_imag);
}

void Phasors::write(amjPacket &p) const{
  p << _name << (int)_p.size();
  p.reserve(p.size()+_p.size()*sizeof(Phasor));
  for(unsigned int i=0;i<_p.size();i++)
    p << _p[i];
}

void Phasors::read(amjPacket &p){
  int n;
  p >> _name >> n;
  _p.resize(n);
  for(unsigned int i=0;i<_p.size();i++)
    p >> _p[i];
}

void PhasorSets::write(amjPacket &p) const{
  p << (int)_s.size();
  for(unsigned int i=0;i<_s.size();i++)
    p << _s[i];
}

void PhasorSets::read(amjPacket &p){
  int n;
  p >> n;
  _s.resize(n);
  for(unsigned int i=0;i<_s.size();i++)
    p >> _s[i];
}

