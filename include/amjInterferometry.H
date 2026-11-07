#ifndef _AMJINTERFEROMETRY_H_
#define _AMJINTERFEROMETRY_H_

#include <cstring>
#include <complex>
#include <vector>
#include <cmath>
#include <cstdint>

#include <iostream>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <functional>
#include <memory>

#include <amjTime.H>

namespace amjInterferometry{

  //
  // Wavelengths
  //
  
  template<typename T=float>
  class Wavelengths{
  public:
    Wavelengths(){};
    template<typename U>
    Wavelengths(const Wavelengths<U> &L){
      _L.resize(L.size());
      for(int i=0;i<L.size();i++)
	_L[i]=L[i];
    }
    ~Wavelengths(){};
    size_t size() const {return _L.size();};
    void resize(size_t n) {_L.resize(n);};
    void name(std::string &name){_name=name;};
    const std::string& name() const {return _name;};
    T& operator[](int i){return _L[i];};
    const T& operator[](int i) const {return _L[i];};
  private:
    std::string _name;
    std::vector<T> _L;
  };

  template<typename T=float>
  Wavelengths<T> makeWavelengths(int nL, T L0, T L1,
                                 const std::string &name=""){
    if(nL < 2 || !std::isfinite(L0) || !std::isfinite(L1) ||
       L0 <= 0 || L1 <= 0)
      throw std::invalid_argument("Need at least two positive finite wavelengths");
    Wavelengths<T> w;
    std::string label=name;
    w.name(label);
    w.resize(static_cast<size_t>(nL));
    const long double b0=1.0L/L0, b1=1.0L/L1;
    for(int i=0;i<nL;i++)
      w[i]=static_cast<T>(1.0L/(b0+(b1-b0)*i/(nL-1)));
    return w;
  }

  //
  // Phasors
  //
  
  template<typename T=float>
  class Phasors{
  public:
    Phasors(size_t n=0){_w.resize(n); _p.resize(n);};
    virtual ~Phasors(){};
    void name(const std::string &n){_name=n;}
    const std::string &name() const {return _name;};
    void time(const amjTime &t){_time=t;};
    const amjTime &time() const {return _time;};
    size_t size() const {return _w.size();};
    void resize(size_t n){_w.resize(n); _p.resize(n);};
    T &operator()(int i){return _w[i];};
    const T &operator()(int i) const {return _w[i];}; 
    std::complex<T>& operator[](int i){return _p[i];};
    const std::complex<T>& operator[](int i) const {return _p[i];};
    template<typename U>
    Phasors<T>& operator+=(const Phasors<U>& other);
    size_t memsize() const {return memsize1()+memsize2();};
    int write(uint8_t *) const;
    size_t memsize1() const {return 2*sizeof(uint32_t);};
    int read1(uint8_t *);
    size_t memsize2() const {return _name.length()+_time.size()
	+_p.size()*3*sizeof(T);};
    int read2(uint8_t *);
  protected:
    std::string _name;
    amjTime _time;
    std::vector<T> _w;
    std::vector<std::complex<T> > _p;
  };

  template<typename T>
  template<typename U>
  Phasors<T>& Phasors<T>::operator+=(const Phasors<U>& other){
    if (_w.size() != other.size())
      throw std::runtime_error("Phasors size mismatch");
    
    for (size_t i = 0; i < other.size(); ++i) {
      _p[i] += other[i];
    }
    
    return *this;
  }
  
  template<typename T>
  int Phasors<T>::write(uint8_t *pos) const{
    uint32_t s1=_name.length(),s2=_p.size();
    uint8_t *p=pos;
    memcpy(p,&s1,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&s2,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,_name.data(),_name.length()); p+=_name.length();
    p+=_time.write(p);
    memcpy(p,_w.data(),sizeof(T)*_w.size()); p+=sizeof(T)*_w.size();
    memcpy(p,_p.data(),2*sizeof(T)*_p.size()); p+=2*sizeof(T)*_p.size();
    return p-pos;
  }
  
  template<typename T>
  int Phasors<T>::read1(uint8_t *pos){
    uint32_t s;
    uint8_t *p=pos;
    memcpy(&s,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    _name.resize(s);
    memcpy(&s,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    _w.resize(s);
    _p.resize(s);    
    return sizeof(uint32_t)+s;
  }
  
  template<typename T>
  int Phasors<T>::read2(uint8_t *pos){
    uint8_t *p=pos;
    memcpy(_name.data(),p,_name.length()); p+=_name.length();
    p+=_time.read(p);
    memcpy(_w.data(),p,sizeof(T)*_w.size()); p+=sizeof(T)*_w.size();
    memcpy(_p.data(),p,2*sizeof(T)*_p.size()); p+=2*sizeof(T)*_p.size();
    return p-pos;
  }

  
  template<typename T=float>
  class PhasorsAverage: public Phasors<T>{
  protected:
    using Phasors<T>::_w;
    using Phasors<T>::_p;
  public:
    PhasorsAverage(){};
    virtual ~PhasorsAverage(){};
    virtual bool add(const Phasors<T> &p);
    virtual double Tavg()=0;
    virtual void Tavg(double)=0;
  };
  
  template<typename T=float>
  class PhasorsBoxcar: public PhasorsAverage<T>{
    using PhasorsAverage<T>::_time;
    using PhasorsAverage<T>::_name;
    using PhasorsAverage<T>::_w;
    using PhasorsAverage<T>::_p;
  public:
    PhasorsBoxcar(const Phasors<T> &p, double tavg=1);
    virtual ~PhasorsBoxcar(){};
    virtual double Tavg() override {return w;}
    virtual void Tavg(double) override;
    virtual bool add(const Phasors<T> &p) override;
  private:
    double w;
    size_t n;
    std::vector<std::vector<std::complex<T>>> _h;
    size_t index;
    bool update(const Phasors<T> &p);
    void update_average();
  };

  template<typename T>
  PhasorsBoxcar<T>::PhasorsBoxcar(const Phasors<T> &p, double tavg):n(1),index(0){
    _name=p.name();
    _time=p.time();
    _w.resize(p.size());
    for(size_t i=0;i<_w.size();i++)
      _w[i]=p(i);
    _p.resize(p.size());
    Tavg(tavg);
  }
  
  template<typename T>
  void PhasorsBoxcar<T>::Tavg(double t){
    w=t;
    if(w<1)
      w=1;
    size_t new_n=std::round(w);
    if(new_n==n)
      return;
    if (_h.size()==n) // history is full; index points to the oldest sample
      std::rotate(_h.begin(),_h.begin()+index,_h.end());
    index=0;
    if (_h.size()>new_n)
      _h.erase(_h.begin(),_h.end()-new_n); // discard oldest, keep newest
    n=new_n;
    update_average();
  }
  
  template<typename T>
  bool PhasorsBoxcar<T>::add(const Phasors<T> &p){
    if(_h.size()<n){
      _h.emplace_back(_w.size());
      index=_h.size()-1;
    }
    for(size_t i=0;i<_w.size();i++)
      _h[index][i]=p[i];
    index=(index+1)%_h.size();
    update_average();
    _time=p.time();
    return true;
  };

  template<typename T>
  void PhasorsBoxcar<T>::update_average(){
    if(_h.size()==0)
      return;
    for(size_t i=0;i<_p.size();i++){
      _p[i]=0;
      for(size_t j=0;j<_h.size();j++)
	_p[i]+=_h[j][i];
      _p[i]/=_h.size();
    }
  }

  template<typename T>
  class PhasorsLeakyBox:public PhasorsAverage<T>{
    using PhasorsAverage<T>::_time;
    using PhasorsAverage<T>::_name;
    using PhasorsAverage<T>::_w;
    using PhasorsAverage<T>::_p;
  public:
    PhasorsLeakyBox(double w){Tavg(w);}
    virtual ~PhasorsLeakyBox(){};
    virtual bool add(const Phasors<T> &p) override;
    virtual double Tavg() override {return w;}
    virtual void Tavg(double t) override {w=t; if(w<1) w=1;}
  private:
    double w;
    bool update(const Phasors<T> &p);
  };
  
  template<typename T>
  bool PhasorsLeakyBox<T>::add(const Phasors<T> &p){
    double alpha=1.0/w;
    alpha=alpha<1?1:alpha;
    for(size_t i=0;i<_w.size();i++)
      _p[i]=p[i]*alpha+_p[i]*(1-alpha);
    _time=p.time();
    return true;
  }
  
  //
  // PowerSpectrum
  //
  
  template<typename T=float>
  class PowerSpectrum{
  public:
    PowerSpectrum(){}
    virtual ~PowerSpectrum(){}
    void name(std::string n){_name=n;}
    const std::string &name() const {return _name;}
    void time(const amjTime &t){_time=t;};
    const amjTime &time() const {return _time;};
    void clear(){_name.clear(); _delay.resize(0); _power.resize(0);}
    T &operator[](int i){return _power[i];}
    const T &operator[](int i) const {return _power[i];}
    T &operator()(int i){return _delay[i];}
    const T &operator()(int i) const {return _delay[i];}
    const std::vector<T> &delay(){return _delay;}
    const std::vector<T> &power(){return _power;}
    size_t size() const {return _delay.size();}
    void resize(size_t s){_delay.resize(s); _power.resize(s);}
    size_t memsize() const {return memsize1()+memsize2();}
    int write(uint8_t *) const;
    size_t memsize1() const {return 2*sizeof(uint32_t);}
    int read1(uint8_t *);
    size_t memsize2() const {return _name.size()+_time.size()
	+sizeof(T)*2*_delay.size();}
    int read2(uint8_t *);
    // Explicit configuration-time check; operator+= assumes matching grids.
    template<typename U>
    bool same_delay_grid(const PowerSpectrum<U> &ps) const;
    // Empty destinations adopt the first spectrum, including name/time.
    // Later additions preserve destination coordinates/name and copy incoming time.
    template<typename U=float>
    void operator+=(const PowerSpectrum<U> &ps);
    template<typename U=float>
    void operator/=(U);
  protected:
    std::string _name;
    amjTime _time;
    std::vector<T> _delay,_power;
  };

  template<typename T>
  int PowerSpectrum<T>::write(uint8_t *pp) const{
    uint32_t s1=_name.size(),s2=_delay.size();
    uint8_t *p=pp;
    memcpy(p,&s1,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&s2,sizeof(uint32_t)); p+=sizeof(T);
    memcpy(p,_name.data(),_name.size()); p+=_name.size();
    p+=_time.write(p);
    memcpy(p,_delay.data(),sizeof(T)*_delay.size()); p+=sizeof(T)*_delay.size();
    memcpy(p,_power.data(),sizeof(T)*_power.size()); p+=sizeof(T)*_power.size();
    return p-pp;
  }

  template<typename T>
  int PowerSpectrum<T>::read1(uint8_t *pp){
    uint8_t *p=pp;
    uint32_t s;    
    memcpy(&s,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    _name.resize(s);
    memcpy(&s,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    _delay.resize(s);
    _power.resize(s);
    return p-pp;
  }

  template<typename T>
  int PowerSpectrum<T>::read2(uint8_t *pp){
    uint8_t *p=pp;
    memcpy(_name.data(),p,_name.size()); p+=_name.size();
    p+=_time.read(p);
    memcpy(_delay.data(),p,sizeof(T)*_delay.size()); p+=sizeof(T)*_delay.size();
    memcpy(_power.data(),p,sizeof(T)*_power.size()); p+=sizeof(T)*_power.size();
    return p-pp;
  }
  
  template<typename T>
  template<typename U>
  bool PowerSpectrum<T>::same_delay_grid(const PowerSpectrum<U> &ps) const{
    if(size()!=ps.size()) return false;
    for(size_t i=0;i<size();i++)
      if(_delay[i]!=ps(i)) return false;
    return true;
  }

  template<typename T>
  template<typename U>
  void PowerSpectrum<T>::operator+=(const PowerSpectrum<U> &ps){
    if(size()==0){
      if(ps.size()==0) return;
      // Reuse capacity retained by clear(); do not allocate temporary vectors.
      resize(ps.size());
      for(size_t i=0;i<ps.size();i++){
        _delay[i]=static_cast<T>(ps(i));
        _power[i]=static_cast<T>(ps[i]);
      }
      _name=ps.name();
      _time=ps.time();
      return;
    }
    if(size()!=ps.size())
      throw std::invalid_argument("PowerSpectrum size mismatch");
    // Stable grids are a caller contract: no coordinate scan or allocation here.
    for(size_t i=0;i<size();i++)
      _power[i]+=ps[i];
    _time=ps.time();
  }

  template<typename T>
  template<typename U>
  void PowerSpectrum<T>::operator/=(U f){
    for(unsigned int i=0;i<_power.size();i++)
      _power[i]/=f;
  }


  //
  // PowerSpectrumAverage
  //
  
  template<typename T=float>
  class PowerSpectrumAverage: public PowerSpectrum<T>{
  protected:
    using PowerSpectrum<T>::_name;
    using PowerSpectrum<T>::_time;
    using PowerSpectrum<T>::_delay;
    using PowerSpectrum<T>::_power;
  public:
    PowerSpectrumAverage():old_N(0),old_d0(0),old_d1(0){}
    virtual ~PowerSpectrumAverage(){}
    virtual bool add(const PowerSpectrum<T>&)=0;
    virtual double Tavg()=0;
    virtual void Tavg(double)=0;
    virtual void reset(){PowerSpectrum<T>::clear();_time=amjTime();old_N=0;}
    bool old(size_t& N,T& d0,T& d1) const {N=old_N;d0=old_d0;d1=old_d1;return old_N>0;}
  protected:
    bool update(const PowerSpectrum<T>& p){
      if(!p.size())throw std::invalid_argument("Empty averaging input");
      if(_delay.size()>0&&p.size()==_delay.size()&&p(0)==_delay[0]&&
         p(p.size()-1)==_delay.back()&&_name==p.name())return false;
      set_old();
      _delay.resize(p.size());_power.assign(p.size(),T(0));
      for(size_t i=0;i<p.size();++i)_delay[i]=p(i);
      _name=p.name();
      return true;
    }
    void set_old(){
      if(!_delay.size())return;
      old_N=_delay.size();old_d0=_delay[0];old_d1=_delay.back();
    }
  private:
    size_t old_N;
    T old_d0,old_d1;
  };

  //
  // PowerSpectrumBoxcar
  //
  
  template<typename T=float>
  class PowerSpectrumBoxcar: public PowerSpectrumAverage<T>{
    using PowerSpectrumAverage<T>::_time;
    using PowerSpectrumAverage<T>::_power;
    using PowerSpectrumAverage<T>::update;
  public:
    PowerSpectrumBoxcar(double w=1):n(1),index(0){Tavg(w);}
    virtual ~PowerSpectrumBoxcar(){};
    double Tavg() override {return w;}
    void Tavg(double);
    bool add(const PowerSpectrum<T>&) override;
    void reset() override {PowerSpectrumAverage<T>::reset();_h.clear();index=0;}
  private:
    double w;
    size_t n,index;
    std::vector<std::vector<T>> _h;
    void update_average();
  };

  template<typename T>
  void PowerSpectrumBoxcar<T>::Tavg(double t){
    if(!std::isfinite(t))throw std::invalid_argument("Boxcar width must be finite");
    w=std::max(1.0,t);
    size_t new_n=std::round(w);
    if(new_n==n)return;
    if(_h.size()==n)std::rotate(_h.begin(),_h.begin()+index,_h.end());
    index=0;
    if(_h.size()>new_n)_h.erase(_h.begin(),_h.end()-new_n);
    n=new_n;update_average();
  }

  template<typename T>
  bool PowerSpectrumBoxcar<T>::add(const PowerSpectrum<T>& p){
    bool changed=update(p);
    if(changed){_h.clear();index=0;}
    if(_h.size()<n){_h.emplace_back(p.size());index=_h.size()-1;}
    for(size_t i=0;i<p.size();++i)_h[index][i]=p[i];
    index=(index+1)%_h.size();
    update_average();_time=p.time();
    return changed;
  }

  template<typename T>
  void PowerSpectrumBoxcar<T>::update_average(){
    if(!_h.size())return;
    for(size_t i=0;i<_power.size();++i){
      _power[i]=0;
      for(size_t j=0;j<_h.size();++j)_power[i]+=_h[j][i];
      _power[i]/=_h.size();
    }
  }

  //
  // PowerSpectrumLeakyBox
  //
  
  template<typename T=float>
  class PowerSpectrumLeakyBox: public PowerSpectrumAverage<T>{
    using PowerSpectrumAverage<T>::_time;
    using PowerSpectrumAverage<T>::_power;
    using PowerSpectrumAverage<T>::update;
  public:
    PowerSpectrumLeakyBox(double w=1){Tavg(w);}
    virtual ~PowerSpectrumLeakyBox(){};
    double Tavg() override {return w;}
    void Tavg(double);
    bool add(const PowerSpectrum<T>&) override;
  private:
    double w;
  };

  template<typename T>
  void PowerSpectrumLeakyBox<T>::Tavg(double t){
    if(!std::isfinite(t))throw std::invalid_argument("Leaky-box width must be finite");
    w=std::max(1.0,t);
  }

  template<typename T>
  bool PowerSpectrumLeakyBox<T>::add(const PowerSpectrum<T>& p){
    bool changed=update(p);
    double alpha=1.0/w;
    if(changed)alpha=1;
    for(size_t i=0;i<_power.size();++i)
      _power[i]=alpha*p[i]+(1-alpha)*_power[i];
    _time=p.time();
    return changed;
  }
  
  //
  // CalcPowerSpectrum
  //
  
  template<typename T=float>
  class CalcPowerSpectrum{
  public:
    CalcPowerSpectrum(){};
    virtual ~CalcPowerSpectrum()=0;
    template<typename U, typename V>
    void operator()(const Phasors<U> &, PowerSpectrum<V> &);
    template<typename U, typename V>
    void operator()(const std::vector<Phasors<U> > &,
		    std::vector<PowerSpectrum<V> > &);
  private:

  };

  template<typename T>
  CalcPowerSpectrum<T>::~CalcPowerSpectrum() {}
  
  //
  // CalcPowerSpectrumFFT
  //
  
  void four1(float data[], unsigned long nn, int isign);
  void four1(double data[], unsigned long nn, int isign);

  template<typename T=float>
  class CalcPowerSpectrumFFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumFFT():wfirst(0),wlast(0),nOver(1),
      validateWavelengths(false),validatedSize(0){};
    virtual ~CalcPowerSpectrumFFT(){};
    void set_nover(int n){
      if(n < 1) throw std::invalid_argument("FFT oversampling must be positive");
      nOver=n;
    };
    // Assume a stable channel grid. When enabled, validate once, then cache.
    // Calling this setter again (even with true) forces a fresh check.
    void set_validate_wavelengths(bool enabled){
      validateWavelengths=enabled;
      validatedSize=0;
    }
    void set_wavelength(T, T);
    template<typename U>
    void set_wavelength(const std::vector<std::pair<U,U> > &w);
    template<typename U, typename V>
    void operator()(const Phasors<U> &, PowerSpectrum<V> &);
    template<typename U, typename V>
    void operator()(const std::vector<Phasors<U> > &,
			    std::vector<PowerSpectrum<V> > &);
  private:
    T wfirst,wlast;
    std::vector<std::pair<T,T> > wavelengths;
    int nOver;
    std::vector<T> fftData;
    bool validateWavelengths;
    size_t validatedSize;
  };

  template<typename T>
  void CalcPowerSpectrumFFT<T>::set_wavelength(T first, T last){
    if(!std::isfinite(first) || !std::isfinite(last) ||
       first <= 0 || last <= 0 || first == last)
      throw std::invalid_argument("FFT needs distinct positive wavelength endpoints");
    wfirst=first;
    wlast=last;
    wavelengths.clear();
    validatedSize=0;
  }

  template<typename T>
  template<typename U>
  void CalcPowerSpectrumFFT<T>::set_wavelength(
      const std::vector<std::pair<U,U> > &){
    // Channel-boundary semantics/resampling are not defined in the old code.
    // Fail here rather than silently ignoring this configuration.
    throw std::logic_error("FFT wavelength-pair configuration is not implemented");
  }

  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumFFT<T>::operator()(const Phasors<U> &phasors,
                                          PowerSpectrum<V> &pspec){
    const size_t n=phasors.size();
    if(n < 2) throw std::invalid_argument("FFT needs at least two phasors");
    if(wfirst <= 0 || wlast <= 0)
      throw std::logic_error("Call FFT set_wavelength(first,last) before calculation");
    const double k0=1.0/static_cast<double>(wfirst);
    const double dk=(1.0/static_cast<double>(wlast)-k0)/(n-1);
    if(validateWavelengths && validatedSize != n){
      // Allow rounding of wavelength coordinates to their declared precision.
      const long double eps=std::max(
          static_cast<long double>(std::numeric_limits<T>::epsilon()),
          static_cast<long double>(std::numeric_limits<U>::epsilon()));
      for(size_t i=0;i<n;i++){
        if(!std::isfinite(phasors(i)) || phasors(i) <= 0)
          throw std::invalid_argument("FFT wavelengths must be positive and finite");
        const long double expected=k0+i*dk;
        if(std::abs(1.0L/phasors(i)-expected) > 8*eps*std::abs(expected))
          throw std::invalid_argument("FFT requires uniform wavenumbers matching its endpoints");
      }
      validatedSize=n; // Cache only a successful check of the entire grid.
    }

    const size_t limit=std::min(fftData.max_size()/2-1,
        static_cast<size_t>(std::numeric_limits<unsigned long>::max()/2));
    if(n > limit/static_cast<size_t>(nOver))
      throw std::length_error("FFT size overflow");
    const size_t required=n*static_cast<size_t>(nOver);
    size_t nFFT=1; // Number of complex samples, including zero padding.
    while(nFFT < required){
      if(nFFT > limit/2) throw std::length_error("FFT size overflow");
      nFFT*=2;
    }
    // four1 is one-based: provide a real dummy element, never data()-1.
    fftData.assign(2*nFFT+1,T(0));
    for(size_t i=0;i<n;i++){
      fftData[2*i+1]=static_cast<T>(phasors[i].real());
      fftData[2*i+2]=static_cast<T>(phasors[i].imag());
    }
    four1(fftData.data(),static_cast<unsigned long>(nFFT),1);
    pspec.resize(nFFT);
    pspec.name(phasors.name());
    pspec.time(phasors.time());
    const double delayStep=1.0/(nFFT*dk);
    for(size_t i=0;i<nFFT;i++){
      // nFFT is a power of two, so wrapping needs no integer division.
      const size_t bin=(i+nFFT/2)&(nFFT-1);
      const double signedBin=static_cast<double>(i)-nFFT/2;
      const T re=fftData[2*bin+1], im=fftData[2*bin+2];
      pspec[i]=static_cast<V>(re*re+im*im);
      pspec(i)=static_cast<V>(signedBin*delayStep);
    }
    // Keep delay coordinates increasing for either wavelength ordering.
    if(dk < 0){
      for(size_t i=0;i<nFFT/2;i++){
        std::swap(pspec[i],pspec[nFFT-1-i]);
        std::swap(pspec(i),pspec(nFFT-1-i));
      }
    }
  }

  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumFFT<T>::operator()(const std::vector<Phasors<U> > &ph,
                                          std::vector<PowerSpectrum<V> > &ps){
    ps.resize(ph.size());
    for(size_t i=0;i<ph.size();i++) operator()(ph[i],ps[i]);
  }
  
  //
  // CalcPowerSpectrumDFT
  //
  
  template<typename T=float>
  class CalcPowerSpectrumDFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumDFT(){};
    ~CalcPowerSpectrumDFT(){};
    template<typename U>
    void set_delay(const std::vector<U> &);
    //template<typename U>
    //void set_delay(const std::vector<std::vector<U> > &);
    template<typename U, typename V>
    void operator()(const Phasors<U> &, PowerSpectrum<V> &);
    template<typename U, typename V>
    void operator()(const std::vector<Phasors<U> > &,
			    std::vector<PowerSpectrum<V> > &);
  private:
    std::vector<T> delay;
  };
    
  template<typename T>
  template<typename U>
  void CalcPowerSpectrumDFT<T>::set_delay(const std::vector<U> &d){
    if(d.empty()) throw std::invalid_argument("DFT delay grid must not be empty");
    std::vector<T> converted(d.size());
    for(size_t i=0;i<d.size();i++){
      converted[i]=static_cast<T>(d[i]);
      if(!std::isfinite(converted[i]))
        throw std::invalid_argument("DFT delays must be finite");
    }
    delay.swap(converted);
  }

  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumDFT<T>::operator()(const Phasors<U> &phasors,
                                          PowerSpectrum<V> &pspec){
    if(delay.empty()) throw std::logic_error("Configure the DFT delay grid first");
    if(phasors.size()==0) throw std::invalid_argument("DFT input must not be empty");
    // Real-time path: callers provide finite phasors and positive wavelengths.
    pspec.name(phasors.name());
    pspec.time(phasors.time());
    pspec.resize(delay.size());
    const double twoPi=2*std::acos(-1.0);
    for(size_t i=0;i<delay.size();i++){
      pspec(i)=static_cast<V>(delay[i]);
      // Retain the original arithmetic: double phase/trigonometry,
      // input-precision rotations, and T-precision accumulation.
      std::complex<T> sum(0,0);
      std::complex<U> rot;
      for(size_t j=0;j<phasors.size();j++){
        rot=std::polar(1.0,twoPi*delay[i]/phasors(j));
        sum+=phasors[j]*rot;
      }
      pspec[i]=static_cast<V>(std::norm(sum));
    }
  }

  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumDFT<T>::operator()(const std::vector<Phasors<U> > &ph,
					   std::vector<PowerSpectrum<V> > &ps){
    ps.resize(ph.size());
    for(size_t i=0;i<ph.size();i++)
      operator()(ph[i],ps[i]);
  }

  //
  // AnalyzePowerSpectrum
  //

  template<typename T=float>
  class AnalyzePowerSpectrum{
  public:
    AnalyzePowerSpectrum():peak_width(40){};
    ~AnalyzePowerSpectrum(){};
    template<typename U, typename V>
    void operator()(const PowerSpectrum<U> &, V &SNR, V &delay);
    template<typename U>
    void set_peak_width(U w){peak_width=w;}
  private:
    T peak_width;
  };
  
  template<typename T>
  template<typename U, typename V>
  void AnalyzePowerSpectrum<T>::operator()(const PowerSpectrum<U> &p, V &snr,
					   V &delay){
    int imax=-1;
    T peak=-1e31;
    
    for(size_t i=0;i<p.size();i++)
      if(p[i]>peak){
	imax=i;
	peak=p[i];
      }
    
    delay=p(imax);
    
    int count=0;
    T v=0,v2=0;
    for(size_t i=0;i<p.size();i++)
      if(p(i)<delay-peak_width/2||p(i)>delay+peak_width/2){
	count++;
	v+=p[i];
	v2+=p[i]*p[i];
      }

    T mean=v/count;
    
    snr=(peak-mean)/sqrt(v2/count-mean*mean);      
  }
  
  //
  // PowerSpectrum smoothing
  //
  
  // Gaussian
  // 1/\sqrt{2pi}/sigma*exp(-(x-mu)^2/2\sigma^2)
  inline double gaussian(double x, double m, double s){
    return 1.0/sqrt(2*M_PI)/s*exp(-(x-m)*(x-m)/2.0/s/s);
  }
  
  // Gaussian smoothing
  template<typename T, typename U>
  PowerSpectrum<T> gsmooth(const PowerSpectrum<T> &p, U sigma){
    if(!std::isfinite(sigma) || sigma <= 0)
      throw std::invalid_argument("Gaussian width must be positive and finite");
    PowerSpectrum<T> result=p;
    for(size_t i=0;i<p.size();i++){
      long double weighted=0, weights=0;
      for(size_t j=0;j<p.size();j++){
        const long double x=(static_cast<long double>(p(j))-p(i))/sigma;
        const long double weight=std::exp(-x*x/2);
        weighted+=weight*p[j];
        weights+=weight;
      }
      result[i]=static_cast<T>(weighted/weights);
    }
    return result;
  }

  // Boxcar smoothing
  template<typename T, typename U>
  PowerSpectrum<T> bsmooth(const PowerSpectrum<T> &p, U w){
    PowerSpectrum<T> _p=p;
    double g,s;
    for(size_t i=0;i<p.size();i++){
      g=0;
      s=0;
      for(size_t j=0;j<p.size();j++){
	if(fabs(p(j)-p(i))<w/2){
	  g++;
	  s+=p[j];
	}
      }
      _p[i]=s/g;
    }
    return _p;
  }


  //
  // DelayMachine
  //
  
  // Widths: sample counts for time/output; delay units for the delay axis.
  // struct DelayMachineConfig {
  //   enum class TemporalMethod { None, Boxcar, Leaky };
  //   enum class AxisMethod { None, Boxcar, Gaussian };
  //   // These select the internally owned default averager.
  //   TemporalMethod time_method=TemporalMethod::Boxcar;
  //   unsigned int time_width=5;
  //   AxisMethod delay_axis_method=AxisMethod::None;
  //   float delay_axis_width=0;
  //   float peak_exclusion_width=40;
  //   TemporalMethod output_method=TemporalMethod::None;
  //   unsigned int output_width=1;
  //   bool operator==(const DelayMachineConfig&) const;
  //   bool operator!=(const DelayMachineConfig& c) const {return !(*this==c);}
  //   void validate() const;
  // };

  // class DelayMachine {
  // public:
  //   using Config=DelayMachineConfig;
  //   using Phasors=amjInterferometry::Phasors<float>;
  //   using PowerSpectrum=amjInterferometry::PowerSpectrum<float>;
  //   using Averager=PowerSpectrumAverage<float>;
  //   using Calculator=std::function<void(const Phasors&,PowerSpectrum&)>;
  //   using Analyzer=std::function<void(const PowerSpectrum&,float&,float&)>;
  //   // Compatibility overload preserves the old zero-filled output startup.
  //   DelayMachine(std::string,Calculator,int nIncoherent,int nSmooth,
  // 		 Analyzer analyzer=AnalyzePowerSpectrum<float>());
  //   DelayMachine(std::string,Calculator,const Config& config=Config());
  //   const std::string& name() const {return _name;}
  //   void load(const Phasors&);
  //   float delay() const {return _delay;}
  //   float snr() const {return _snr;}
  //   const amjTime& time() const {return pSpec.time();}
  //   const PowerSpectrum& p() const {return pSpec;}
  //   const Config& config() const {return config_;}
  //   void set_config(const Config&);
  //   void set_calculator(Calculator);
  //   // The caller owns an injected averager and configures its mode and width.
  //   void set_averager(Averager&);
  //   void use_default_averager();
  //   // Custom analyzers replace the default; peak_exclusion_width then has no effect.
  //   void set_analyzer(Analyzer);
  //   void use_default_analyzer();
  //   void reset();
  // private:
  //   void configure_default_averager();
  //   void build_axis_weights();
  //   std::string _name;
  //   Config config_;
  //   Calculator calculateSpectrum;
  //   Analyzer analyzeSpectrum;
  //   bool customAnalyzer_=false, legacyOutputStartup_=false;
  //   PowerSpectrum raw_,pSpec;
  //   std::unique_ptr<Averager> ownedAverager_;
  //   Averager* averager_=nullptr;
  //   std::vector<float> peakHistory,snrHistory;
  //   size_t nextOutput_=0,outputCount_=0;
  //   float outputDecay_=0,_delay=0,_snr=0;
  //   // Sparse normalized kernel, built once per configuration/grid size.
  //   std::vector<size_t> axisOffsets_,axisIndices_;
  //   std::vector<double> axisWeights_;
  // };
  
  //
  // Obsolete
  //
  
  template<typename T>
  int findMax(std::vector<T> &p);
  
  template<typename T>
  int findMax(std::vector<T> &p){
    int j=p.size()/2;
    T max=0;
    
    for(int i=0;i<p.size();i++)
      if(p[i]>max){
	j=i;
	max=p[i];
      }
    
    return j;
  }
  
}

#endif /* _AMJINTERFEROMETRY */
