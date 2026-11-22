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
#include <format>
#include <numbers>

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
    Phasors(const Phasors &)=default;
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

  //
  // PhasorsAverage
  //
  
  template<typename T=float>
  class PhasorsAverage: public Phasors<T>{
  protected:
    using Phasors<T>::_w;
    using Phasors<T>::_p;
  public:
    PhasorsAverage()=default;
    explicit PhasorsAverage(const Phasors<T> &p):Phasors<T>(p){};
    virtual ~PhasorsAverage(){};
    virtual bool add(const Phasors<T> &p)=0;
    virtual double tavg()=0;
    virtual void tavg(double)=0;
  };
  
  //
  // PhasorsBoxcar
  //

  template<typename T=float>
  class PhasorsBoxcar: public PhasorsAverage<T>{
    using PhasorsAverage<T>::_time;
    using PhasorsAverage<T>::_name;
    using PhasorsAverage<T>::_w;
    using PhasorsAverage<T>::_p;
  public:
    PhasorsBoxcar(const Phasors<T> &p, double w=1):
      PhasorsAverage<T>(p),n(1),index(0){tavg(w);}
    virtual ~PhasorsBoxcar(){};
    virtual double tavg() override {return w;}
    virtual void tavg(double) override;
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
  void PhasorsBoxcar<T>::tavg(double t){
    w=std::max(1.0,t);
    const auto new_n=static_cast<std::size_t>(std::round(w));
    if(new_n==n)
      return;

    if(_h.size()==n && index!=0)
      std::rotate(_h.begin(),_h.begin()+index,_h.end());
    
    if(_h.size()>new_n)
      _h.erase(_h.begin(),_h.end()-new_n);
    
    n=new_n;
    index=0;
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
    if(_h.empty())
      return;
    for(size_t i=0;i<_p.size();i++){
      _p[i]=0;
      for(size_t j=0;j<_h.size();j++)
	_p[i]+=_h[j][i];
      _p[i]/=_h.size();
    }
  }

  //
  // PhasorsLeakyBox
  //
  
  template<typename T=float>
  class PhasorsLeakyBox:public PhasorsAverage<T>{
    using PhasorsAverage<T>::_time;
    using PhasorsAverage<T>::_name;
    using PhasorsAverage<T>::_w;
    using PhasorsAverage<T>::_p;
  public:
    PhasorsLeakyBox(const Phasors<T> &p, double w=1):
      PhasorsAverage<T>(p){tavg(w);}
    virtual ~PhasorsLeakyBox(){};
    virtual bool add(const Phasors<T> &p) override;
    virtual double tavg() override {return w;}
    virtual void tavg(double t) override {w=t; if(w<1) w=1;}
  private:
    double w;
    bool update(const Phasors<T> &p);
  };
  
  template<typename T>
  bool PhasorsLeakyBox<T>::add(const Phasors<T> &p){
    T alpha=1.0/w;
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
    bool update(const PowerSpectrum<T>&);
  public:
    PowerSpectrumAverage(){}
    virtual ~PowerSpectrumAverage(){}
    virtual PowerSpectrumAverage<T> &operator()(const PowerSpectrum<T>&)=0;
    virtual double tavg()=0;
    virtual void tavg(double)=0;
  };

  template<typename T>
  bool PowerSpectrumAverage<T>::update(const PowerSpectrum<T>& p){
    if(_name==p.name()&&this->same_delay_grid(p))
      return false;
    _delay.resize(p.size());
    for(size_t i=0;i<_delay.size();i++)
      _delay[i]=p(i);
    _power.resize(p.size());
    _name=p.name();
    return true;
  }

  //
  // PowerSpectrumBoxcar
  //
  
  template<typename T=float>
  class PowerSpectrumBoxcar: public PowerSpectrumAverage<T>{
    using PowerSpectrumAverage<T>::_name;
    using PowerSpectrumAverage<T>::_time;
    using PowerSpectrumAverage<T>::_delay;
    using PowerSpectrumAverage<T>::_power;
  public:
    PowerSpectrumBoxcar(double w=1){tavg(w);}
    virtual ~PowerSpectrumBoxcar(){};
    double tavg() override {return w;}
    void tavg(double) override;
    PowerSpectrumBoxcar<T> &operator()(const PowerSpectrum<T>&) override;
  private:
    double w;
    size_t n=1,index=0;
    std::vector<std::vector<T>> _h;
    void update_average();
    //    bool update(const PowerSpectrum<T> &);
  };

  template<typename T>
  void PowerSpectrumBoxcar<T>::tavg(double t){
    if(!std::isfinite(t))
      throw std::invalid_argument("Boxcar width must be finite");
    w=std::max(1.0,t);
    size_t new_n=std::round(w);
    if(new_n==n)return;
    if(_h.size()==n)
      std::rotate(_h.begin(),_h.begin()+index,_h.end());
    index=0;
    if(_h.size()>new_n)
      _h.erase(_h.begin(),_h.end()-new_n);
    n=new_n;
    update_average();
  }

  template<typename T>
  PowerSpectrumBoxcar<T> &PowerSpectrumBoxcar<T>::operator()(const PowerSpectrum<T>& p){
    if(this->update(p)){
      _h.clear();
      index=0;
    }
    if(_h.size()<n){
      _h.emplace_back(p.size());
      index=_h.size()-1;
    }
    for(size_t i=0;i<p.size();++i)
      _h[index][i]=p[i];
    index=(index+1)%_h.size();
    update_average();
    _time=p.time();
    return *this;
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

  // template<typename T>
  // bool PowerSpectrumBoxcar<T>::update(const PowerSpectrum<T> &p){
  //   if(p.size()==_delay.size()&&p(0)==_delay[0]
  //      &&p(p.size())==_delay[_delay.size()]&&_name==p.name())
  //     return false;
  //   _delay.resize(p.size());
  //   for(size_t i=0;i<_delay.size();i++)
  //     _delay[i]=p(i);
  //   _power.resize(p.size());
  //   _name=p.name();
  //   return true;
  // }
  
  //
  // PowerSpectrumLeakyBox
  //
  
  template<typename T=float>
  class PowerSpectrumLeakyBox: public PowerSpectrumAverage<T>{
    using PowerSpectrumAverage<T>::_name;
    using PowerSpectrumAverage<T>::_time;
    using PowerSpectrumAverage<T>::_delay;
    using PowerSpectrumAverage<T>::_power;
  public:
    PowerSpectrumLeakyBox(double w=1){tavg(w);}
    virtual ~PowerSpectrumLeakyBox(){};
    double tavg() override {return w;}
    void tavg(double) override;
    PowerSpectrumLeakyBox<T> &operator()(const PowerSpectrum<T>&) override;
  private:
    double w;
    //    bool update(const PowerSpectrum<T> &);
  };

  template<typename T>
  void PowerSpectrumLeakyBox<T>::tavg(double t){
    if(!std::isfinite(t))
      throw std::invalid_argument("Leaky-box width must be finite");
    w=std::max(1.0,t);
  }

  template<typename T>
  PowerSpectrumLeakyBox<T> &PowerSpectrumLeakyBox<T>::operator()(const PowerSpectrum<T>& p){
    if(this->update(p)||w==1){
      for(size_t i=0;i<_power.size();++i)
        _power[i]=p[i];
    }
    else{
      const double alpha=1.0/w;
      for(size_t i=0;i<_power.size();++i)
        _power[i]=alpha*p[i]+(1-alpha)*_power[i];
    }
    _time=p.time();
    return *this;
  }

  // template<typename T>
  // bool PowerSpectrumLeakyBox<T>::update(const PowerSpectrum<T> &p){
  //   if(p.size()==_delay.size()&&p(0)==_delay[0]
  //      &&p(p.size())==_delay[_delay.size()]&&_name==p.name())
  //     return false;
  //   _delay.resize(p.size());
  //   for(size_t i=0;i<_delay.size();i++)
  //     _delay[i]=p(i);
  //   _power.resize(p.size());
  //   _name=p.name();
  //   return true;
  // }

  //
  // CalcPowerSpectrum
  //
  
  template<typename T=float>
  class CalcPowerSpectrum{
  public:
    CalcPowerSpectrum()=default;
    virtual ~CalcPowerSpectrum()=default;
    virtual void operator()(const Phasors<T> &, PowerSpectrum<T> &)=0;
    virtual void operator()(const std::vector<Phasors<T>> &,
			    std::vector<PowerSpectrum<T>> &)=0;
  private:
  };

  
  //
  // CalcPowerSpectrumFFT
  //
  
  void four1(float data[], unsigned long nn, int isign);
  void four1(double data[], unsigned long nn, int isign);

  template<typename T=float>
  class CalcPowerSpectrumFFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumFFT():wfirst(0),wlast(0),nOver(1){};
    CalcPowerSpectrumFFT(const Phasors<T> &, int);
    virtual ~CalcPowerSpectrumFFT(){};
    void set_nover(int n){nOver=static_cast<size_t>(std::max(1,n)); conf();}
    void set_nphasors(int n){nPhasors=static_cast<size_t>(std::max(2,n));
      conf();}
    void set_wavelength(T w0, T w1){wfirst=w0; wlast=21; conf();}
    void operator()(const Phasors<T> &, PowerSpectrum<T> &) override;
    void operator()(const std::vector<Phasors<T>> &,
		    std::vector<PowerSpectrum<T>> &) override;
  private:
    size_t nPhasors;
    T wfirst,wlast;
    size_t nOver;
    size_t nFFT=1;
    std::vector<T> fftData;
    void conf();
  };

  template<typename T>
  CalcPowerSpectrumFFT<T>::CalcPowerSpectrumFFT(const Phasors<T> &p, int n){
    nPhasors=p.size();
    nOver=static_cast<size_t>(n);
    wfirst=p(0);
    wlast=p(p.size()-1);
    conf();
  }
  
  template<typename T>
  void CalcPowerSpectrumFFT<T>::operator()(const Phasors<T> &ph,
					   PowerSpectrum<T> &ps){
    // const size_t n=phasors.size();
    // if(n < 2) throw std::invalid_argument("FFT needs at least two phasors");
    // if(wfirst <= 0 || wlast <= 0)
    //   throw std::logic_error("Call FFT set_wavelength(first,last) before calculation");
    // const double k0=1.0/static_cast<double>(wfirst);
    // const double dk=(1.0/static_cast<double>(wlast)-k0)/(n-1);
    // if(validateWavelengths && validatedSize != n){
    //   // Allow rounding of wavelength coordinates to their declared precision.
    //   const long double eps=std::max(
    //       static_cast<long double>(std::numeric_limits<T>::epsilon()),
    //       static_cast<long double>(std::numeric_limits<U>::epsilon()));
    //   for(size_t i=0;i<n;i++){
    //     if(!std::isfinite(phasors(i)) || phasors(i) <= 0)
    //       throw std::invalid_argument("FFT wavelengths must be positive and finite");
    //     const long double expected=k0+i*dk;
    //     if(std::abs(1.0L/phasors(i)-expected) > 8*eps*std::abs(expected))
    //       throw std::invalid_argument("FFT requires uniform wavenumbers matching its endpoints");
    //   }
    //   validatedSize=n; // Cache only a successful check of the entire grid.
    // }

    // const size_t limit=std::min(fftData.max_size()/2-1,
    //     static_cast<size_t>(std::numeric_limits<unsigned long>::max()/2));
    // if(n > limit/static_cast<size_t>(nOver))
    //   throw std::length_error("FFT size overflow");
    // const size_t required=n*static_cast<size_t>(nOver);
    // size_t nFFT=1; // Number of complex samples, including zero padding.
    // while(nFFT < required){
    //   if(nFFT > limit/2) throw std::length_error("FFT size overflow");
    //   nFFT*=2;
    // }
    // // four1 is one-based: provide a real dummy element, never data()-1.
    // fftData.assign(2*nFFT+1,T(0));
    T dk=(1.0/wlast-1.0/wfirst)/(nPhasors-1);
    if(ph.size()!=nPhasors)
      throw std::runtime_error(std::format("ph.size={} but nPhasors={}",
					   ph.size(),nPhasors));
			       
    for(size_t i=0;i<nPhasors;i++){
      fftData[2*i+1]=static_cast<T>(ph[i].real());
      fftData[2*i+2]=static_cast<T>(ph[i].imag());
    }
    four1(fftData.data(),static_cast<unsigned long>(nFFT),1);
    ps.resize(nFFT);
    ps.name(ph.name());
    ps.time(ph.time());
    const double delayStep=1.0/(nFFT*dk);
    for(size_t i=0;i<nFFT;i++){
      // nFFT is a power of two, so wrapping needs no integer division.
      const size_t bin=(i+nFFT/2)&(nFFT-1);
      const double signedBin=static_cast<double>(i)-nFFT/2;
      const T re=fftData[2*bin+1], im=fftData[2*bin+2];
      ps[i]=static_cast<T>(re*re+im*im);
      ps(i)=static_cast<T>(signedBin*delayStep);
    }
    // Keep delay coordinates increasing for either wavelength ordering.
    if(dk < 0){
      for(size_t i=0;i<nFFT/2;i++){
        std::swap(ps[i],ps[nFFT-1-i]);
        std::swap(ps(i),ps(nFFT-1-i));
      }
    }
  }

  template<typename T>
  void CalcPowerSpectrumFFT<T>::operator()(const std::vector<Phasors<T> > &ph,
					   std::vector<PowerSpectrum<T> > &ps){
    ps.resize(ph.size());
    for(size_t i=0;i<ph.size();i++) operator()(ph[i],ps[i]);
  }

  template<typename T>
  void CalcPowerSpectrumFFT<T>::conf(){
    size_t nrequired=nPhasors*nOver;
    while(nFFT<nrequired)
      nFFT*=2;
    fftData.assign(2*nFFT+1,T(0));
  }
  
  //
  // CalcPowerSpectrumDFT
  //
  
  template<typename T=float>
  class CalcPowerSpectrumDFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumDFT(){};
    CalcPowerSpectrumDFT(const std::vector<T> &d){set_delay(d);}
    ~CalcPowerSpectrumDFT(){};
    template<typename U>
    void set_delay(const std::vector<U> &);
    virtual void operator()(const Phasors<T> &, PowerSpectrum<T> &) override;
    void operator()(const std::vector<Phasors<T>> &,
		    std::vector<PowerSpectrum<T>> &) override;
  private:
    std::vector<T> delay;
  };
    
  template<typename T>
  template<typename U>
  void CalcPowerSpectrumDFT<T>::set_delay(const std::vector<U> &d){
    delay.resize(d.size());
    for(size_t i=0;i<d.size();i++)
      delay[i]=static_cast<T>(d[i]);
  }

  template<typename T>
  void CalcPowerSpectrumDFT<T>::operator()(const Phasors<T> &ph,
					   PowerSpectrum<T> &ps){
    if(delay.empty())
      throw std::logic_error("Configure the DFT delay grid first");
    if(ph.size()==0)
      throw std::invalid_argument("DFT input must not be empty");
    ps.name(ph.name());
    ps.time(ph.time());
    ps.resize(delay.size());
    const double twoPi=2*std::numbers::pi;
    for(size_t i=0;i<delay.size();i++){
      ps(i)=static_cast<T>(delay[i]);
      std::complex<T> sum(0,0);
      std::complex<T> rot;
      for(size_t j=0;j<ph.size();j++){
        rot=std::polar(1.0,twoPi*delay[i]/ph(j));
        sum+=ph[j]*rot;
      }
      ps[i]=static_cast<T>(std::norm(sum));
    }
  }
  
  template<typename T>
  void CalcPowerSpectrumDFT<T>::operator()(const std::vector<Phasors<T> > &ph,
					   std::vector<PowerSpectrum<T> > &ps){
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
    template<typename U>
    void operator()(const PowerSpectrum<U> &);
    T delay(){return _delay;}
    T snr(){return _snr;}
    template<typename U>
    void set_peak_width(U w){peak_width=w;}
  private:
    T peak_width;
    T _snr,_delay;
  };
  
  template<typename T>
  template<typename U>
  void AnalyzePowerSpectrum<T>::operator()(const PowerSpectrum<U> &p){
    int imax=-1;
    T peak=-1e31;
    
    for(size_t i=0;i<p.size();i++)
      if(p[i]>peak){
	imax=i;
	peak=p[i];
      }
    
    _delay=p(imax);
    
    int count=0;
    T v=0,v2=0;
    for(size_t i=0;i<p.size();i++)
      if(p(i)<_delay-peak_width/2||p(i)>_delay+peak_width/2){
	count++;
	v+=p[i];
	v2+=p[i]*p[i];
      }

    T mean=v/count;
    
    _snr=(peak-mean)/sqrt(v2/count-mean*mean);      
  }
  
  //
  // PowerSpectrum smoothing
  //
  template<typename T=float>
  class PowerSpectrumSmooth: public PowerSpectrum<T>{
  public:
    PowerSpectrumSmooth(){};
    virtual ~PowerSpectrumSmooth(){};
    virtual PowerSpectrumSmooth<T> &operator()(const PowerSpectrum<T> &)=0;
    virtual void width(T)=0;
    virtual T width() const =0;
  };

  //
  // PowerSpectrumSmoothNone
  template<typename T=float>
  class PowerSpectrumSmoothNone: public PowerSpectrumSmooth<T>{
  public:
    PowerSpectrumSmoothNone(){};
    ~PowerSpectrumSmoothNone(){};
    virtual PowerSpectrumSmoothNone<T> &operator()(const PowerSpectrum<T> &) override;
    virtual void width(T w) override {_width=w;}
    virtual T width() const override {return _width;}
  private:
    T _width;
  };

  template<typename T>
  PowerSpectrumSmoothNone<T> &PowerSpectrumSmoothNone<T>::operator()(const PowerSpectrum<T> &p){
    for(size_t i=0;i<p.size();i++)
      this->_power[i]=p[i];
    return *this;
  }
  
  //
  // PowerSpectrumGaussian
  //
  template<typename T=float>
  class PowerSpectrumSmoothGaussian: public PowerSpectrumSmooth<T>{
  public:
    PowerSpectrumSmoothGaussian()=default;
    PowerSpectrumSmoothGaussian(T w):_width(w){};
    ~PowerSpectrumSmoothGaussian()=default;
    virtual PowerSpectrumSmoothGaussian<T> &operator()(const PowerSpectrum<T> &) override;
    virtual void width(T w) override {_width=w;}
    virtual T width() const override {return _width;}
  private:
    T _width;
  };
    
  template<typename T>
  PowerSpectrumSmoothGaussian<T> &PowerSpectrumSmoothGaussian<T>::operator()(const PowerSpectrum<T> &p){
    for(size_t i=0;i<p.size();i++){
      T weighted=0,weights=0;
      for(size_t j=0;j<p.size();j++){
	const T x=(p(j)-p(i))/this->_width/2;
	const T weight=std::exp(-x*x/2);
	weighted+=weight*p[j];
	weights+=weight;
      }
      this->_power[i]=weighted/weights;
    }
    return *this;
  }
  
  //
  // PowerSpectrumSBoxcar
  //
  
  template<typename T=float>
  class PowerSpectrumSmoothBoxcar: public PowerSpectrumSmooth<T>{
  public:
    PowerSpectrumSmoothBoxcar()=default;
    PowerSpectrumSmoothBoxcar(T w):_width(w){};
    virtual ~PowerSpectrumSmoothBoxcar()=default;
    virtual PowerSpectrumSmoothBoxcar<T> &operator()(const PowerSpectrum<T> &);
    virtual void width(T w) override {_width=w;}
    virtual T width() const override {return _width;}
  private:
    T _width;
  };

  template<typename T>
  PowerSpectrumSmoothBoxcar<T> &PowerSpectrumSmoothBoxcar<T>::operator()(const PowerSpectrum<T> &p){
    for(size_t i=0;i<p.size();i++){
      T sum=0;
      size_t n=0;
      for(size_t j=0;j<p.size();j++){
	if(std::abs(p(j)-p(i))<this->_width/2){
	  sum+=p[j];
	  n++;
	}
      }
      this->_power[i]=sum/n;
    }
    return *this;
  }
  
  // Gaussian
  // 1/\sqrt{2pi}/sigma*exp(-(x-mu)^2/2\sigma^2)
  // inline double gaussian(double x, double m, double s){
  //   return 1.0/sqrt(2*M_PI)/s*exp(-(x-m)*(x-m)/2.0/s/s);
  // }
  
  // Gaussian smoothing
  // template<typename T=float, typename U=float>
  // PowerSpectrum<T> gsmooth(const PowerSpectrum<T> &p, U sigma){
  //   if(!std::isfinite(sigma) || sigma <= 0)
  //     throw std::invalid_argument("Gaussian width must be positive and finite");
  //   PowerSpectrum<T> result=p;
  //   for(size_t i=0;i<p.size();i++){
  //     long double weighted=0, weights=0;
  //     for(size_t j=0;j<p.size();j++){
  //       const long double x=(static_cast<long double>(p(j))-p(i))/sigma;
  //       const long double weight=std::exp(-x*x/2);
  //       weighted+=weight*p[j];
  //       weights+=weight;
  //     }
  //     result[i]=static_cast<T>(weighted/weights);
  //   }
  //   return result;
  // }

  // Boxcar smoothing
  // template<typename T=float, typename U=float>
  // PowerSpectrum<T> bsmooth(const PowerSpectrum<T> &p, U w){
  //   PowerSpectrum<T> _p=p;
  //   double g,s;
  //   for(size_t i=0;i<p.size();i++){
  //     g=0;
  //     s=0;
  //     for(size_t j=0;j<p.size();j++){
  // 	if(fabs(p(j)-p(i))<w/2){
  // 	  g++;
  // 	  s+=p[j];
  // 	}
  //     }
  //     _p[i]=s/g;
  //   }
  //   return _p;
  // }


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
