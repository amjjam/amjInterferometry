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

#include <amjTime.H>

namespace amjInterferometry{
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
    virtual ~Wavelengths(){};
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
  private:
    std::string _name;
    amjTime _time;
    std::vector<T> _w;
    std::vector<std::complex<T> > _p;
  };

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

  template<typename T=float>
  class CalcPowerSpectrumFFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumFFT():wfirst(0),wlast(0),nOver(1),
      validateWavelengths(false),validatedSize(0){};
    ~CalcPowerSpectrumFFT(){};
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
  int findMax(std::vector<T> &p);
  
  
  
  // Wavelengths

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

  // Phasors

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

  // PowerSpectrum

  // Legacy one-based interleaved FFT API: data[0] is unused.
  // Allocate 2*nn+1 elements; nn is a positive power of two.
  void four1(float data[], unsigned long nn, int isign);
  void four1(double data[], unsigned long nn, int isign);


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

  // Gaussian
  // 1/\sqrt{2pi}/sigma*exp(-(x-mu)^2/2\sigma^2)
  inline double gaussian(double x, double m, double s){
    return 1.0/sqrt(2*M_PI)/s*exp(-(x-m)*(x-m)/2.0/s/s);
  }
  
  // SmoothPowerSpectrum
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

  template<typename T, typename U>
  PowerSpectrum<T> bsmooth(const PowerSpectrum<T> &p, U w){
    PowerSpectrum<T> _p=p;
    double t,g,s;
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

#if 0 // Unfinished temporal smoothing classes; not yet a supported API.
  template<typename T>
  class SmoothPowerSpectrum<T>: public PowerSpectrum<T>{
  public:
    SmoothPowerSpectrum(){};
    virtual ~SmoothPowerSpectrum(){};
    template<typename U>
    SmoothPowerSpectrum<T> &add(PowerSpectrum<U> &)=0;
  private:
    
  };

  template<typename T>
  class BoxcarSmoothPowerSpectrum<T>: public SmoothPowerSpectrum<T>{
  public:
    BoxcarSmoothPowerSpectrum(size_t n):n(n);
    virtual ~BoxcarSmoothPowerSpectrum(){};
    template<typename U>
    BoxcarSmoothPowerSpectrum<T> &add(PowerSpectrum<U> &);
  private:
    size_t n;
  };
  
  template<typename T>
  class LeakyboxSmoothPowerSpectrum<T>: public SmoothPowerSpectrum<T>{
  public:
    SmoothLeakyPowerSpectrum(size_t n):n(n){};
    virtual SmoothLeakyPowerSpectrum(){};
    template<typename U>
    LeakyboxSmoothPowerSpectrum<T> &add(PowerSpectrum<U> &);
  private:
    size_t n;
  };

  template<typename T>
  template<typename U>
  LeakyboxSmoothPowerSpecrum<T> &
  LeakyboxSmoothPowerSpectrum<T>::add(PowerSpectrum<U> &){

    return *this;
  }
    
#endif

  // CalcPowerSpectrum

  template<typename T>
  CalcPowerSpectrum<T>::~CalcPowerSpectrum() {}

  
  // CalcPowerSpectrumFFT

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

  // CalcPowerSpectrumDFT

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
  


// One concrete phasors -> power spectrum -> delay/SNR pipeline.
// Configuration changes must not run concurrently with load().
class DelayMachine {
public:
  using Phasors = amjInterferometry::Phasors<float>;
  using PowerSpectrum = amjInterferometry::PowerSpectrum<float>;
  using Calculator = std::function<void(const Phasors&, PowerSpectrum&)>;
  // The analyzer argument order matches AnalyzePowerSpectrum: SNR, then delay.
  using Analyzer = std::function<void(const PowerSpectrum&, float&, float&)>;

  DelayMachine(std::string name, Calculator calculator,
               int nIncoherent, int nSmooth,
               Analyzer analyzer=amjInterferometry::AnalyzePowerSpectrum<float>());

  const std::string& name() const { return _name; }
  void load(const Phasors& phasors);
  float delay() const { return _delay; }
  float snr() const { return _snr; }
  const amjTime& time() const { return pSpec.time(); }
  const PowerSpectrum& p() const { return pSpec; }

  // Both setters reset history. Configure callable objects before installation.
  void set_calculator(Calculator calculator);
  void set_analyzer(Analyzer analyzer);
  void reset();

private:
  std::string _name;
  Calculator calculateSpectrum;
  Analyzer analyzeSpectrum;
  PowerSpectrum pSpec;
  std::vector<PowerSpectrum> pSpecHistory;
  size_t ipSpecHistory=0;
  std::vector<float> peakHistory, snrHistory;
  size_t iPeakHistory=0;
  float _delay=0, _snr=0;
};

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
  
  template<typename T>
  template<typename U, typename V>
  void AnalyzePowerSpectrum<T>::operator()(const PowerSpectrum<U> &p, V &snr,
					   V &delay){
    int imax=-1;
    T peak=-1e31;
    
    for(int i=0;i<p.size();i++)
      if(p[i]>peak){
	imax=i;
	peak=p[i];
      }
    
    delay=p(imax);
    
    int count=0;
    T v=0,v2=0;
    for(int i=0;i<p.size();i++)
      if(p(i)<delay-peak_width/2||p(i)>delay+peak_width/2){
	count++;
	v+=p[i];
	v2+=p[i]*p[i];
      }

    T mean=v/count;
    
    snr=(peak-mean)/sqrt(v2/count-mean*mean);      
  }
}

#endif /* _AMJINTERFEROMETRY */
