#ifndef _AMJINTERFEROMETRY_H_
#define _AMJINTERFEROMETRY_H_

#include <cstring>
#include <complex>
#include <vector>
#include <cmath>
#include <cstdint>

//#include <amjCom/amjPacket.H>
//#include <amjCom/amjPacketRW.H>

#include <iostream>

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
    Phasors(size_t n=0){_p.resize(n);};
    virtual ~Phasors(){};
    const std::string &name(){return _name;};
    void name(std::string name){_name=name;};
    size_t size() const {return _w.size();};
    void resize(size_t n){_w.resize(n); _p.resize(n);};
    T &operator()(int i){return _w[i];};
    const T &operator()(int i) const {return _w[i];}; 
    std::complex<T>& operator[](int i){return _p[i];};
    const std::complex<T>& operator[](int i) const {return _p[i];};
    size_t memsize() const {return memsize1()+memsize2();};
    int write(uint8_t *) const;
    size_t memsize1() const {return 2*sizeof(uint32_t);};
    int read1(uint8_t *);
    size_t memsize2() const {return _name.length()+_p.size()*3*sizeof(T);};
    int read2(uint8_t *);
  private:
    std::string _name;
    std::vector<T> _w;
    std::vector<std::complex<T> > _p;
  };

  template<typename T=float>
  class PowerSpectrum{
  public:
    PowerSpectrum(){}
    virtual ~PowerSpectrum(){};
    void clear(){_delay.resize(0); _power.resize(0);};
    T &operator[](int i){return _power[i];};
    const T &operator[](int i) const {return _power[i];};
    T &operator()(int i){return _delay[i];};
    const T &operator()(int i) const {return _delay[i];}
    const std::vector<T> &delay(){return _delay;}
    const std::vector<T> &power(){return _power;}
    size_t size() const {return _delay.size();};
    void resize(size_t s){_delay.resize(s); _power.resize(s);};
    size_t memsize() const {return memsize1()+memsize2();};
    int write(uint8_t *) const;
    size_t memsize1() const {return sizeof(uint32_t);};
    int read1(uint8_t *);
    size_t memsize2() const {return sizeof(T)*2*_delay.size();}
    int read2(uint8_t *);
    template<typename U=float>
    void operator+=(const PowerSpectrum<U> &ps);
    template<typename U=float>
    void operator/=(U);
  private:
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
    CalcPowerSpectrumFFT():nOver(1){};
    ~CalcPowerSpectrumFFT(){};
    void set_nover(int n){nOver=n;};
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
  };

  template<typename T=float>
  class CalcPowerSpectrumDFT: public CalcPowerSpectrum<T>{
  public:
    CalcPowerSpectrumDFT(){};
    ~CalcPowerSpectrumDFT(){};
    template<typename U>
    void set_delay(const std::vector<U> &);
    template<typename U>
    void set_delay(const std::vector<std::vector<U> > &);
    template<typename U, typename V>
    void operator()(const Phasors<U> &, PowerSpectrum<V> &);
    template<typename U, typename V>
    void operator()(const std::vector<Phasors<U> > &,
			    std::vector<PowerSpectrum<V> > &);
  private:
    std::vector<T> delay;
  };
    
  template<typename T=float>
  class DelayMachine{
  public:
    DelayMachine(){};
    template<typename U>
    DelayMachine(std::string name, const Wavelengths<U> &L):_name(name),L(L){};
    virtual ~DelayMachine()=0;
    std::string name() const {return _name;}
    long unsigned int size() const {return L.size();}
    template<typename U>
    void load(Phasors<U> &p){};
    virtual std::pair<float,float> delay()=0;
  protected:
    std::string _name;
    Wavelengths<T> L;
  };

  template<typename T=float>
  class DelayMachineGD: public DelayMachine<T>{
  public:
    DelayMachineGD(){};
    template<typename U>
    DelayMachineGD(std::string name, const Wavelengths<U> &L, int nOver,
		   int nIncoherent, int nSmooth);
    virtual ~DelayMachineGD(){};
    template<typename U>
    void load(const Phasors<U> &p);
    virtual std::pair<T,T> delay() const {return _delay;};
    PowerSpectrum<T> &p(){return pSpec;};
  private:
    unsigned int nOver;
    unsigned int nIncoherent;
    unsigned int nSmooth;
    
    std::vector<T> fftData;
    PowerSpectrum<T> pSpec;
    std::vector<PowerSpectrum<T> > pSpecHistory;
    unsigned int ipSpecHistory;
    
    PowerSpectrum<T> pSpecX;
    std::vector<T> peakHistory;
    std::vector<T> snrHistory;
    unsigned int iPeakHistory;
    
    //float scale;
    
    std::pair<T,T> _delay;
    void computePSpec();
    void computeAveragePSpec();
    void computeAveragePSpecPeak();
    void computeAveragePSpecPeakSmooth();
  };


  template<typename T>
  int findMax(std::vector<T> &p);
  
  
  
  // Wavelengths

  template<typename T=float>
  Wavelengths<T> makeWavelengths(int nL, T L0, T L1, std::string &name=""){
    Wavelengths w;
    w.name(name);
    w.resize(nL);
    float B0=1/L0,B1=1/L1;
    float dB=(B1-B0)/(nL-1);
    for(unsigned int i=0;i<w.size();i++)
      w[i]=1/(B0+i*dB);
  }

  // Phasors

  template<typename T>
  int Phasors<T>::write(uint8_t *pos) const{
    uint32_t s1=_name.length(),s2=_p.size();
    uint8_t *p=pos;
    p=memcpy(p,&s1,sizeof(uint32_t))+sizeof(uint32_t);
    p=memcpy(p,&s2,sizeof(uint32_t))+sizeof(uint32_t);
    p=memcpy(p,&_w.data(),sizeof(T)*_w.size())+sizeof(T)*_w.size();
    p=memcpy(p,&_p.data(),2*sizeof(T)*_p.size())+2*sizeof(T)*_p.size();
  }
  
  template<typename T>
  int Phasors<T>::read1(uint8_t *pos){
    uint32_t s;
    uint8_t *p=pos;
    p=memcpy(&s,p,sizeof(uint32_t))+sizeof(uint32_t);
    _name.resize(s);
    p=memcpy(&s,p,sizeof(uint32_t))+sizeof(uint32_t);
    _w.resize(s);
    _p.resize(s);    
    return sizeof(uint32_t)+s;
  }
  
  template<typename T>
  int Phasors<T>::read2(uint8_t *pos){
    uint8_t *p=pos;
    p=memcpy(_name.data(),p,_name.length())+_name.length();
    p=memcpy(_w.data(),p,sizeof(T)*_w.size())+sizeof(T)*_w.size();
    p=memcpy(_p.data(),p,2*sizeof(T)*_p.size())+2*sizeof(T)*_p.size();
    return p-pos;
  }

  // PowerSpectrum

  void four1(float data[], unsigned long nn, int isign);


  template<typename T>
  int PowerSpectrum<T>::write(uint8_t *pp) const{
    uint32_t tmp=_delay.size();
    uint8_t *p=(uint8_t *)memcpy(pp,&tmp,sizeof(uint32_t))+sizeof(T);
    p=(uint8_t *)memcpy(p,_delay.data(),sizeof(T)*_delay.size())
      +sizeof(T)*_delay.size();
    return (uint8_t *)memcpy(p,_power.data(),sizeof(T)*_power.size())
      +sizeof(T)*_power.size()-pp;
  }

  template<typename T>
  int PowerSpectrum<T>::read1(uint8_t *pp){
    uint32_t s;
    memcpy(&s,pp,sizeof(uint32_t));
    _delay.resize(s);
    _power.resize(s);
    return sizeof(uint32_t);
  }

  template<typename T>
  int PowerSpectrum<T>::read2(uint8_t *pp){
    uint8_t *p=(uint8_t *)memcpy(_delay.data(),p,sizeof(T)*_delay.size())
      +sizeof(T)*_delay.size();
    return (uint8_t *)memcpy(_power.data(),p,sizeof(T)*_power.size())
      +sizeof(T)*_power.size()-pp;
  }
  
  template<typename T>
  template<typename U>
  void PowerSpectrum<T>::operator+=(const PowerSpectrum<U> &ps){
    // Resize to the longest power spectrum
    int n=ps.size()>size()?ps.size():size();
    _power.resize(n,0);
    for(unsigned int i=0;i<_power.size();i++)
      _power[i]+=ps._power[i];
  }
  
  template<typename T>
  template<typename U>
  void PowerSpectrum<T>::operator/=(U f){
    for(unsigned int i=0;i<_power.size();i++)
      _power[i]/=f;
  }

  // CalcPowerSpectrum
  template<typename T>
  CalcPowerSpectrum<T>::~CalcPowerSpectrum() {}

  
  // CalcPowerSpectrumFFT
  template<typename T>
  void CalcPowerSpectrumFFT<T>::set_wavelength(T _wfirst, T _wlast){
    wfirst=_wfirst;
    wlast=_wlast;
    wavelengths.resize(0);
  }
  
  template<typename T>
  template<typename U>
  void CalcPowerSpectrumFFT<T>::set_wavelength(const
						std::vector<std::pair<U,U> > &w){
    wavelengths=w;
  }
  
  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumFFT<T>::operator()(const Phasors<U> &phasors,
					   PowerSpectrum<V> &pspec){
    // The length of fftData must be the smallest power of 2 which is
    // at least 2*nOver*phasors.size()
    size_t nFFT=1;
    while(nFFT<2*nOver*phasors.size()) nFFT*=2;
    fftData.clear();
    fftData.resize(nFFT,0);
    
    // Load X and Y into array and pad with zeros
    for(unsigned int i=0;i<pspec.size();i++){
      fftData[2*i]=phasors[i].real();
      fftData[2*i+1]=phasors[i].imag();
    }

    // Power spectrum size
    pspec.resize(fftData.size()/2);

    // Compute the delays
    double ndk=(1/wlast-1/wfirst);
    //std::cout << "ndk=" << ndk << std::endl;
    //std::cout << "pspec.size()=" << pspec.size() << std::endl;
    
    // Compute the Fourier transform
    std::vector<float> fftDataNR(fftData.size() + 1);
    std::copy(fftData.begin(), fftData.end(), fftDataNR.begin() + 1);
    
    four1(fftDataNR.data(), fftData.size()/2, 1);

    std::copy(fftDataNR.begin() + 1, fftDataNR.end(), fftData.begin());

    //four1(fftData.data()-1,fftData.size()/2,1);
    for(int i=0;i<(int)pspec.size()/2;i++){
      pspec[i+pspec.size()/2]=fftData[2*i]*fftData[2*i]
	+fftData[2*i+1]*fftData[2*i+1];
      pspec(i+pspec.size()/2)=i/ndk;
      pspec[i]=fftData[2*i+pspec.size()]*fftData[2*i+pspec.size()]+
	fftData[2*i+1+pspec.size()]*fftData[2*i+1+pspec.size()];
      pspec(i)=((double)i-(double)pspec.size()/2)/ndk;
      //pspec(i)=-10;
      //std::cout << i << " " << ndk << " " << pspec.size() << " " << pspec(i) <<  " " << pspec(i+pspec.size()/2) << std::endl;
    }  
  }
  
  
  // CalcPowerSpectrumDFT
  template<typename T>
  template<typename U>
  void CalcPowerSpectrumDFT<T>::set_delay(const std::vector<U> &d){
    delay.resize(d.size());
    for(unsigned int i=0;i<d.size();i++)
      delay[i]=d[i];
  }
  
  template<typename T>
  template<typename U>
  void CalcPowerSpectrumDFT<T>::set_delay(const std::vector<std::vector<U> > &d){
    
  }

  template<typename T>
  template<typename U, typename V>
  void CalcPowerSpectrumDFT<T>::operator()(const Phasors<U> &phasors,
					   PowerSpectrum<V> &pspec){
    pspec.resize(delay.size());
    for(unsigned int i=0;i<delay.size();i++){
      pspec(i)=delay[i];
      std::complex<T> sum=0;
      double angle;
      std::complex<T> rotation;
      for(unsigned int j=0;j<phasors.size();j++){
	angle=2*M_PI*delay[i]/phasors(j);
	rotation=std::polar(1.0,angle);
	sum+=phasors[j]*rotation;
      }
      pspec[i]=std::norm(sum);
    }
  }
  
  
  // DelayMachine

  // DelayMachineGD

  template<typename T>
  template<typename U>
  DelayMachineGD<T>::DelayMachineGD(std::string name,
				    const Wavelengths<U> &L, int nOver,
				    int nIncoherent, int nSmooth):
    DelayMachine<T>(name,L),nOver(nOver),nIncoherent(nIncoherent),nSmooth(nSmooth),
    ipSpecHistory(0),peakHistory(nSmooth,0),snrHistory(nSmooth,0),
    iPeakHistory(0){
    
    // // Compute the size of the FFT to be a power of 2
    // unsigned int nFFT=1;
    // while(nFFT<nOver*L.size()) nFFT*=2;
    // fftData.resize(2*nFFT);
    // PowerSpectrumDelays delays(L,nFFT); // Compute delay values
    pSpec=PowerSpectrum();//L.name(),delays);
    pSpecHistory.clear();
    for(int i=0;i<nIncoherent;i++)
      pSpecHistory.push_back(PowerSpectrum());//L.name(),delays));
  }

  template<typename T>
  void initialize(){ // Reset after any change of parameeters or dimensions
    
    
  }
  
  template<typename T>
  template<typename U>
  void DelayMachineGD<T>::load(const Phasors<U> &phasors){
    // Add power spectrum to history
    pSpecHistory[ipSpecHistory](phasors);
    ipSpecHistory=(ipSpecHistory+1)%pSpecHistory.size();
    
    // Compute the average of the history
    pSpec.clear();
    for(unsigned int i=0;i<nIncoherent;i++)
      pSpec+=pSpecHistory[i];
    pSpec/=nIncoherent;
    
    // Find the peak location and snr
    //computeAveragePSpecPeak();
    computeAveragePSpecPeakSmooth();

    pSpecX=pSpec;
    
    // Find the peak
    unsigned int iMax=findMax(pSpecX.power());
    T max=pSpecX.power()[iMax];
    
    peakHistory[iPeakHistory]=pSpecX.delays()[iMax];
    
    // Zero the 10% largest value
    for(unsigned int i=0;i<pSpecX.size()/10;i++)
      pSpecX.power()[findMax(pSpecX.power())]=0;
    
    // Find the average of the 90% low values
    int nn=0;
    T sum=0;
    for(unsigned int i=0;i<pSpecX.size();i++)
      if(pSpecX.power()[i]>0){
	sum+=pSpecX.power()[i];
	nn++;
      }
    T avg=sum/nn;
    
    // The RMS of the 90% low values around their average is the noise
    T tmp;
    sum=0;
    for(unsigned int i=0;i<pSpecX.size();i++)
      if(pSpecX.power()[i]>0){
	tmp=pSpecX.power()[i]-avg;
	sum+=tmp*tmp;
      }

    // This is the SNR
    snrHistory[iPeakHistory]=(max-avg)/sqrt(sum/nn);
    iPeakHistory=(iPeakHistory+1)%nSmooth;

    // Compute smoothed peak
    _delay={0,0};
    for(unsigned int i=0;i<nSmooth;i++){
      _delay.first+=peakHistory[i];
      _delay.second+=snrHistory[i];
    }
    _delay.first/=nSmooth;
    _delay.second/=nSmooth;
  }
  
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
