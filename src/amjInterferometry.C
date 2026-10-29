#include "../include/amjInterferometry.H"

#include <iostream>

namespace amjInterferometry{
  // Begin stolen from Numerical Recipes
#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
  
  template<typename T>
  void four1_impl(T data[], unsigned long nn, int isign)
  {
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    T tempr,tempi;
    if(!data || nn==0 || (nn & (nn-1)) ||
       nn > std::numeric_limits<unsigned long>::max()/2 ||
       (isign != 1 && isign != -1))
      throw std::invalid_argument("Invalid four1 arguments");
    
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
      if (j > i) {
	SWAP(data[j],data[i]);
	SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m) {
	j -= m;
	m >>= 1;
      }
      j += m;
    }
    mmax=2;
    while (n > mmax) {
      istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2) {
	for (i=m;i<=n;i+=istep) {
	  j=i+mmax;
	  tempr=wr*data[j]-wi*data[j+1];
	  tempi=wr*data[j+1]+wi*data[j];
	  data[j]=data[i]-tempr;
	  data[j+1]=data[i+1]-tempi;
	  data[i] += tempr;
	  data[i+1] += tempi;
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
    }
  }
  void four1(float data[], unsigned long nn, int isign){
    four1_impl(data,nn,isign);
  }
  void four1(double data[], unsigned long nn, int isign){
    four1_impl(data,nn,isign);
  }
#undef SWAP
  // End stolen from Numerical Recipes

DelayMachine::DelayMachine(std::string name, Calculator calculator,
                           int nIncoherent, int nSmooth, Analyzer analyzer)
    : _name(std::move(name)), calculateSpectrum(std::move(calculator)),
      analyzeSpectrum(std::move(analyzer)) {
  if(nIncoherent<1 || nSmooth<1)
    throw std::invalid_argument("DelayMachine history lengths must be positive");
  if(!calculateSpectrum || !analyzeSpectrum)
    throw std::invalid_argument("DelayMachine requires a calculator and analyzer");
  pSpecHistory.resize(static_cast<size_t>(nIncoherent));
  peakHistory.resize(static_cast<size_t>(nSmooth),0);
  snrHistory.resize(static_cast<size_t>(nSmooth),0);
}

void DelayMachine::set_calculator(Calculator calculator) {
  if(!calculator) throw std::invalid_argument("Empty spectrum calculator");
  calculateSpectrum=std::move(calculator);
  reset();
}

void DelayMachine::set_analyzer(Analyzer analyzer) {
  if(!analyzer) throw std::invalid_argument("Empty spectrum analyzer");
  analyzeSpectrum=std::move(analyzer);
  reset();
}

void DelayMachine::reset() {
  for(auto& spectrum : pSpecHistory) spectrum.clear();
  pSpec.clear();
  pSpec.time(amjTime());
  std::fill(peakHistory.begin(),peakHistory.end(),0);
  std::fill(snrHistory.begin(),snrHistory.end(),0);
  ipSpecHistory=0;
  iPeakHistory=0;
  _delay=0;
  _snr=0;
}

void DelayMachine::load(const Phasors& phasors) {
  calculateSpectrum(phasors,pSpecHistory[ipSpecHistory]);
  if(pSpecHistory[ipSpecHistory].size()==0)
    throw std::invalid_argument("Calculator returned an empty spectrum");

  // Preserve the existing summation order, including startup with partial history.
  pSpec=pSpecHistory[ipSpecHistory];
  size_t n=1; // The newest spectrum is already included.
  for(size_t i=0;i<pSpecHistory.size();i++) {
    if(i==ipSpecHistory || pSpecHistory[i].size()!=pSpec.size()) continue;
    pSpec+=pSpecHistory[i];
    ++n;
  }
  ipSpecHistory=(ipSpecHistory+1)%pSpecHistory.size();
  pSpec/=n;
  // Adding history entries copies their timestamps; report the newest input.
  pSpec.time(phasors.time());

  analyzeSpectrum(pSpec,snrHistory[iPeakHistory],peakHistory[iPeakHistory]);
  iPeakHistory=(iPeakHistory+1)%peakHistory.size();

  // Preserve existing smoothing: unused startup entries are zero, and the
  // divisor is always the configured smoothing length.
  _delay=0;
  _snr=0;
  for(size_t i=0;i<peakHistory.size();i++) {
    _delay+=peakHistory[i];
    _snr+=snrHistory[i];
  }
  _delay/=peakHistory.size();
  _snr/=peakHistory.size();
}

}
