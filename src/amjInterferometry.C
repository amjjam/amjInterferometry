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

bool DelayMachineConfig::operator==(const DelayMachineConfig& c) const {
  return time_method==c.time_method && time_width==c.time_width &&
    delay_axis_method==c.delay_axis_method && delay_axis_width==c.delay_axis_width &&
    peak_exclusion_width==c.peak_exclusion_width &&
    output_method==c.output_method && output_width==c.output_width;
}
void DelayMachineConfig::validate() const {
  using T=TemporalMethod; using A=AxisMethod;
  auto valid=[](T m){return m==T::None || m==T::Boxcar || m==T::Leaky;};
  if(!valid(time_method)||!valid(output_method)||
     (delay_axis_method!=A::None && delay_axis_method!=A::Boxcar && delay_axis_method!=A::Gaussian)||
     !time_width || !output_width || !std::isfinite(delay_axis_width)||
     delay_axis_width<0 || (delay_axis_method!=A::None && delay_axis_width==0)||
     !std::isfinite(peak_exclusion_width)||peak_exclusion_width<0)
    throw std::invalid_argument("Invalid DelayMachine configuration");
}
SmoothPowerSpectrum::SmoothPowerSpectrum(Method m,unsigned int w){configure(m,w);}
void SmoothPowerSpectrum::configure(Method m,unsigned int w){
  DelayMachineConfig c;c.time_method=m;c.time_width=w;c.validate();
  method_=m;width_=w;decay_=std::exp(-1.0f/w);
  history_.resize(m==Method::Boxcar?w:0);reset();
}
void SmoothPowerSpectrum::reset(){
  clear();time(amjTime());next_=count_=0;
  for(auto& p:history_)p.clear();
}
void SmoothPowerSpectrum::add(const PowerSpectrum<float>& p){
  if(!p.size())throw std::invalid_argument("Empty smoothing input");
  if(size() && size()!=p.size())throw std::invalid_argument("Smoothing size changed; reset first");
  if(method_==Method::None){PowerSpectrum<float>::operator=(p);return;}
  if(method_==Method::Leaky){
    if(!size())PowerSpectrum<float>::operator=(p);
    else {
      for(size_t i=0;i<size();++i)(*this)[i]=decay_*(*this)[i]+(1-decay_)*p[i];
      name(p.name());time(p.time());
    }
    return;
  }
  history_[next_]=p;
  PowerSpectrum<float>::operator=(p);
  count_=std::min(count_+1,static_cast<size_t>(width_));
  for(size_t j=0;j<count_;++j)if(j!=next_)
    for(size_t i=0;i<size();++i)(*this)[i]+=history_[j][i];
  (*this)/=count_;
  next_=(next_+1)%width_;
}
DelayMachine::DelayMachine(std::string name,Calculator calc,const Config& c)
 : _name(std::move(name)),calculateSpectrum(std::move(calc)) {
  if(!calculateSpectrum)throw std::invalid_argument("Empty spectrum calculator");
  c.validate();config_=c;
  temporal_.configure(c.time_method,c.time_width);
  peakHistory.resize(c.output_width);snrHistory.resize(c.output_width);
  outputDecay_=std::exp(-1.0f/c.output_width);
  use_default_analyzer();
}
DelayMachine::DelayMachine(std::string name,Calculator calc,int ni,int ns,Analyzer analyzer)
 : DelayMachine(std::move(name),std::move(calc),[&]{
     if(ni<1||ns<1)throw std::invalid_argument("History lengths must be positive");
     Config c;c.time_width=ni;c.output_method=Config::TemporalMethod::Boxcar;c.output_width=ns;return c;
   }()) {
  if(!analyzer)throw std::invalid_argument("Empty spectrum analyzer");
  // Recognize the default callable so peak width remains configurable.
  customAnalyzer_=analyzer.target<AnalyzePowerSpectrum<float>>()==nullptr;
  analyzeSpectrum=std::move(analyzer);legacyOutputStartup_=true;
}
void DelayMachine::set_config(const Config& c){
  c.validate();if(c==config_)return;
  const bool timeChange=c.time_method!=config_.time_method||c.time_width!=config_.time_width;
  const bool axisChange=c.delay_axis_method!=config_.delay_axis_method||c.delay_axis_width!=config_.delay_axis_width;
  const bool analysisChange=c.peak_exclusion_width!=config_.peak_exclusion_width;
  if(timeChange)temporal_.configure(c.time_method,c.time_width);
  if(timeChange||axisChange){axisOffsets_.clear();pSpec.clear();}
  // Any changed upstream weighting invalidates the final output history.
  peakHistory.assign(c.output_width,0);snrHistory.assign(c.output_width,0);
  nextOutput_=outputCount_=0;_delay=_snr=0;
  config_=c;legacyOutputStartup_=false;outputDecay_=std::exp(-1.0f/c.output_width);
  if(analysisChange&&!customAnalyzer_){
    AnalyzePowerSpectrum<float> a;a.set_peak_width(c.peak_exclusion_width);analyzeSpectrum=a;
  }
}
void DelayMachine::set_calculator(Calculator c){
  if(!c)throw std::invalid_argument("Empty spectrum calculator");
  calculateSpectrum=std::move(c);reset();
}
void DelayMachine::set_analyzer(Analyzer a){
  if(!a)throw std::invalid_argument("Empty spectrum analyzer");
  analyzeSpectrum=std::move(a);customAnalyzer_=true;reset();
}
void DelayMachine::use_default_analyzer(){
  AnalyzePowerSpectrum<float> a;a.set_peak_width(config_.peak_exclusion_width);
  analyzeSpectrum=a;customAnalyzer_=false;reset();
}
void DelayMachine::reset(){
  raw_.clear();temporal_.reset();pSpec.clear();pSpec.time(amjTime());
  axisOffsets_.clear();axisIndices_.clear();axisWeights_.clear();
  std::fill(peakHistory.begin(),peakHistory.end(),0);
  std::fill(snrHistory.begin(),snrHistory.end(),0);
  nextOutput_=outputCount_=0;_delay=_snr=0;
}
void DelayMachine::build_axis_weights(){
  axisOffsets_.clear();axisIndices_.clear();axisWeights_.clear();
  const double width=config_.delay_axis_width;
  for(size_t i=0;i<temporal_.size();++i){
    axisOffsets_.push_back(axisIndices_.size());double sum=0;
    for(size_t j=0;j<temporal_.size();++j){
      double x=(double(temporal_(j))-temporal_(i))/width;
      double w=config_.delay_axis_method==Config::AxisMethod::Boxcar
        ? (std::abs(x)<0.5?1:0) : std::exp(-4*x*x);
      if(w>0){axisIndices_.push_back(j);axisWeights_.push_back(w);sum+=w;}
    }
    for(size_t k=axisOffsets_.back();k<axisIndices_.size();++k)axisWeights_[k]/=sum;
  }
  axisOffsets_.push_back(axisIndices_.size());
}
void DelayMachine::load(const Phasors& ph){
  calculateSpectrum(ph,raw_);
  temporal_.add(raw_);pSpec=temporal_;pSpec.time(ph.time());
  if(config_.delay_axis_method!=Config::AxisMethod::None){
    if(axisOffsets_.empty())build_axis_weights();
    for(size_t i=0;i<pSpec.size();++i){
      double sum=0;
      for(size_t k=axisOffsets_[i];k<axisOffsets_[i+1];++k)
        sum+=axisWeights_[k]*temporal_[axisIndices_[k]];
      pSpec[i]=float(sum);
    }
  }
  float snr,delay;analyzeSpectrum(pSpec,snr,delay);
  using M=Config::TemporalMethod;
  if(config_.output_method==M::None){_delay=delay;_snr=snr;return;}
  if(config_.output_method==M::Leaky){
    if(!outputCount_){_delay=delay;_snr=snr;outputCount_=1;}
    else {_delay=outputDecay_*_delay+(1-outputDecay_)*delay;
          _snr=outputDecay_*_snr+(1-outputDecay_)*snr;}
    return;
  }
  peakHistory[nextOutput_]=delay;snrHistory[nextOutput_]=snr;
  nextOutput_=(nextOutput_+1)%peakHistory.size();
  outputCount_=std::min(outputCount_+1,peakHistory.size());
  _delay=_snr=0;
  for(size_t i=0;i<outputCount_;++i){_delay+=peakHistory[i];_snr+=snrHistory[i];}
  const size_t divisor=legacyOutputStartup_?peakHistory.size():outputCount_;
  _delay/=divisor;_snr/=divisor;
}

}
