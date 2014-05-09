#ifndef TH2INFO_H
#define TH2INFO_H

#include <map>
#include <string>
#include <TFile.h>

using namespace std;

enum TH2List {
  ABCDana_2D_, 
  ABCDval_2D_, 
  TH2_Size_ 
};

class TH2Info {
  public: 
    TH2Info (
        bool        output, 
        std::string name, 
        std::string title, 
        std::string xtitle, 
        std::string ytitle, 
        std::string xunit, 
        std::string yunit, 
        int 	      xbin, 
        double      xmin, 
        double      xmax, 
        int 	      ybin, 
        double      ymin, 
        double      ymax  
        ) : 
      _output (output), 
      _name   (name), 
      _title  (title), 
      _xtitle (xtitle), 
      _ytitle (ytitle), 
      _xunit  (xunit), 
      _yunit  (yunit), 
      _xbin   (xbin), 
      _xmin   (xmin), 
      _xmax   (xmax), 
      _ybin   (ybin), 
      _ymin   (ymin), 
      _ymax   (ymax) { } 
    bool        _output;
    std::string _name;
    std::string _title;
    std::string _xtitle;
    std::string _ytitle;
    std::string _xunit;
    std::string _yunit;
    int 	      _xbin;
    double      _xmin;
    double      _xmax;
    int 	      _ybin;
    double      _ymin;
    double      _ymax;
};

#endif  

#ifdef TH2INFO_H
#ifndef TH2INFOCLASS_H
#define TH2INFOCLASS_H

template<typename TH2_type> 
class TH2InfoClass{
  public:
    TH2InfoClass() {};
    void CreateTH2();
    void CreateTH2(edm::Service<TFileService> f);
    void CreateTH2( TFile* f, std::string dir_name ); 
    void SetTitles(); 
    void Sumw2();
    TH2_type* GetTH2(std::string _name_);
    TH2Info GetVar(std::string _name_);
    TH2Info GetVar(int index);

  private:
    map<std::string, TH2_type*> mapTH2;
    map<std::string, int> indexTH2;
    static TH2Info Var[TH2_Size_];
};

template<typename TH2_type> 
TH2Info TH2InfoClass<TH2_type>::Var[TH2_Size_] = {
  TH2Info( 1, "ABCDana_2D",			"",	"",	"Pruned Mass(Higgs)", 		"", "GeV/c^{2}",	2, 0, 2, 151, 0, 151), 
  TH2Info( 1, "ABCDval_2D",			"",	"",	"Pruned Mass(Higgs)", 		"", "GeV/c^{2}",	2, 0, 2, 151, 0, 151)  
} ; 

template<typename TH2_type> 
void TH2InfoClass<TH2_type>::CreateTH2(){
  for(int i=0; i<TH2_Size_; i++){ 
    mapTH2[Var[i]._name] = new TH2(Var[i]._name.c_str(),"",Var[i]._xbin, Var[i]._xmin, Var[i]._xmax, Var[i]._ybin, Var[i]._ymin, Var[i]._ymax);
  }

}
template<typename TH2_type> 
void TH2InfoClass<TH2_type>::CreateTH2(edm::Service<TFileService> f){
  for(int i=0; i<TH2_Size_; i++){ 
    mapTH2[Var[i]._name] = f->make<TH2_type>(Var[i]._name.c_str(),"",Var[i]._xbin, Var[i]._xmin, Var[i]._xmax, Var[i]._ybin, Var[i]._ymin, Var[i]._ymax);
  }
}
template<typename TH2_type> 
void TH2InfoClass<TH2_type>::CreateTH2( TFile* f, std::string dir_name="" ){
  for(int i=0; i<TH2_Size_; i++){ 
    mapTH2[Var[i]._name] = (TH2_type*)f->Get( (dir_name+Var[i]._name).c_str());
  }
}

template<typename TH2_type> 
void TH2InfoClass<TH2_type>::SetTitles(){
  for(int i=0; i<TH2_Size_; i++){ 
    mapTH2[Var[i]._name]->SetXTitle( Var[i]._xtitle.c_str());
    mapTH2[Var[i]._name]->SetYTitle( Var[i]._ytitle.c_str());
  }
}

template<typename TH2_type> 
void TH2InfoClass<TH2_type>::Sumw2(){
  for(int i=0; i<TH2_Size_; i++){ 
    mapTH2.find(Var[i]._name)->second->Sumw2();
  }
}

template<typename TH2_type> 
TH2_type* TH2InfoClass<TH2_type>::GetTH2(std::string _name_){
  return mapTH2.find(_name_)->second;
}

template<typename TH2_type> 
TH2Info TH2InfoClass<TH2_type>::GetVar(std::string _name_){
  return Var[indexTH2.find(_name_)->second];
}
template<typename TH2_type> 
TH2Info TH2InfoClass<TH2_type>::GetVar(int index){
  return Var[index];
}

#endif
#endif
