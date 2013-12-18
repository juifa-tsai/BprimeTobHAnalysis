//////////////////////////////////////////////////////////////
// Usage: 
// TH1D* hnum   = "Some histogram initialization"
// TH1D* hdenom = "Some histogram initialization"
// TGraphAsymmErrors *grEff = DefineGraph(hdenom, hnum) 
/////////////////////////////////////////////////////////////"w

double error1 (double n1 , double n2)
{
  if (n1 == 0)
  {
    return (0);
  }
  else
  {
    return ( (n2/n1 + 0.5/n1 + sqrt(n2/pow(n1,2)*(1-n2/n1) + 0.25/pow(n1,2)))/(1+1.0/n1) - n2/n1 ) ;
  }
}


double error2 (double n1 , double n2)
{
  if (n1 == 0)
  {
    return (0);
  }
  else
  {
    return ( n2/n1 - (n2/n1 + 0.5/n1 - sqrt(n2/pow(n1,2)*(1-n2/n1) + 0.25/pow(n1,2)))/(1+1.0/n1) ) ;
  }
}


//----- User Given Numbers To Draw The TurnOn ----------------//
const int ntrig = 4, neta1 = 2, neta2 = 4, npu1 = 0, npu2 = 2 ;


TGraphAsymmErrors *DefineGraph(TH1D *h1, TH1D *h2) { 

  TGraphAsymmErrors *gr;
  int n1 = h1->GetNbinsX(); int n2 = h2->GetNbinsX();
  if(n1 != n2)
    cout<<"Warning !!!! The two Input Histograms have Different Binning. "<<endl;

  TH1D *h3; h3 = (TH1D*)h2->Clone();
  h3->Divide(h1);

  double *ratio, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  ratio=new double[n1]; errorU=new double[n1]; errorD=new double[n1];
  XerrorU=new double[n1+1]; XerrorD=new double[n1+1]; Xval=new double[n1+1]; 

  for (int i=0; i<n1; i++)
  {
    ratio[i] = h3->GetBinContent(i+1);
    errorU[i] = error1(h1->GetBinContent(i+1), h2->GetBinContent(i+1));
    errorD[i] = error2(h1->GetBinContent(i+1), h2->GetBinContent(i+1));
    Xval[i] = h1->GetBinCenter(i+1);
  } // for (int i=0; i<n1; i++)

  XerrorD[0] = 0; XerrorU[n1] = 0;

  for (int i=0; i<n1; i++)
  {
    XerrorD[i] =  h1->GetBinWidth(i+1)/2;
    XerrorU[i] =  h1->GetBinWidth(i+1)/2;         
  }

  gr = new TGraphAsymmErrors(n1, Xval, ratio, XerrorD, XerrorU, errorD, errorU);

  return gr;

  delete [] ratio; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
  ratio = NULL; errorU = NULL; errorD = NULL; XerrorU = NULL; XerrorD = NULL; Xval = NULL;

} // end of function

