#include <iostream>

int FuncLoadTree(double *in, double *out, double *dpar, char *option);

int main()
{
  double in[15];
  double out[15];
  double dpar[15];
  char option[120] = "~/Workspace/FullReco_Mainz/test_nnL_dpi_newgeo_w25cm.root";

  for(int i=0;i<15;++i)
    {
      in[i] = 0.;
      out[i] = 0.;
      dpar[i] = 0.;
    }
  int a = FuncLoadTree(in, out, dpar, option);
  std::cout<<"A "<<a<<std::endl;
  for(int i=0;i<15;++i)
    std::cout<<out[i]<<" ";
  std::cout<<std::endl;

  a = FuncLoadTree(in, out, dpar, option);
  for(int i=0;i<15;++i)
    std::cout<<out[i]<<" ";
  std::cout<<std::endl;

  
  return 0;

}
