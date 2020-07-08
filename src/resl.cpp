//This is a program for calculating Eex resolution
//2020/07/06

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <string.h>

using namespace std;

#define OP_C 0
/* -c -> case 1,2 or 3 */
#define NUM_OP 1

struct OptTable{
  char *option;
  int n;
}OptTable[]={
  {"-c",0}  
};

int CheckOption(char *arg){
  int jj=0;
  for(;jj<NUM_OP;jj++){
    if(strcmp(arg,OptTable[jj].option)==0) return(OptTable[jj].n);
  }
  return -1;
}

/*------------------------------------------------*/
int main(int argc, char *argv[])
{
  double Eex=60., Mbeam=52000.; 
  int num_case=1;//default
  
  if((argc==2*NUM_OP+1 || argc == 1)!=1){
    cout<<"Usage : ./simu_resl -c [case opt]"<<endl;
    cout<<"default usage is ./simul_resl"<<endl;
  } 
  
  for(int i=1;i<argc;i++){
    if(*argv[i] == '-'){
      if(CheckOption(argv[i])==OP_C){
	num_case = atoi(argv[i+1]);
      }
    }
  }

  cout<< "Energy resolution [case "<<num_case<<"] "
      <<": INPUT each RESOLUTION"<<endl;

  double dEex;
  double dTheta_p, dKE_p;//protn resolution
  double dTheta_b, dKE_b;//beam resolution
  
  /* for case 1 (conventional ESPRI) */
  const double ll = 1.0; //meter 
  double dx_b, dx_p, dz, dtheta_mlt;
  /* for case 2,3 (2 sillicon or +MWDC) */
  const double ll1[2] = {0.15,0.5}; //meter{for case2, case3} 
  const double ll2[2] = {0.15,0.5}; //meter 
  double dx_s1, dx_s2, dtheta_mlt1, dtheta_mlt2;

  if(num_case==1){
    cout<<"L = "<<ll<<" m"<<endl;
  }else if(num_case==2){
    cout<<"L1,L2 = "<<ll1[num_case-2]<<","<<ll2[num_case-2]<<" m"<<endl;
  }else if(num_case==3){
    cout<<"L1,L2 = "<<ll1[num_case-2]<<","<<ll2[num_case-2]<<" m"<<endl;
  }else{
    cout<<"ERROR"<<endl;
    return 1;
  }


  
  cout<<"dKE_p (MeV) : "; cin>>dKE_p;
  cout<<"dKE_b (MeV) : "; cin>>dKE_b;
  cout<<"dTheta_b (mrad): "; cin>>dTheta_b;
  cout<<"-------------------------"<<endl;
  
  switch(num_case){
  case 1:
    cout<<"dx_b (mm): "; cin>>dx_b;
    cout<<"dx_p (mm): "; cin>>dx_p;
    cout<<"dz (mm): "; cin>>dz;
    cout<<"dtheta_mlt (mrad) : "; cin>>dtheta_mlt;
    dTheta_p = pow(dTheta_b,2.0)
      + pow(dx_b/ll,2.0)
      + pow(dx_p/ll,2.0)
      + pow(dz/ll,2.0)
      + pow(dtheta_mlt,2.0);
    dTheta_p = sqrt(dTheta_p);
    break;
    
  case 2:
  case 3:
    cout<<"dx_s1 (mm): "; cin>>dx_s1;
    cout<<"dx_s2 (mm): "; cin>>dx_s2;
    cout<<"dtheta_mlt1 (mrad) @SHT : "; cin>>dtheta_mlt1;
    cout<<"dtheta_mlt2 (mrad) @1st-det: "; cin>>dtheta_mlt2;
    
    dTheta_p = pow(dTheta_b,2.0)
      + pow(dx_s1/ll1[num_case-2],2.0)
      + pow(dx_s2/ll2[num_case-2],2.0)
      + pow(dtheta_mlt1,2.0)
      + pow(dtheta_mlt2,2.0);
    dTheta_p = sqrt(dTheta_p);
    break;
  }
    

  dEex = (pow(dTheta_p,2.)+pow(dKE_p,2.0)+pow(dKE_b,2.0))
  	  /pow(Eex+Mbeam,2.0);
  dEex = sqrt(dEex);

  
  cout<<"\n-------------------\n"<<endl;
  cout<<"Eres = "<<dEex<<endl;
  //   cout<<"d-Theta_p (mrad) = "<<dTheta_p<<endl;

  return 0;

}
