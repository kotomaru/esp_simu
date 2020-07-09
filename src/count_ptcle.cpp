//This is a program for counting the Number of protons
//2020/07/08

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
#include<bits/stdc++.h>

using namespace std;

#define NUM_BEAM 1E+3

const double dev_theta = 1.0;//degree
const double min_theta = 30.0;//degree for recoil proton
const double max_theta = 85.0;//degree for recoil proton

ifstream ifs("./data/52Ca_RMF_295.OBS");
ofstream ofs("./data/52Ca_RMF_295_lab_calc.OBS",ios_base::trunc);

double dCmToLab(double xcm){//degree CM to LAB
  double xlab = 0.0008*pow(xcm,2.0)-0.6664*xcm+90.063;
  return xlab;
}
double cCmToLab(double xcm,double ccm){//cross section CM to LAB
  double clab = 2.*sqrt(2.)*sqrt(1.+cos(xcm*M_PI/180.))*ccm;
  return clab;
}


/*------------------------------------------------*/
int main(int argc, char *argv[])
{
  if(ifs.fail()){
    cerr<<"Failed to open input file"<<endl;
    return -1;
  }
  if(ofs.fail()){
    cerr<<"Failed to open output file"<<endl;
    return -1;
  }
  vector<double> deg_cm;
  vector<double> crs_cm;
  vector<double> deg_lab;
  vector<double> crs_lab;
  double tmp1,tmp2,tmp3,tmp4;
  while(ifs >> tmp1 >> tmp2 >> tmp3 >> tmp4){
    deg_cm.push_back(tmp1);
    crs_cm.push_back(tmp2);
    deg_lab.push_back(dCmToLab(tmp1));
    crs_lab.push_back(cCmToLab(tmp1,tmp2));

    //    ofs << deg_lab[i]<<" "<<crssec_lab[i]<<endl;
    //  cout<<i<<":"<<deg_lab[i]<<":"<<crssec_lab[i]<<endl;
  }

  vector<double> theta_rcl;//recoil degree (about 30~85)
  vector<double> crs_rcl;//calculated cross section from the nearest values

  int i=0; double tmp_deg=0;
  double n_tmp_deg; //search near value of tmp_deg
  
  while(tmp_deg<max_theta){
    tmp_deg = min_theta + dev_theta*i;\
    double tmp_diff = deg_lab.front() - tmp_deg;//under
    for(int j=1;j<deg_lab.size();j++){
      double x = deg_lab.at(j) - tmp_deg;
      if(x>=0 && x < tmp_diff){
	tmp_diff = x;
	n_tmp_deg = j;
      }else{break;}
    }
    if(n_tmp_deg==0){cerr<<"ERROR @ theta search"<<endl;return -1;}
    
    theta_rcl.push_back(tmp_deg);
    double x1,x2,y;
    x1 = tmp_deg - deg_lab.at(n_tmp_deg+1);
    x2 = deg_lab.at(n_tmp_deg) - tmp_deg;
    y = (x2*crs_lab.at(n_tmp_deg+1) + x1*crs_lab.at(n_tmp_deg))
      /(x1+x2);
    crs_rcl.push_back(y);
    ofs << theta_rcl[i]<<" "<<crs_rcl[i]<<endl;
    i++;
  }
  cout<<theta_rcl[40]<<":"<<crs_rcl[40]<<endl;
  ifs.close();// ofs.close();
  
  return 0;

}
