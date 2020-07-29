//This is a program for counting the Number of protons
//2020/07/08

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <bits/stdc++.h>

using namespace std;

#define NUM_BEAM 259200E+4 // 3days 1E+3 cps beam
#define NUM_AVO 6.02214086E+23

const double thick_sht = 0.1;//cm
const double dens_sht = 0.0868;// g/cm^3
const double eff_det = 0.7;

const double dev_theta = 1.0;//degree
const double min_theta = 30.0;//degree for recoil proton
const double max_theta = 85.0;//degree for recoil proton

const double dev_phi = 0.05;//degree
const double min_phi = -90.0;//degree for recoil proton
const double max_phi = 90.0;//degree for recoil proton

/*RDC position in lab frame (rotated)*/
const double RDCx = 100.0; //cm
const double RDCymin = -16.0; 
const double RDCymax = 16.0;
const double RDCzmin = 5.0; 
const double RDCzmax = 38.0; 

ifstream ifs("./data/52Ca_RMF_295.OBS");
//ofstream ofs("./data/52Ca_RMF_295_lab_calc.OBS",ios_base::trunc);
ofstream ofs("./result/N_recoil.dat",ios_base::trunc);

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
  vector<double> theta_cm;//dgree
  vector<double> crs_cm;
  vector<double> theta_lab;
  vector<double> crs_lab;
  double tmp1,tmp2,tmp3,tmp4;
  while(ifs >> tmp1 >> tmp2 >> tmp3 >> tmp4){
    theta_cm.push_back(tmp1);
    crs_cm.push_back(tmp2);
    theta_lab.push_back(dCmToLab(tmp1));//cm to lab
    crs_lab.push_back(cCmToLab(tmp1,tmp2));//cross section at lab angle
    //    ofs << deg_lab[i]<<" "<<crssec_lab[i]<<endl;
    //    cout<<ii<<":"<<theta_lab[ii]<<":"<<crs_lab[ii]<<endl;
   }
  /**  reslut data ***/
  vector<double> theta_rcl;//recoil degree (about 30~85)
  vector<double> crs_rcl;//calculated cross section from the nearest values
  vector<double> phi_rcl;//phi covered with 1 RDC at each theta
  vector<double> n_rcl; //number of recoil particle (proton)

  int nn = 100;
  //  int nn = (int)(max_theta - min_theta)/dev_theta +1;
  vector<vector<double>> phi_min(nn);
  vector<vector<double>> phi_max(nn);
  
  int i=0;
  double tmp_theta=0;
  double n_tmp_theta; //search near value of tmp_theta

  while(tmp_theta<max_theta){
    tmp_theta = min_theta + dev_theta*i;//theta which we want to calculate
    theta_rcl.push_back(tmp_theta);

    /* calc. number of detected protons at each theta*/    
    //search the nearest point from the input data set
    double tmp_diff = theta_lab.front() - tmp_theta;//under
    for(int j=1;j<theta_lab.size();j++){
      double x = theta_lab.at(j) - tmp_theta;
      if(x>=0 && x < tmp_diff){
	tmp_diff = x;
	n_tmp_theta = j;
      }else{break;}
    }
    if(n_tmp_theta==0){cerr<<"ERROR @ theta search"<<endl;return -1;}

    //decide cross section by taking internally dividing piont 
    double x1,x2,y;
    x1 = tmp_theta - theta_lab.at(n_tmp_theta+1);
    x2 = theta_lab.at(n_tmp_theta) - tmp_theta;
    y = (x2*crs_lab.at(n_tmp_theta+1) + x1*crs_lab.at(n_tmp_theta))
      /(x1+x2);
    crs_rcl.push_back(y);
    //    ofs << theta_rcl[i]<<" "<<crs_rcl[i]<<endl;


    /* calculating Phi covered with 1 RDC */
    int iphi=0; int flug_phi = -1;//flug_phi 0:out of RDC, 1:inside of RDC
    double tmp_phi=0;
    double yph,zph;//(y,z) from target. z=beam
    double tph; //parameter for three dimentional line of proton's trajectory
    flug_phi = 0;
    while(tmp_phi<max_phi){
      tmp_phi = min_phi + dev_phi*iphi;
      //(x,y,z) = (0,0,0)+tph(sin(the)cos(phi),sin(the)sin(phi),cos(the))
      //x->t-parameter
      tph = RDCx/sin(tmp_theta/180.*M_PI)/cos(tmp_phi/180.*M_PI) ;
      //      yph = tph*sin(tmp_theta/180.*M_PI)*sin(tmp_phi/180.*M_PI);
      yph = RDCx*tan(tmp_phi/180.*M_PI);
      zph = tph*cos(tmp_theta/180.*M_PI);

      //whether inside or outside of RDC 
      if(RDCymin<=yph&&yph<=RDCymax && RDCzmin<=zph&&zph<=RDCzmax){
    	if(flug_phi==0){
	  cout<<iphi<<":"<<tph<<" "<<yph<<" "<<zph<<endl;
    	  phi_min.at(i).push_back(tmp_phi);
    	  flug_phi=1;
    	}else flug_phi=1;
      }else{
    	if(flug_phi==1){//outside of RDC & fisrt outside point 
    	  phi_max.at(i).push_back(tmp_phi-dev_phi);
    	  flug_phi=0;
    	}else {
    	  flug_phi=0;
    	}
      }
      iphi++;
    }
    phi_rcl.emplace_back();
    for(int j=0;j<phi_min.at(i).size();j++){
      phi_rcl.at(i) += phi_max.at(i).at(j)-phi_min.at(i).at(j);
    }
    //    cout<<tmp_theta<<":"<<phi_rcl.at(i)<<endl;
    // if(phi_min.at(i).size()>=1 && phi_max.at(i).size()>=1){
    //   cout<<nn<<":"<<i<<":"<<tmp_theta<<" : "<<phi_min.at(i).at(0)
    // 	  <<" "<<phi_max.at(i).at(0)<<endl;
    // }
    i++;
  }



  
  ofs<<"# degree  Ncount (Nbeam="<<NUM_BEAM<<endl;
  double n_tgt;
  //  double phi = 2*M_PI/5.;
  n_tgt = dens_sht * thick_sht*sqrt(2) / 2. *2. * NUM_AVO;
  for(int i=0;i<crs_rcl.size()-1;i++){
    double tmp_n;
    tmp_n = NUM_BEAM * n_tgt
      * ((crs_rcl.at(i)*sin(theta_rcl.at(i)/180.*M_PI)
	  + crs_rcl.at(i+1)*sin(theta_rcl.at(i+1)/180.*M_PI))*1.E-27
	 *dev_theta/180.*M_PI /2.)
      * (phi_rcl.at(i)+phi_rcl.at(i+1))/2./180.*M_PI*2.
      * eff_det;
    // tmp_n = NUM_BEAM * n_tgt
    //   * ((crs_rcl.at(i) + crs_rcl.at(i+1))*1.E-27
    // 	 *dev_theta/2.*sin((theta_rcl.at(i)+theta_rcl.at(i+1))/2./180.*M_PI))
    //   * (phi_rcl.at(i)+phi_rcl.at(i+1))/2./180.*M_PI*2.
    //   * eff_det;
    n_rcl.push_back(tmp_n);
    ofs<<theta_rcl.at(i)+dev_theta/2.<<" "<<n_rcl.at(i)<<endl;
    //    cout<<n_rcl.at(i)<<endl;
  }
  cout<<"  ----- RESULT   ----------  "<<endl;
  cout<<"theta(lab)  ProtonNumber"<<endl;
  for(int i=0;i<n_rcl.size();i++){
    cout<<theta_rcl.at(i)+dev_theta/2.<<" "<<n_rcl.at(i)<<endl;
  }
  
  
  ifs.close();// ofs.close();
  
  return 0;

}
