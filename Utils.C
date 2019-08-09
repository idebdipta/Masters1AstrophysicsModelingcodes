#include <iostream>
#include <algorithm> 
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>

#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TArrow.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

//Utilities without any physics/////////////////////////////////////////////////////////////////////////////////////////////
double utils_Integral(std::vector<double> vx,std::vector<double> vy){
	double I=0.;

	int n=vx.size();
	int ny=vy.size();
	if(ny!=n) std::cout<<"Problem in the integration: x and y don't have the same size"<<std::endl;
	else{
		for(int i=0;i<n-1;i++) I+=(vx.at(i+1)-vx.at(i))*(vy.at(i+1)+vy.at(i))/2.;
	}
	return I;
}

double utils_proba_2_sigma(double proba){
	return abs(sqrt(2)*TMath::ErfcInverse(1.-proba));
}

//Fill a vector from xmin to xmax = xmin + n*dx/////////////////////////////////
std::vector<double> utils_InitializeVector(double xmin, double xmax, double dx){
	std::vector<double> vRes;
	int n = 1+ceil((xmax-xmin)/dx);
	for(int i=0; i<n; i++) vRes.push_back(xmin+i*dx);
	return vRes;
}

//Computes the average//////////////////////////////////////////////////////////
double utils_Average(std::vector< double > vx){
	return std::accumulate(vx.begin(),vx.end(),0.)/vx.size();
}

//Computes the variance/////////////////////////////////////////////////////////
double utils_Variance(std::vector< double > vx, double mean){
	return (std::inner_product(vx.begin(), vx.end(), vx.begin(), 0.0)/vx.size()-mean*mean);
}

double utils_Variance(std::vector< double > vx){
	return utils_Variance(vx,utils_Average(vx));
}

//Used in interpolations to find the element just before x/////////////////////
int utils_GetIndex(double x, std::vector< double > vx){

	std::vector< double >::iterator it_lower_bound = std::lower_bound (vx.begin(), vx.end(), x);//first element larger than x
	if(it_lower_bound==vx.begin())
	{
		it_lower_bound++;
		double relative_distance = (vx.front()-x)/(vx.back()-vx.front())*(vx.size()-1.);
		if(relative_distance>0.1) std::cout<<"Value out of range (>xmax) in interpolation - extrapolation performed, distance to xmin / step = "<<relative_distance*100<<"%"<< std::endl;
	}
	else if(it_lower_bound==vx.end())
	{
		double relative_distance = (x-vx.back())/(vx.back()-vx.front())*(vx.size()-1.);
		if(relative_distance>0.1) std::cout<<"Value out of range (>xmax) in interpolation - extrapolation performed, distance to xmax / step = "<<relative_distance*100<<"%"<< std::endl;
	}
	int iup = it_lower_bound-vx.begin();
	if(iup>=(int)vx.size()) iup=vx.size()-1.;//In case of off-bound events

	return iup;
}

//Interpolation tools////////////////////////////////////////////////////////////
//Interpolate a vy, vx, in x0
double utils_Interpolate(std::vector<double> vy, std::vector<double> vx, double x){
	int iup = utils_GetIndex(x,vx);
	return vy[iup-1] + (vy[iup]-vy[iup-1])*(x-vx[iup-1])/(vx[iup]-vx[iup-1]);
}	

//Interpolate a function = vector< vector < > >[i_x][i_y] at x
std::vector<double> utils_Interpolate(std::vector< std::vector<double> > vData, std::vector<double> vx, double x){
	std::vector<double> vInterpol;
	int iup = utils_GetIndex(x,vx);
	for(int j=0; j<(int)vData[iup].size(); j++)	vInterpol.push_back(vData[iup-1][j] + (vData[iup][j]-vData[iup-1][j])*(x-vx[iup-1])/(vx[iup]-vx[iup-1]));
	return vInterpol;
}	

//Interpolate a z[i_x][i_y] at x, y
double utils_Interpolate(std::vector< std::vector<double> > vz, std::vector<double> vx, double x, std::vector<double> vy, double y){
	std::vector<double> vInterpol;
	int iupx = utils_GetIndex(x,vx);
	int iupy = utils_GetIndex(y,vy);

	double res = vz[iupx-1][iupy-1]*(vx[iupx]-x)*(vy[iupy]-y);
				 res+= vz[iupx][iupy-1]*(x-vx[iupx-1])*(vy[iupy]-y);
				 res+= vz[iupx-1][iupy]*(vx[iupx]-x)*(y-vy[iupy-1]);
				 res+= vz[iupx][iupy]*(x-vx[iupx-1])*(y-vy[iupy-1]);
	res /= (vx[iupx]-vx[iupx-1])*(vy[iupy]-vy[iupy-1]);

	return res;
}	


//Load a TH1D from a vector/////////////////////////////////////////////////
TH1D *utils_LoadHisto(std::string hname, int nbin, std::vector< double > vx){

	double xmin =*std::min_element(vx.begin(),vx.end());
	double xmax =*std::max_element(vx.begin(),vx.end());
	double dx = (xmax-xmin)/(nbin-1.);

	TH1D *h = new TH1D(hname.c_str(),"",nbin+2,xmin-1.5*dx,xmax+1.5*dx);
	for(int i=0; i<(int)vx.size(); i++) h->Fill(vx[i]);

	return h;
}


//Load a TGraph from two vectors/////////////////////////////////////////////////
TGraph *utils_LoadGraph(std::vector< double > vx, std::vector< double > vy){

	if(vx.size()!=vy.size()){
		std::cout<<"Not the same number of x and y points - return an empty graph"<<std::endl;
		TGraph *G = new TGraph();
		return G;
	}

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	for(int i=0; i<n; i++){
		x[i]=vx[i];
		y[i]=vy[i];
	}
	TGraph *G = new TGraph(n,x,y);

	return G;
}

//Load a TGraphErrors from three vectors/////////////////////////////////////////
TGraphErrors *utils_LoadGraphErrors(std::vector< double > vx, std::vector< double > vy, std::vector< double > vey){

	if(vx.size()!=vy.size() || vx.size()!=vey.size()){
		std::cout<<"Not the same number of x and y points - return an empty graph errors"<<std::endl;
		TGraphErrors *G = new TGraphErrors();
		return G;
	}

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *ey = new double[n];
	for(int i=0; i<n; i++){
		x[i]=vx[i];
		y[i]=vy[i];
		ey[i]=vey[i];
	}
	TGraphErrors *G = new TGraphErrors(n,x,y,0,ey);

	return G;
}

//Load a TGraphErrors from four vectors//////////////////////////////////////////
TGraphErrors *utils_LoadGraphErrors(std::vector< double > vx, std::vector< double > vy, std::vector< double > vex, std::vector< double > vey){

	if(vx.size()!=vy.size() || vx.size()!=vex.size() || vx.size()!=vey.size()){
		std::cout<<"Not the same number of x and y points - return an empty graph errors"<<std::endl;
		TGraphErrors *G = new TGraphErrors();
		return G;
	}

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *ex = new double[n];
	double *ey = new double[n];
	for(int i=0; i<n; i++){
		x[i]=vx[i];
		y[i]=vy[i];
		ex[i]=vex[i];
		ey[i]=vey[i];
	}
	TGraphErrors *G = new TGraphErrors(n,x,y,ex,ey);

	return G;
}


//Load a TGraphErrors from four vectors/////////////////////////////////////////
TGraphAsymmErrors *utils_LoadGraphAsymmErrors(std::vector< double > vx, std::vector< double > vy, std::vector< double > veyl, std::vector< double > veyu){

	if(vx.size()!=vy.size() || vx.size()!=veyl.size()  || vx.size()!=veyu.size()){
		std::cout<<"Not the same number of x and y points - return an empty graph asymm"<<std::endl;
		TGraphAsymmErrors *G = new TGraphAsymmErrors();
		return G;
	}

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *eyl = new double[n];
	double *eyu = new double[n];
	for(int i=0; i<n; i++){
		x[i]=vx[i];
		y[i]=vy[i];
		eyl[i]=veyl[i];
		eyu[i]=veyu[i];
	}
	TGraphAsymmErrors *G = new TGraphAsymmErrors(n,x,y,0,0,eyl,eyu);

	return G;
}


//Load a TGraphErrors from six vectors//////////////////////////////////////////
TGraphAsymmErrors *utils_LoadGraphAsymmErrors(std::vector< double > vx, std::vector< double > vy, std::vector< double > vexl, std::vector< double > vexu, std::vector< double > veyl, std::vector< double > veyu){

	if(vx.size()!=vy.size() 
		|| vx.size()!=vexl.size()  || vx.size()!=vexu.size()
		|| vx.size()!=veyl.size()  || vx.size()!=veyu.size()){
		std::cout<<"Not the same number of x and y points - return an empty graph asymm"<<std::endl;
		TGraphAsymmErrors *G = new TGraphAsymmErrors();
		return G;
	}

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *exl = new double[n];
	double *exu = new double[n];
	double *eyl = new double[n];
	double *eyu = new double[n];
	for(int i=0; i<n; i++){
		x[i]=vx[i];
		y[i]=vy[i];
		exl[i]=vexl[i];
		exu[i]=vexu[i];
		eyl[i]=veyl[i];
		eyu[i]=veyu[i];
	}
	TGraphAsymmErrors *G = new TGraphAsymmErrors(n,x,y,exl,exu,eyl,eyu);

	return G;
}

//Load a TGraphAsymmErrors from a 4 column file
TGraphAsymmErrors *utils_LoadGraphAsymm(char *filename, bool verbose = false){

	std::vector<double> vx, vy, veyl, veyu;
	std::ifstream myfile;	
	myfile.open(filename);
	if(myfile.is_open()){
		while(myfile.good()){
			std::string line;
			std::getline (myfile,line);
			std::stringstream ss(line);
			std::vector< double > vbuf;
			double buf;
			while(ss>>buf) vbuf.push_back(buf);
			if(((vbuf.size()<3) || (vbuf.size()>4))){
				if(verbose) std::cout<<"Line "<<line<<" unread in file "<<filename<<std::endl;
			}
			else{
				vx.push_back(vbuf[0]);
				vy.push_back(vbuf[1]);
				veyl.push_back(vbuf[2]);
				if(vbuf.size()==3) veyu.push_back(vbuf[2]);
				else if(vbuf.size()==4)	veyu.push_back(vbuf[3]);
			}
		}
		myfile.close();
	}
	else std::cout<<"Could not open file: "<<filename<<std::endl;

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *eyl = new double[n];
	double *eyu = new double[n];
	for(int i=0; i<n; i++){
		x[i] = vx[i];
		y[i] = vy[i];
		eyl[i] = veyl[i];
		eyu[i] = veyu[i];
	}	
	TGraphAsymmErrors *G = new TGraphAsymmErrors(n,x,y,0,0,eyl,eyu);

	return G;
}

//Load a TGraphAsymmErrors from a 6 column file
TGraphAsymmErrors *utils_LoadGraphAsymmWithEX(string filename, bool verbose = false){

	std::vector<double> vx, vy, vexl, vexu, veyl, veyu;
	std::ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		while(myfile.good()){
			std::string line;
			std::getline (myfile,line);
			std::stringstream ss(line);
			std::vector< double > vbuf;
			double buf;
			while(ss>>buf) vbuf.push_back(buf);
			if(vbuf.size()!=6){
				if(verbose) std::cout<<"Line "<<line<<" unread in file "<<filename<<std::endl;
			}
			else{
				vx.push_back(vbuf[0]);
				vy.push_back(vbuf[1]);
				vexl.push_back(vbuf[2]);
				vexu.push_back(vbuf[3]);
				veyl.push_back(vbuf[4]);
				veyu.push_back(vbuf[5]);
			}
		}
		myfile.close();
	}
	else std::cout<<"Could not open file: "<<filename<<std::endl;

	int n = vx.size();
	double *x = new double[n];
	double *y = new double[n];
	double *eyl = new double[n];
	double *eyu = new double[n];
	double *exl = new double[n];
	double *exu = new double[n];
	for(int i=0; i<n; i++){
		x[i] = vx[i];
		y[i] = vy[i];
		exl[i] = vexl[i];
		exu[i] = vexu[i];
		eyl[i] = veyl[i];
		eyu[i] = veyu[i];
	}	
	TGraphAsymmErrors *G = new TGraphAsymmErrors(n,x,y,exl,exu,eyl,eyu);

	return G;
}

//returns the min [0] and max [1] y points of a graph
std::vector< double > utils_min_and_max_points(TGraph *G0){
	std::vector< double > vres;
	
	double min = G0->GetY()[0];
	double max = G0->GetY()[0];
	for(int i=1; i<G0->GetN(); i++){
		double phi = G0->GetY()[i];
		if(phi>max) max = phi;
		else if(phi<min) min = phi;
	}

	vres.push_back(min);
	vres.push_back(max);
	return vres;
}

//Loads arrows for UL points
std::vector< TArrow* > utils_ArrowsUpperLimits(TGraphAsymmErrors *G0, double size = 1.5){
	std::vector< TArrow* > vRes;

	int n = G0->GetN();	    
	double *x  = G0->GetX();
	double *y = G0->GetY();		
	for(int i = 0;i<n;i++){
		if(G0->GetErrorYlow(i)==0){
			TArrow *arr = new TArrow(x[i],y[i],x[i],y[i]/size,0.01,"|>");
			arr->SetLineWidth(2); arr->SetAngle(60);
			arr->SetLineColor(G0->GetMarkerColor()); 
			arr->SetFillColor(G0->GetMarkerColor());
			vRes.push_back(arr);
		}
	}

	return vRes;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Utilities containing some physics/////////////////////////////////////////////////////////////////////////////////////////
//General
#define c_pi 3.14159265359
#define c_sigmaT 6.6524586E-25//cm2	
#define c_speed 2.99792458E10//cm/s
#define c_m_e_eV 511.E3//eV
#define c_h_eV 4.1356673E-15//eV s
#define c_eV_2_Joule 1.60217657E-19//Joule
#define c_kB_eV_div_K 8.617343E-5//eV / K

//Conversion
#define c_Mpc_2_cm 3.08568E24// cm
#define c_W_2_nW 1.E9//nW
#define c_eV_2_TeV 1E-12//TeV
#define c_um_2_cm 1E-4//cm
#define c_cm2_2_m2 1E-4//m2
#define c_m2_2_cm2 1E4//cm2
#define c_km_2_cm 1.E5//cm
#define c_TeV_2_erg 1.60217657//erg

//Convert a density to a nuInu at nu///////////////////////////////////////////////////
double utils_Density2Intensity(
													double n,//cm-3
													double nu//Hz
){
	return (c_W_2_nW*c_eV_2_Joule/c_cm2_2_m2)*(c_speed/(4.*c_pi))*(c_h_eV*nu)*n;
}

//Convert a nuInu to a density at nu///////////////////////////////////////////////////
double utils_Intensity2Density(
													double nuInu,//nW m-2 sr-1
													double nu//Hz
){
	return nuInu/((c_W_2_nW*c_eV_2_Joule/c_cm2_2_m2)*(c_speed/(4.*c_pi))*(c_h_eV*nu));
}

//Convert log(E/m_e c^2) to E[TeV]///////////////////////////////////////////////////
double utils_e2E_TeV(double e){
	return exp(e)* c_m_e_eV * c_eV_2_TeV;
}

//Convert E[TeV] to a log(E/m_e c^2)///////////////////////////////////////////////////
double utils_E_TeV2e(double E){
	return log(E/(c_m_e_eV * c_eV_2_TeV));
}

//Convert eV to 10-9 Joule/////////////////////////////////////////////////////////////

double utils_eV2nJ(double eps){
	return eps*c_eV_2_Joule*c_W_2_nW;
}

//Convert a log(E/m_e c^2) to nu[Hz]///////////////////////////////////////////////////
double utils_e2nu_Hz(double e){
	return exp(e) * c_m_e_eV / c_h_eV;
}


//Convert E[TeV] to nu[Hz]/////////////////////////////////////////////////////////////
double utils_E_TeV2nu_Hz(double E){
	return E/(c_eV_2_TeV*c_h_eV);
}

//Convert a log(E/m_e c^2) to lambda[um]///////////////////////////////////////////////
double utils_e2lambda_um(double e){
	return c_speed*c_h_eV / (exp(e) * c_m_e_eV * c_um_2_cm);
}

//Convert a lambda[um] to log(E/m_e c^2)///////////////////////////////////////////////
double utils_lambda_um2e(double lambda){
	return log(c_speed*c_h_eV / (lambda * c_m_e_eV * c_um_2_cm));
}

//Convert nu[Hz] to lambda [um]////////////////////////////////////////////////////////
double utils_nu2lambda_um(double nu){
	return c_speed/(nu*c_um_2_cm);
}

//Convert lambda [um] to a nu[Hz]////////////////////////////////////////////////////////
double utils_lambda2nu_Hz(double lambda){
	return c_speed/(lambda*c_um_2_cm);
}

//Returns the EBL intensity based on an SED stored in two vectors///////////////////
double utils_DensityEBL(double e, std::vector< double > vlambda, std::vector< double > vnuInu){
	double res = 0;
	double lambda = utils_e2lambda_um(e);
	double nu = utils_lambda2nu_Hz(lambda);
	if((lambda > 0.05) && (lambda < 1E3)) res = utils_Intensity2Density(utils_Interpolate(vnuInu,vlambda,lambda),nu) ;//Intensity nuInu
	return res;
}

//Returns the EBL intensity based on an SED stored in a TF1*////////////////////////
double utils_DensityEBL(double e, TF1* fintensity){
	double res = 0;
	double lambda = utils_e2lambda_um(e);
	if((lambda > 0.05) && (lambda < 1E3)) res = utils_Intensity2Density(fintensity->Eval(lambda),utils_lambda2nu_Hz(lambda)) ;//Intensity nuInu
	return res;
}

//Finds the the maximum EBL wavelength//////////////////////////////////////////////
double utils_MaximumWavelength(std::vector< double > vEmax, std::vector< double > vz){
	double lambda_max = 0;

	if(vz.size()!=vEmax.size())	std::cout<<"Problem computing the maximum wavelength, not the same number of redshifts and spectra!"<<std::endl;
	else{
		std::vector<double> vlambda;
		for(int i=0; i<(int)vz.size(); i++){
			vlambda.push_back(utils_e2lambda_um(-2.*log(1.+vz[i])-utils_E_TeV2e(vEmax[i])));
//			std::cout<<"Maximum lambda: "<<vlambda.back()<<std::endl;
		}
		lambda_max = *std::max_element(vlambda.begin(),vlambda.end());
	}

	return lambda_max;
}

//Finds the minimum     EBL wavelength//////////////////////////////////////////////
double utils_MinimumWavelength(double lambda0, double lambda_max, double dln_lambda){
	int n0 = round(log(lambda_max/lambda0)/dln_lambda);
	return lambda_max* exp(-n0*dln_lambda);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Gilmore+12 SED//////////////////////////////////////////////////////
std::vector< std::vector< double > > utils_LoadG12_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/Intensity_z=0_Gilmore2012_fixed", bool verbose = false){
	std::vector< std::vector< double > > vRes;
	vRes.push_back(std::vector< double > ());//vLambda
	vRes.push_back(std::vector< double > ());//vnuInu

	std::ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		while(myfile.good()){
			std::string line;
			std::getline (myfile,line);
			std::stringstream ss(line);
			std::vector< double > vbuf;
			double buf;
			while(ss>>buf) vbuf.push_back(buf);
			if((vbuf.size()!=2)){
				if(verbose) std::cout<<"Line "<<line<<" unread in file "<<filename<<std::endl;
			}
			else{
				vRes[0].push_back(vbuf[0]);
				vRes[1].push_back(vbuf[1]);
			}
		}
		myfile.close();
	}
	else std::cout<<"Could not open file: "<<filename<<std::endl;

	return vRes;
}

//Khaire and Srianand 2014
std::vector< std::vector< double > > utils_LoadKS_mw_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KhaireSrianand_EBL_z0/SpectraKS_gal_dustmed_mw_0.00_proper_units.txt"){	return utils_LoadG12_EBL_SED(filename); }
std::vector< std::vector< double > > utils_LoadKS_cal_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KhaireSrianand_EBL_z0/SpectraKS_gal_dustmed_cal_0.00_proper_units.txt"){	return utils_LoadG12_EBL_SED(filename); }
std::vector< std::vector< double > > utils_LoadKS_lmc_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KhaireSrianand_EBL_z0/SpectraKS_gal_dustmed_lmc_0.00_proper_units.txt"){	return utils_LoadG12_EBL_SED(filename); }
std::vector< std::vector< double > > utils_LoadKS_lmc2_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KhaireSrianand_EBL_z0/SpectraKS_gal_dustmed_lmc2_0.00_proper_units.txt"){	return utils_LoadG12_EBL_SED(filename); }
std::vector< std::vector< double > > utils_LoadKS_smc_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KhaireSrianand_EBL_z0/SpectraKS_gal_dustmed_smc_0.00_proper_units.txt"){	return utils_LoadG12_EBL_SED(filename); }

//Finke 2010
std::vector< std::vector< double > > utils_LoadFi10_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/Finke2010_EBL_modelC_total_intensity_z0.00.csv"){	return utils_LoadG12_EBL_SED(filename); }

//Dominguez 2011
std::vector< std::vector< double > > utils_LoadD11_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/Dominguez_2011"){	return utils_LoadG12_EBL_SED(filename); }

//Kneiske & Dole 2010 (MRF LL) // 2010A&A...515A..19K
std::vector< std::vector< double > > utils_LoadKD10_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/KneiskeDole"){	return utils_LoadG12_EBL_SED(filename); }

//Stecker 2006 //2006ApJ...648..774S
std::vector< std::vector< double > > utils_LoadS06_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/Stecker"){	return utils_LoadG12_EBL_SED(filename); }

//Franceschini+08 SED//////////////////////////////////////////////////////
std::vector< std::vector< double > > utils_LoadFR08_EBL_SED(std::string filename = "../../Data/EBL_data/Other_models/Intensity_z=0_Franceschini2008", bool verbose = false){
	std::vector< std::vector< double > > vRes;
	vRes.push_back(std::vector< double > ());//vLambda
	vRes.push_back(std::vector< double > ());//vnuInu

	std::ifstream myfile;	
	myfile.open(filename.c_str());
	if(myfile.is_open()){
		while(myfile.good()){
			std::string line;
			std::getline (myfile,line);
			std::stringstream ss(line);
			std::vector< double > vbuf;
			double buf;
			while(ss>>buf) vbuf.push_back(buf);
			if((vbuf.size()!=2)){
				if(verbose) std::cout<<"Line "<<line<<" unread in file "<<filename<<std::endl;
			}
			else{
				double nu = pow(10.,vbuf[0])/c_h_eV;
				vRes[0].push_back(utils_nu2lambda_um(nu));
				vRes[1].push_back(utils_Density2Intensity(pow(10.,vbuf[1]),nu));
			}
		}
		myfile.close();
	}
	else std::cout<<"Could not open file: "<<filename<<std::endl;

	return vRes;
}

//Franceschini+08 optical depth///////////////////////////////////////////////
std::vector< std::vector< double > > utils_LoadFR08_OpticalDepth(std::vector<double> vz0, std::vector<double> vEgamma0, string filename = "../../Data/EBL_data/Other_models/OpticalDepth_z<1_Franceschini2008", bool verbose = false){

	std::vector< double > vz;
	std::vector< std::vector< std::vector< double > > > vEBL;

  std::ifstream myfile;
	myfile.open(filename.c_str());
  if(myfile.is_open()){
    while(myfile.good()){
		  std::string line;
      getline (myfile,line);
			std::stringstream ss(line);
			double buf;
			std::vector< double > vbuf;
			while(ss>>buf) vbuf.push_back(buf);
			if(vbuf.size()==1){
				vEBL.push_back(std::vector< std::vector< double > >());
				vz.push_back(vbuf[0]);
      } 
			else if(vbuf.size()==4){
				if(vEBL.back().size()==0){
					vEBL.back().push_back(std::vector< double > ());
					vEBL.back().push_back(std::vector< double > ());
				}
				vEBL.back()[0].push_back(vbuf[0]);
				vEBL.back()[1].push_back(vbuf[2]);
			} 
			else if(verbose) std::cout<<"Line unread in EBL file: "<<line<<std::endl;
    }
    myfile.close();
  }
	else std::cout<<"EBL File not found"<<std::endl;

	std::vector< std::vector< double > > vOpacity;
	for(int k=0; k<(int)vz0.size(); k++){//redshift
		vOpacity.push_back(std::vector< double >());
		for(int j=0; j<(int)vEgamma0.size(); j++){//energy
			int iz = utils_GetIndex(vz0[k],vz);//Given that vz is at 10-3, no need to interpolate
			vOpacity[k].push_back(utils_Interpolate(vEBL[iz][1],vEBL[iz][0],vEgamma0[j]));
		}
	}	

	return vOpacity;
}


//Gilmore+12 optical depth///////////////////////////////////////////////
std::vector< std::vector< double > > utils_LoadG12_OpticalDepth(std::vector<double> vz0, std::vector<double> vEgamma0, string filename = "../../Data/EBL_data/Other_models/OpticalDepth_All_Gilmore2012_fixed", bool verbose = false){
	std::vector< std::vector< double > > vOpacity;

	std::vector< double > vz;
	std::vector< double > vE;
	std::vector< std::vector< double > > vtau;

  std::ifstream myfile;
	myfile.open(filename.c_str());
  if(myfile.is_open()){
    while(myfile.good()){
		  std::string line;
      getline (myfile,line);
			std::stringstream ss(line);
			double buf;
			if(vz.size()==0) while(ss>>buf) vz.push_back(buf);
			else{
				std::vector< double > vbuf;
				while(ss>>buf) vbuf.push_back(buf);
				if(vbuf.size()==53){
					vE.push_back(vbuf.front());	
					vbuf.erase(vbuf.begin());
					vtau.push_back(vbuf);
				}
				else if(verbose) std::cout<<"Line unread in EBL file: "<<line<<std::endl;
      } 
		}
    myfile.close();
	}
	else std::cout<<"EBL File not found"<<std::endl;

	for(int k=0; k<(int)vz0.size(); k++){//redshift
		vOpacity.push_back(std::vector< double >());
		for(int j=0; j<(int)vEgamma0.size(); j++){//energy
			vOpacity[k].push_back(utils_Interpolate(vtau,vE,vEgamma0[j],vz,vz0[k]));
		}
	}	

	return vOpacity;
}


//Dominguez & Prada 2013//////////////////////////////////////////////////
std::vector< std::vector< double > > utils_LoadD13_OpticalDepth(std::vector<double> vz0, std::vector<double> vEgamma0, string filename = "../../Data/EBL_data/Other_models/DominguezPrada13_H0/hr_tau_hr_ebl_dominguez11_h0.700_WM0.28_w0-1.00.out", bool verbose = false){
	std::vector< std::vector< double > > vOpacity;

	std::vector< double > vz;
	std::vector< double > vE;
	std::vector< std::vector< double > > vtau;

  std::ifstream myfile;
	myfile.open(filename.c_str());
  if(myfile.is_open()){
    while(myfile.good()){
		  std::string line;
      getline (myfile,line);
			std::stringstream ss(line);
			double buf;
			if(vz.size()==0) while(ss>>buf) vz.push_back(buf);
			else{
				std::vector< double > vbuf;
				while(ss>>buf) vbuf.push_back(buf);
				if(vbuf.size()==50){
					vE.push_back(vbuf.front());	
					vbuf.erase(vbuf.begin());
					vtau.push_back(vbuf);
				}
				else if(verbose) std::cout<<"Line unread in EBL file: "<<line<<std::endl;
      } 
		}
    myfile.close();
	}
	else std::cout<<"EBL File not found"<<std::endl;

	for(int k=0; k<(int)vz0.size(); k++){//redshift
		vOpacity.push_back(std::vector< double >());
		for(int j=0; j<(int)vEgamma0.size(); j++){//energy
			vOpacity[k].push_back(utils_Interpolate(vtau,vE,vEgamma0[j],vz,vz0[k]));
//			std::cout<<vE.front()<<" "<<vEgamma0[j]<<" "<<vE.back()<<" - "<<vz.front()<<" "<<vz0[k]<<" "<<vz.back()<<std::endl;
		}
	}	

	return vOpacity;
}
