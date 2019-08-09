#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <iostream>
#include <sstream>

#include <string>
#include <vector>

#include <cmath>
#include <vector>
#include <fstream>

#include <TCanvas.h>
#include "TLegend.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLine.h"

#include "../../Code/Base/Utils.C"
#include "../../Code/Base/Model.C"

//Conversion 

#define c_eV_2_TeV 1E-12//TeV
#define c_um_2_cm 1E-4//cm
#define c_TeV_2_eV 1E12//eV


using namespace std;

//Convert lambda [um] to a e_eV////////////////////////////////////////////////////////

double utils_lambda2e_eV(double lambda)
{
	return c_speed*c_h_eV/(lambda*c_um_2_cm); //output in eV
}

double eGamma(double eGamma_eV)

{
	return eGamma_eV/c_TeV_2_eV;
}

double utils_lambda_um_E_TeV_2_e0(double lambda, double eGamma)//um, TeV
{
	
	double e_eV=utils_lambda2e_eV(lambda);
	double eGamma_eV=eGamma*c_TeV_2_eV;
	double e0= e_eV*eGamma_eV/(c_m_e_eV*c_m_e_eV); //e_not is in unit of (cm*TeV)/(eV*um)  

	return log(e0);
	return e0;
}

double ve(double e, double eGamma)
{
	double c_tmp_m=c_m_e_eV*c_m_e_eV;
	double c_tmp_cspeed=c_speed*c_speed*c_speed*c_speed;
	double ve= log((e*eGamma)/(c_tmp_m*c_tmp_cspeed));
	return ve;
}

//Convert TeV into log(E/m_ec2)
double veGamma(double eGamma)//TeV
{
	double veGamma=log(eGamma*c_TeV_2_eV/c_m_e_eV);
	return veGamma;
}


/////////////////////////////////////////////////////////////////////////////////////
//Opacity as a function of veGamma[i] = log(E/mc^2) based on a template SED: vnu, vnuInu//
std::vector< double > Compute_Opacity(
	std::vector<double> vlambda, 
	std::vector<double> vnuInu,
	std::vector<double> vKernel_fixed_z,  //OmegaM, OmegaLambda, f_evol
	double H0,//km s-1 Mpc-1
	std::vector< double > ve,             //ve[i] = log(epsilonEBL*E/m2c4) with energies as of today
	std::vector< double > veGamma        //veGamma[i] = log(Ei/mc^2) where Ei is the gamma-ray energy
){
	std::vector< double > vOpacity;

	double R_H = c_speed/(H0*c_km_2_cm/c_Mpc_2_cm);//cm
	double normalization = 0.75*c_sigmaT*R_H;
	for(int i=0; i<(int)veGamma.size(); i++){
		double integral_over_e = 0;
		double step_e = (ve.back()-ve.front())/ve.size();
		for(int j=0; j<(int)ve.size(); j++){//Integration over log(e0)
			double nEBL = utils_DensityEBL(ve[j]-veGamma[i], vlambda, vnuInu);
			double integrand = normalization*nEBL*exp(-2*ve[j])*vKernel_fixed_z[j]*step_e;
			if((j==0) || (j==(int)ve.size()-1)) integrand *= 0.5;//Trapezoidal integration with a uniform grid in ve
			integral_over_e += integrand;
		}
		vOpacity.push_back(integral_over_e);
	}

	return vOpacity;
}


//Opacity as a function of veGamma[i] = log(E/mc^2) based on a template SED: vnu, vnuInu//

//Returning to initials for different f_evol, which is two different Kernel basically

std::vector< double > Compute_Opacity12(
	std::vector<double> vlambda, 
	std::vector<double> vnuInu,
	std::vector<double> vKernel12_fixed_z1,  //OmegaM, OmegaLambda, f_evol1
	std::vector<double> vKernel12_fixed_z2,  //OmegaM, OmegaLambda, f_evol2
	double H0,//km s-1 Mpc-1
	std::vector< double > ve,		//ve[i] = log(epsilonEBL*E/m2c4) with energies as of today
	std::vector< double > veGamma,		//veGamma[i] = log(Ei/mc^2) where Ei is the gamma-ray energy
	double lambda0
){
	std::vector<double> vOpacity12;

	double R_H = c_speed/(H0*c_km_2_cm/c_Mpc_2_cm);//cm
	double normalization = 0.75*c_sigmaT*R_H;
	for(int i=0; i<(int)veGamma.size(); i++){
		double Egamma = c_m_e_eV*exp(veGamma[i])/c_TeV_2_eV; //Gamma-energy in TeV
		double e0 = utils_lambda_um_E_TeV_2_e0(lambda0, Egamma); 

		double integral_over_e = 0;
		double step_e = (ve.back()-ve.front())/ve.size();
		for(int j=0; j<(int)ve.size(); j++){//Integration over log(e0)
			double nEBL = utils_DensityEBL(ve[j]-veGamma[i], vlambda, vnuInu);
			double kernel12 = 0.;
			if(ve[j]>e0) kernel12 = vKernel12_fixed_z1[j];
			else kernel12 = vKernel12_fixed_z2[j];

			double integrand = normalization*nEBL*exp(-2*ve[j])*kernel12*step_e;
			if((j==0) || (j==(int)ve.size()-1)) integrand *= 0.5;//Trapezoidal integration with a uniform grid in ve
			integral_over_e += integrand;
		}
		vOpacity12.push_back(integral_over_e);
	}

	return vOpacity12;
}

/////////////////////////////////////////////////////////////////////////////////////
//root function returning the total chi2/////////////////////////////////////////////
void ShowOpticalDepth( 
	double lambda0 = 10,//10um
	double ztest1 = 0.05,//redshift at which the optical depth is evaluated
 	double ztest2 = 0.2,
	double ztest3 = 1,
	double f_evol = 2,//Evolution parameter
	double f_evol1_1 = 0,//Evolution parameter
	double f_evol1_2 = 2,//Evolution parameter
	double f_evol1_3 = 4,//Evolution parameter
	double f_evol2 = 2,//Evolution parameter
	double zmax = 2.0//Maximum redshift up to which the optical depth is computed
){

	//////////// INITIALIZATION /////////////////////////////////////////////////////////////

	//Define the cosmological constants
	double H0 = 70;//Hubble constant in km s-1 Mpc-1
	double omegaM = 0.3;//Omega Matter
	double omegaLambda = 1.-omegaM;//Omega Lambda

	//Loads the EBL model at z=0
	std::vector< std::vector< double > > vSED = utils_LoadFR08_EBL_SED();
	std::string sname_model = "Franceschini 2008";
	std::vector< double > vlambda = vSED[0];//lambda [micro-meter]
	std::vector< double > vnuInu = vSED[1];//\nuI\nu [nW m-2 sr-1] 

	/* Alternative models
	utils_LoadG12_EBL_SED(); //-> Gilmore 2012
	"Gilmore 2012";

	utils_LoadD11_EBL_SED(); //-> Dominguez 2011
	"Dominguez 2011";

	utils_LoadKD10_EBL_SED(); //-> Kneiske & Dole 2010
	"Kneiske & Dole 2010";

	utils_LoadKS_lmc2_EBL_SED(); //-> Khaire and Srianand 2014
	"Khaire and Srianand 2014 - lmc2";

	utils_LoadFi10_EBL_SED(); //-> Finke 2010
	"Finke 2010";

	utils_LoadS06_EBL_SED(); //-> Kneiske & Dole 2010
	"Stecker 2006";
	*/
	/////////////////////////////////////////////////////////////////////////////////////////


	//////////// COMPUTATION ////////////////////////////////////////////////////////////////
	//Initialize the Kernel for computation of the optical depth up to zmax
	std::vector<double> veGamma, ve, vzsteps;//used for the integration
	std::vector< std::vector<double> > vKernel;//Kernel defined in Biteau and Williams 2015

	std::vector<double> vEgamma;//Gamma-ray energy
	Model_InitializeKernel(omegaM, omegaLambda, f_evol, zmax, veGamma, vEgamma, ve,  vzsteps, vKernel);

//////////// COMPUTATION ////////////////////////////////////////////////////////////////
	//Initialize the Kernel12 for computation of the optical depth up to zmax


	std::vector< std::vector<double> > vKernel1_1;//Kernel12 defined in Biteau and Williams 2015
	Model_InitializeKernel(omegaM, omegaLambda, f_evol1_1, zmax, veGamma, vEgamma, ve,  vzsteps, vKernel1_1);

	std::vector< std::vector<double> > vKernel1_2;//Kernel12 defined in Biteau and Williams 2015
	Model_InitializeKernel(omegaM, omegaLambda, f_evol1_2, zmax, veGamma, vEgamma, ve,  vzsteps, vKernel1_2);

	std::vector< std::vector<double> > vKernel1_3;//Kernel12 defined in Biteau and Williams 2015
	Model_InitializeKernel(omegaM, omegaLambda, f_evol1_3, zmax, veGamma, vEgamma, ve,  vzsteps, vKernel1_3);




	std::vector< std::vector<double> > vKernel2;//Kernel12 defined in Biteau and Williams 2015
	Model_InitializeKernel(omegaM, omegaLambda, f_evol2, zmax, veGamma, vEgamma, ve,  vzsteps, vKernel2);



	//Load the optical depth at ztest for fevol***********************
	std::vector< double > vOpacity_1 = Compute_Opacity(vlambda,vnuInu, utils_Interpolate(vKernel,vzsteps,ztest1),H0, ve, veGamma);//Optical depth vs z vs E
std::vector< double > vOpacity_2 = Compute_Opacity(vlambda,vnuInu, utils_Interpolate(vKernel,vzsteps,ztest2),H0, ve, veGamma);//Optical depth vs z vs E
std::vector< double > vOpacity_3 = Compute_Opacity(vlambda,vnuInu, utils_Interpolate(vKernel,vzsteps,ztest3),H0, ve, veGamma);//Optical depth vs z vs E
	/////////////////////////////////////////////////////////////////////////////////////////****************end-of-operation


	//Load the optical depth at ztest for fevol1**********************************

	//z=z1=0.05//////////////////////////////vOpacity12 for fevol1 or fevol2, (_1)= fevol1, _2=fevol2 and so on, the last term implies the z =1,2,3 values. at the last, ev1 is fevol 1
	std::vector< double > vOpacity12_11ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_1,vzsteps,ztest1), utils_Interpolate(vKernel2,vzsteps,ztest1),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_21ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_2,vzsteps,ztest1), utils_Interpolate(vKernel2,vzsteps,ztest1),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_31ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_3,vzsteps,ztest1), utils_Interpolate(vKernel2,vzsteps,ztest1),H0, ve, veGamma, lambda0);


	//z=z1=0.2//////////////////////////////
	std::vector< double > vOpacity12_12ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_1,vzsteps,ztest2), utils_Interpolate(vKernel2,vzsteps,ztest2),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_22ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_2,vzsteps,ztest2), utils_Interpolate(vKernel2,vzsteps,ztest2),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_32ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_3,vzsteps,ztest2), utils_Interpolate(vKernel2,vzsteps,ztest2),H0, ve, veGamma, lambda0);


	//z=z1=1//////////////////////////////
	std::vector< double > vOpacity12_13ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_1,vzsteps,ztest3), utils_Interpolate(vKernel2,vzsteps,ztest3),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_23ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_2,vzsteps,ztest3), utils_Interpolate(vKernel2,vzsteps,ztest3),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_33ev1 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel1_3,vzsteps,ztest3), utils_Interpolate(vKernel2,vzsteps,ztest3),H0, ve, veGamma, lambda0);//Optical depth vs z vs E
	/////////////////////////////////////////////////////////////////////////////////////////****************end-of-operation fevol1 (3 values of fevol1)




//**************************************
		//Load the optical depth at ztests for fevol2
	std::vector< double > vOpacity12_1ev2 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel2,vzsteps,ztest1), utils_Interpolate(vKernel2,vzsteps,ztest1),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_2ev2 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel2,vzsteps,ztest2), utils_Interpolate(vKernel2,vzsteps,ztest2),H0, ve, veGamma, lambda0);

	std::vector< double > vOpacity12_3ev2 = Compute_Opacity12(vlambda,vnuInu, utils_Interpolate(vKernel2,vzsteps,ztest3), utils_Interpolate(vKernel2,vzsteps,ztest3),H0, ve, veGamma, lambda0);//Optical depth vs z vs E
	/////////////////////////////////////////////////////////////////////



//************************************* end of evol2

//////////// PLOTTING-fevol ///////////////////////////////////////////////////////////////////
	//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E1 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity_1[0]);
		Gopt_vs_E1->SetLineWidth(4);
		Gopt_vs_E1->SetLineStyle(1);
		Gopt_vs_E1->SetLineColor(kBlue);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E2 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity_2[0]);
		Gopt_vs_E2->SetLineWidth(4);
		Gopt_vs_E2->SetLineStyle(1);
		Gopt_vs_E2->SetLineColor(kBlue);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E3 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity_3[0]);
		Gopt_vs_E3->SetLineWidth(4);
		Gopt_vs_E3->SetLineStyle(1);
		Gopt_vs_E3->SetLineColor(kBlue);

//////////// PLOTTING-fevol1 ///////////////////////////////////////////////////////////////////
	//Load the graph of the optical depth vs energy (in log-log space)


//////****fevol1_1///////////
//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E1z1 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_11ev1[0]);
		Gopt_vs_E1z1->SetLineWidth(4);
		Gopt_vs_E1z1->SetLineStyle(2);
		Gopt_vs_E1z1->SetLineColor(kGreen);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E2z1 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_21ev1[0]);
		Gopt_vs_E2z1->SetLineWidth(4);
		Gopt_vs_E2z1->SetLineStyle(10);
		Gopt_vs_E2z1->SetLineColor(kRed);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E3z1 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_31ev1[0]);
		Gopt_vs_E3z1->SetLineWidth(4);
		Gopt_vs_E3z1->SetLineStyle(8);
		Gopt_vs_E3z1->SetLineColor(kCyan);

//////****fevol1_2///////////

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E1z2 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_12ev1[0]);
		Gopt_vs_E1z2->SetLineWidth(4);
		Gopt_vs_E1z2->SetLineStyle(2);
		Gopt_vs_E1z2->SetLineColor(kGreen);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E2z2 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_22ev1[0]);
		Gopt_vs_E2z2->SetLineWidth(4);
		Gopt_vs_E2z2->SetLineStyle(10);
		Gopt_vs_E2z2->SetLineColor(kRed);


//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E3z2 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_32ev1[0]);
		Gopt_vs_E3z2->SetLineWidth(4);
		Gopt_vs_E3z2->SetLineStyle(8);
		Gopt_vs_E3z2->SetLineColor(kCyan);

//////****fevol1_3///////////

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E1z3 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_13ev1[0]);
		Gopt_vs_E1z3->SetLineWidth(4);
		Gopt_vs_E1z3->SetLineStyle(2);
		Gopt_vs_E1z3->SetLineColor(kGreen);

//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E2z3 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_23ev1[0]);
		Gopt_vs_E2z3->SetLineWidth(4);
		Gopt_vs_E2z3->SetLineStyle(10);
		Gopt_vs_E2z3->SetLineColor(kRed);


//Load the graph of the optical depth vs energy (in log-log space)
	TGraph* Gopt_vs_E3z3 = new TGraph(vEgamma.size(), &vEgamma[0], &vOpacity12_33ev1[0]);
		Gopt_vs_E3z3->SetLineWidth(4);
		Gopt_vs_E3z3->SetLineStyle(8);
		Gopt_vs_E3z3->SetLineColor(kCyan);

	//Histogram for plotting multiple curves
	TH1D *h0 =  new TH1D("h0","",1000,0.05,60.);//from 50 GeV up to 2 TeV
		h0->SetMinimum(0.02); 
		h0->SetMaximum(5.0);
		h0->GetXaxis()->SetTitleSize(0.04);
		h0->GetYaxis()->SetTitleSize(0.04);
		h0->GetXaxis()->SetLabelSize(0.04);
		h0->GetYaxis()->SetLabelSize(0.04);
		h0->GetXaxis()->SetTitle("#bf{Energy}  [TeV]");
		h0->GetXaxis()->SetTitleOffset(1.2);
		h0->GetXaxis()->CenterTitle();
		h0->GetYaxis()->SetTitle("#bf{Optical depth}  #tau");
		h0->GetYaxis()->SetTitleOffset(1.3);
		h0->GetYaxis()->CenterTitle();
		h0->SetStats(0); 
		h0->SetDirectory(0);

		TLatex tl;
		tl.SetTextAlign(12);
		tl.SetTextSize(0.04);
		tl.SetNDC();

	//Legend for different curves
	TLegend *leg1 = new TLegend(0.7,0.13,0.88,0.25);
		leg1->SetLineColor(kWhite); 
		leg1->SetFillColor(kWhite);
		leg1->SetMargin(0.3); 
		leg1->SetTextSize(0.04);

//////////for fevol and fevol1_1, fevol1_2, fevol1_3, total of 9 combination as there are also 3 z values************///////////////
	
		

		std::stringstream ssk1_1;
		ssk1_1<<"f_{evol 1} = "<<f_evol1_1;
		leg1->AddEntry(Gopt_vs_E1z1,ssk1_1.str().c_str(),"l");


		std::stringstream ssk1_2;
		ssk1_2<<"f_{evol 1} = "<<f_evol1_2;
		leg1->AddEntry(Gopt_vs_E2z1,ssk1_2.str().c_str(),"l");

		std::stringstream ssk1_3;
		ssk1_3<<"f_{evol 1} = "<<f_evol1_3;
		leg1->AddEntry(Gopt_vs_E3z1,ssk1_3.str().c_str(),"l");

/////////////////////end of fevol and fevol1 stuff////////////

	//Draw the canvas
	TCanvas *c1 = new TCanvas("c1","c1");
		c1->cd()->SetLogx();
//		c1->cd()->SetLogy();
		h0->Draw();
		leg1->Draw();

		tl.DrawLatex(0.72,0.4,"z = 0.05");
		tl.DrawLatex(0.5,0.5,"z = 0.2");
		tl.DrawLatex(0.29,0.71,"z = 1");

		tl.DrawLatex(0.4,0.94,"f_{evol 2} = 2,  Varying f_{evol 1}, #lambda_{cut} = 10 #mum");


		Gopt_vs_E1z1->Draw("same l");	
		Gopt_vs_E2z1->Draw("same l");	
		Gopt_vs_E3z1->Draw("same l");

		Gopt_vs_E1z2->Draw("same l");	
		Gopt_vs_E2z2->Draw("same l");	
		Gopt_vs_E3z2->Draw("same l");

		Gopt_vs_E1z3->Draw("same l");	
		Gopt_vs_E2z3->Draw("same l");	
		Gopt_vs_E3z3->Draw("same l");
	/////////////////////////////////////////////////////////////////////////////////////////
}
/////////////////////////////////////////////////////////////////////////////////////
