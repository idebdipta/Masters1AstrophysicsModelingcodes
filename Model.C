#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <gsl/gsl_sf_dilog.h>

#include "TMatrixD.h"
#include "TMath.h"

//Compute the luminosity distance///////////////////////////////////////////////////////////////////////////////////////////
double lumi_distance_Kernel(double z, double omegaM, double omegaLambda){
	return 1./sqrt(omegaLambda+omegaM*(1.+z)*(1.+z)*(1.+z));
}

std::vector<double> lumi_LoadDistances(std::vector< double > vz,	double omegaM, double omegaLambda){// [z]
	std::vector<double> vLum;
	double step_z = (vz.back()-vz.front())/vz.size();

	//Initialization for z=0
	vLum.push_back(0.);

	//Loop on z
	for(int i=1; i<(int)vz.size(); i++){
		double kz_z =  lumi_distance_Kernel(vz[i],omegaM,omegaLambda);
		double kz_zm1 =  lumi_distance_Kernel(vz[i-1],omegaM,omegaLambda); 
		vLum.push_back(vLum[i-1]+0.5*step_z*(kz_z+kz_zm1));//trapezoidal rule
	}

	return vLum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Compute the EBL optical depth/////////////////////////////////////////////////////////////////////////////////////////////
//Cosmology kernel//////////////////////////////////////////////////////////////
double eblmodel_Redshift_Kernel(double z, double omegaM, double omegaLambda, double f_evol){
	return 1./(pow(1.+z,2.+f_evol)*sqrt(omegaLambda+omegaM*(1.+z)*(1.+z)*(1.+z)));
}

//Particle-physics kernel///////////////////////////////////////////////////////
double eblmodel_ParticlePhysics_Kernel(double x){	
	double ln1px = log(1.+x);
	double ln1mx = log(1.-x);
	double ln2 = log(2);

	double res = ln2*ln2 -c_pi*c_pi/6. + 2.*gsl_sf_dilog(0.5*(1.-x));
	res+= (ln1px-2.*ln2)*ln1mx + 0.5*(ln1mx*ln1mx - ln1px*ln1px);
	res+= (0.5*(ln1px-ln1mx)*(1.+x*x*x*x)-(x+x*x*x))/(1.-x*x);

	return  res;
}

//Particle-physics cosmology kernel x e0^3///////////////////////////////////////
std::vector< std::vector<double> > eblmodel_e0_cube_times_Kernel(std::vector< double > vz,	std::vector< double > ve,	double omegaM, double omegaLambda, double f_evol){// [z][E]
	std::vector< std::vector<double> > vKernel;
	double step_z = (vz.back()-vz.front())/vz.size();

	//Initialization for z=0
	vKernel.push_back( std::vector<double>() );
	for(int j=0; j<(int)ve.size(); j++) vKernel[0].push_back(0.);

	//Loop on z
	for(int i=1; i<(int)vz.size(); i++){
		vKernel.push_back( std::vector<double>() );
		double kz_z =  eblmodel_Redshift_Kernel(vz[i],omegaM,omegaLambda,f_evol);
		double kz_zm1 =  eblmodel_Redshift_Kernel(vz[i-1],omegaM,omegaLambda,f_evol); 
		for(int j=0; j<(int)ve.size(); j++){
			
			double ke_z =  0.; 
			double ke_zm1 = 0.; 
			if(((1.+vz[i-1])>exp(-0.5*ve[j]))){
				ke_z =  eblmodel_ParticlePhysics_Kernel( sqrt( 1.-exp(-ve[j])/((1.+vz[i])*(1.+vz[i])) ) );
				ke_zm1 = eblmodel_ParticlePhysics_Kernel( sqrt( 1.-exp(-ve[j])/((1.+vz[i-1.])*(1.+vz[i-1.])) ) ); 
			}
			vKernel[i].push_back( vKernel[i-1][j] + 0.5*step_z*(kz_zm1*ke_zm1+kz_z*ke_z) );
		}
	}

	return vKernel;
}

//Initialization of the Kernel and of the useful vectors/////////////////////////
void Model_InitializeKernel(double omegaM, double omegaLambda, double f_evol, double zmax, std::vector<double>& veGamma, std::vector<double>& vEgamma, std::vector<double>& ve, std::vector<double>& vz, std::vector< std::vector<double> >& vKernel){
	//Initialization of the vectors used in the EBL integration///////////////
	double Emin_simu = 0.03;//TeV
	double Emax_simu = 25.;//TeV

	//veGamma[i] = log(Ei/mc^2) where Ei is the gamma-ray energy
	double eGamma_min = log(Emin_simu/0.511E-6);// 50 GeV / 511 keV
	double eGamma_max = log(Emax_simu/0.511E-6);// 50 TeV / 511 keV
	double deGamma = 0.1;// Typical energy resolution
	veGamma.clear();
	veGamma = utils_InitializeVector(eGamma_min, eGamma_max, deGamma);
	
	//Energy
	vEgamma.clear();
	for(int i=0; i<(int)veGamma.size(); i++) vEgamma.push_back( exp(veGamma[i])*c_m_e_eV*c_eV_2_TeV );

	//vz
	double dz = 0.001;
	vz.clear();
	vz = utils_InitializeVector(0.,zmax,dz);

	//ve[i] = log(epsilonEBL*E/m2c4) with energies as of today
	double emin = -2*log(1.+zmax);
	double emax = log(10000.);// integration up to exp(emax) times the threshold energy
	double de = 0.01;// de should be smaller energy resolution in Egamma or in Eebl
	ve.clear();
	ve = utils_InitializeVector(emin,emax, de); 

	//kernel
	vKernel.clear();
	vKernel = eblmodel_e0_cube_times_Kernel(vz, ve, omegaM, omegaLambda, f_evol);
}

//Reduced optical depth of a single Gaussian in the N-Gaussian model
double Model_ComputeSingleOpacity(double ei, double dln_lambda, double eGamma, std::vector<double> vKernel_fixed_z, std::vector< double > ve, double H0){
	double integral = 0;

	double sigma_e = dln_lambda/(2.*sqrt(2.*log(2.)));
	double RH = c_speed/(H0*c_km_2_cm/c_Mpc_2_cm);//cm
	double normalization = 3.*c_pi*(c_sigmaT/c_speed)*RH*exp(eGamma)/(c_m2_2_cm2*utils_eV2nJ(c_m_e_eV));
	double step_e = (ve.back()-ve.front())/ve.size();
	for(int j=0; j<(int)ve.size(); j++){//Integration over log(e0)
		double gaussian = exp(-0.5*(eGamma + ei - ve[j])*(eGamma + ei - ve[j])/(sigma_e*sigma_e));
		double integrand = normalization*gaussian*exp(-3*ve[j])*vKernel_fixed_z[j]*step_e;
		if((j==0) || (j==(int)ve.size()-1)) integrand *= 0.5;//Trapezoidal integration with a uniform grid in ve
		integral += integrand;
	}

	return integral;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Load the EBL SED//////////////////////////////////////////////////////////////////////////////////////////////////////////
//model with a sum of logparabolas
double ebl_model_N_logPar(double *x, double *p){
	double res = 0;

	int nbump = (int)p[0];
	int npar_lp = 3;
	for(int i=0; i<nbump; i++){
		double log_x = log(x[0]) - p[2+i*npar_lp];
		res+= p[1+i*npar_lp]*exp(-log_x*log_x*p[3+i*npar_lp]);
	}

	return res;
}

//relative uncertainty at a given flux for the sum of logparabolas
double ebl_model_relative_unc_N_logpar(TF1 *febl_model, double lambda){
	double res = 1.0;

	int nbump = (int)febl_model->GetParameter(0);
	int npar_lp = 3;
	
	double phi=0;
	double sigma_phi2=0;
	for(int i=0; i<nbump; i++){
		double a = febl_model->GetParameter(1+i*npar_lp);
		double ea = febl_model->GetParError(1+i*npar_lp);
		double b = febl_model->GetParameter(2+i*npar_lp);
		double eb = febl_model->GetParError(2+i*npar_lp);
		double c = febl_model->GetParameter(3+i*npar_lp);
		double ec = febl_model->GetParError(3+i*npar_lp);
	
		double log_x = log(lambda) - b;
		double phi0 = a*exp(-log_x*log_x*c);

		double da = phi0/a;
		double db = -2.*log_x*c*phi0;
		double dc = -log_x*log_x*phi0;

		phi += phi0;	
		sigma_phi2 += da*da*ea*ea + db*db*eb*eb + dc*dc*ec*ec;
	}

	if(phi!=0) res = sqrt(sigma_phi2)/abs(phi);
	
	return res;
}

//Loads the model with a sum of logparabolas
TF1 *Model_LoadMultiBumpAndError(std::string filename, double lambda_min, double lambda_max, bool verbose=false){
	int npar_lp = 3;
	std::vector< std::vector< double > > vpar;
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
			if((int)vbuf.size()==2*npar_lp) vpar.push_back(vbuf);
			else if(verbose) std::cout<<"Line "<<line<<" unread in EBL file "<<filename<<std::endl;
		}
		myfile.close();
	}
	else std::cout<<"Could not open file: "<<filename<<std::endl;

	int nbump = (int)vpar.size();
	int npar = 1 + nbump*npar_lp;
	TF1 *fintensity =  new TF1("flogN",ebl_model_N_logPar,lambda_min,lambda_max,npar);
	fintensity->SetParameter(0,nbump);
	for(int i=0; i<nbump; i++){
		for(int j=0; j<npar_lp; j++){
			fintensity->SetParameter(1+j+npar_lp*i,vpar[i][j]);
			fintensity->SetParError(1+j+npar_lp*i,vpar[i][3+j]);
		}
	}

	return fintensity;	
}

//Sets the Gaussian weights of the EBL model //0: n, 1:dln_lambda, [2;n+1]:ai, [n+2;2n+1]:xi
std::vector< std::vector< double > > Model_InitializeGaussianWeights(double H0, double lambda_min, double lambda_max, double dln_lambda, TF1* febl_model = new TF1()){// [i][0]: value - [i][1]: error

	double H0_ref = 67.;

	std::vector< std::vector< double > > vRes;
	vRes.push_back(std::vector< double >());//[0] -> param
	vRes.push_back(std::vector< double >());//[1] -> uncertainty

	int npts = round(log(lambda_max/lambda_min)/dln_lambda);
	double ln_lambda_start = log(lambda_max) - (npts - 0.5)*dln_lambda;

	vRes[0].push_back(npts); vRes[1].push_back(0.1);
	vRes[0].push_back(dln_lambda); vRes[1].push_back(0.1);

	//first push-back 0 for the a_i
	for(int i=0; i<npts; i++){
		vRes[0].push_back(0.);
		vRes[1].push_back(0.1);
	}

	//then sets the xi=ln(lambda_i)	
	for(int i=0; i<npts; i++){
		vRes[0].push_back(ln_lambda_start+i*dln_lambda);
		vRes[1].push_back(0.1);
	}

	if(febl_model->GetNpar()>0){
		//fills the matrix that takes into account the leakage from one EBL bin to the other
		TMatrixD mat_i0i(npts,npts);
		double sl2 = sqrt(log(2));
		double const_norm = 0.25*sqrt(c_pi)/sl2;
		for(int i0=0; i0<npts; i0++){
			for(int i=0; i<npts; i++) mat_i0i[i0][i] = const_norm*(TMath::Erf(sl2*(2*(i0-i)+1)) - TMath::Erf(sl2*(2*(i0-i)-1)));
		}
		//Invert the matrix to fill the weigths based on the average flux in -dx/2 +dx/2
		TMatrixD mat_ii0 = mat_i0i.Invert();
		for(int i=0; i<npts; i++){
			double ai = 0.;
			for(int i0=0; i0<npts; i0++){
				ai += mat_ii0[i][i0]*febl_model->Eval(exp(vRes[0][i0+npts+2]));
	 		} 
			ai*=H0/H0_ref;//So that it shifts up and down the EBL
			vRes[0][i+2]=ai;
			vRes[1][i+2]=abs(ai)*ebl_model_relative_unc_N_logpar(febl_model,exp(vRes[0][i+npts+2]));//Crude estimate of the uncertainty based on the uncertainty on the flux at this wavelength
		}
	}
	else for(int i=0; i<npts; i++) vRes[0][i+2]=1.;

	return vRes;
}

//Loads the EBL model based on a sum of Gaussians//par: 0: n, 1:dln_lambda, [2;n+1]:ai, [n+2;2n+1]:xi
double ebl_model_MultiGauss(double *x, double *p){
	double res = 0;

	int npts = p[0];
	double dln_lambda=p[1];
	double sigma = dln_lambda/(2.*sqrt(2.*log(2.)));//so that gauss(x0 + dx/2) + gauss(x1 - dx/2) = gaus(x0) =  gaus(x1);
	for(int i=2; i<npts+2; i++)	res+=p[i]*TMath::Gaus(log(x[0]),p[i+npts],sigma);

	return res;
}

//Loads the EBL model based on a sum of Gaussians//par: 0: n, 1:dln_lambda, [2;n+1]:ai, [n+2;2n+1]:xi
TF1 *Model_fMultiGauss(double* par, double lambda_min, double lambda_max){
	
	//Initialization of the function
	int npts = par[0];
	int npar = 2*npts+2;

	TF1 *fintensity =  new TF1("flogN",ebl_model_MultiGauss,lambda_min,lambda_max,npar);
	fintensity->FixParameter(0,npts);
	fintensity->FixParameter(1,par[1]);
	for(int i=2; i<npts+2; i++){
		fintensity->SetParameter(i,par[i]);//ai
		fintensity->FixParameter(i+npts,par[i+npts]);//xi
	}

	return fintensity;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
