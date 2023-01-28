#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>



using namespace std;

using dvector = vector<double>;

void initialize_phi(dvector &phi, dvector &rphi, dvector &dz_rphi);
void initialize_coeff(dvector &coeff);
void set_value(dvector &vec, const double &value);
void num_integrate_z(dvector &output, const dvector &integrand);
void deriv_z(dvector &output, const dvector &input);
void deriv_s(dvector &output, const dvector &vec1, const dvector &vec2);
void print_vec(const dvector &vec);
void time_step_rphi(dvector &dz_rphi_next, const dvector &dz_rphi, dvector &ddz_rphi,dvector &rphi_next,  const dvector &rphi, const dvector &coeff,
		const dvector &z_array, const double &curr_s);
void phi_from_rphi(dvector &phi, const dvector &rphi, const dvector &dz_rphi, const double &curr_s);
void iteration_loop(dvector &phi, dvector &rphi, dvector &dz_rphi, const dvector &coeff, const dvector &z_array);

void initialize_CSV(const dvector &phi, const dvector &dz_rphi, const dvector &ddz_rphi, const dvector &rphi, const dvector &dz_phi, const dvector &ds_rphi);
void update_CSV(const dvector &phi, const dvector &dz_rphi, const dvector &ddz_rphi, const dvector &rphi, const dvector &dz_phi,const dvector &ds_rphi);
void max_value(double &max, dvector &vec);


double ds = .0001; //time step
double dz = .0005; //spatial step
double total_s =10; double total_z = 1; //total time and space unit
const int ns = total_s/ds + 1; const int nz = total_z/dz +1;
int total_recorded = 500; int recordinterval = ns / total_recorded; //Will record 500 data points
double damped_constant = .00001; // optional damping

double k = .1; //wave parameter



int main (){
	dvector phi(nz,0); dvector rphi(nz,0); dvector dz_rphi(nz,0); dvector coeff(nz,0); dvector zeros(nz,0);
	dvector ddz_rphi(nz,0);

	dvector z_array(nz,0);
	for (int i =0; i <= nz-1 ; i++){
		z_array[i] = -1*total_z + i * dz;
	};

	initialize_phi(phi,rphi,dz_rphi); initialize_coeff(coeff);
	deriv_z(ddz_rphi, dz_rphi);

	initialize_CSV(phi, dz_rphi, ddz_rphi, rphi, zeros, zeros);

	iteration_loop(phi, rphi, dz_rphi, coeff, z_array);


	return 0;
};

void iteration_loop(dvector &phi, dvector &rphi, dvector &dz_rphi, const dvector &coeff, const dvector &z_array){
	dvector phi_next(nz,0); dvector rphi_next(nz,0);  dvector dz_rphi_next(nz,0); dvector dz_phi(nz,0);
	dvector ddz_rphi(nz,0); dvector ds_rphi(nz,0); dvector ds_phi(nz,0);
  double curr_s = 0;

	int counter = 0; double max =0;

  for (int i =0; i <= ns; i++){

		if (counter % recordinterval == 0){ // we'll record first
			update_CSV(phi, dz_rphi, ddz_rphi, rphi, dz_phi,ds_rphi);
			max_value(max, dz_rphi);
			cout << "current value of s is " << curr_s << "\n";

			if(max < damped_constant){
				cout << "solution has dispersed by time "<< curr_s << "\n";
				break;
			};
			if (dz_rphi[nz-1] > 200){
  			cout << "simulation broke at s value " << curr_s << "\n";
  			break;
  		};
		};

      time_step_rphi(dz_rphi_next, dz_rphi, ddz_rphi, rphi_next, rphi, coeff, z_array,curr_s);

  		num_integrate_z(rphi_next, dz_rphi_next);
  		phi_from_rphi(phi_next, rphi_next, dz_rphi_next, curr_s);
  		deriv_s(ds_rphi, rphi_next, rphi);
  		deriv_s(ds_phi,phi_next,phi);


  		dz_rphi = dz_rphi_next; rphi = rphi_next; phi = phi_next;
  		deriv_z(dz_phi,phi);

  		counter = counter + 1;
      curr_s = curr_s + ds;

  };
};

void time_step_rphi(dvector &dz_rphi_next, const dvector &dz_rphi, dvector &ddz_rphi,dvector &rphi_next,  const dvector &rphi, const dvector &coeff, const dvector &z_array, const double &curr_s){
	deriv_z(ddz_rphi, dz_rphi); double temp = 0;
	double growth = 0;
	double nonlin = 0;

	double oscillations = 0;
	double coeff_dec = 0;
	double slow_timescale=k; // if you want to kill this effect, change to =k


	for (int i=1; i <= nz-2; i++){
		temp =  -1*(1-k)*(slow_timescale*1/k)*z_array[i]*ddz_rphi[i] + exp(-curr_s*coeff_dec)*(1+oscillations * cos(curr_s))*(coeff[i]+nonlin*pow(dz_rphi[i],2))*rphi[i] + growth*dz_rphi[i] ;
		dz_rphi_next[i] = .5*(dz_rphi[i-1] + dz_rphi[i+1]) + ds * temp;
	};
  dz_rphi_next[0] = dz_rphi_next[1];
	dz_rphi_next[nz-1] = dz_rphi[nz-1] + ds * exp(-curr_s*coeff_dec)*(1+oscillations * cos(curr_s))*(coeff[nz-1]+nonlin*pow(dz_rphi[nz-1],2))*rphi[nz-1] + ds*dz_rphi[nz-1]*growth;
};

void initialize_phi(dvector &phi, dvector &rphi, dvector &dz_rphi){
	set_value(phi,0); set_value(rphi,0); double abs_z = 0;
	for (int i =0; i <= nz - 1; i++){
		abs_z = total_z - i*dz;
		dz_rphi[i] = -1*(.3*pow(1+abs_z,.5)+.5*cos(200*abs_z)*exp(-400*pow(abs_z,2)) + .5*cos(111*abs_z)*exp(-400*pow(abs_z-1,2)) - .1*cos(3.1*abs_z)*cos(3.2*abs_z)+.6*sin(6*abs_z));
		//sample oscillatory initial data
	};
	num_integrate_z(rphi,dz_rphi);
};

void initialize_coeff(dvector &coeff){
	double abs_z = 0;
	for (int i =0; i <= nz - 1; i++){
			abs_z = total_z - i*dz;
			coef[i] = -k; //sample constant potential
	};
};

void set_value(dvector &vec, const double &value){
	int len = vec.size();
	for (int i =0; i <= len-1; i++){
		vec[i] = value;
	}
}

void num_integrate_z(dvector &output, const dvector &integrand){ // uses trapezoid rule, assumes 1) integrating from z=-1, 2) zero boundary conditions
	set_value(output,0); double temp = 0;
	for (int i =1; i<= nz-1; i++){
		temp = .5 * dz * (integrand[i-1] + integrand[i]);
		output[i] = output[i-1] + temp;
	};
};

void deriv_z(dvector &output, const dvector &input){
	for (int i =1; i <= nz - 2; i++){
		output[i] = 1/(2*dz)*(input[i+1] - input[i-1]);
	};
	output[0] = 1/(2*dz)*(-3*input[0] + 4*input[1] - input[2]);
	output[nz-1] = 1/(2*dz)*(3*input[nz-1] - 4*input[nz-2] + input[nz-3]);
};

void deriv_s(dvector &output, const dvector &vec1, const dvector &vec2){
	for (int i =0; i <= nz -1; i++){
		output[i] = 1/ds * (vec2[i] - vec1[i]);
	}
};

void print_vec(const dvector &vec){
	int len = vec.size();
	for (int i =0; i<= len-1; i++){
		cout << vec[i] << ", " << "\n";
	}
}

void phi_from_rphi(dvector &phi, const dvector &rphi, const dvector &dz_rphi,const double &curr_s){
	for (int i =1; i<= nz - 1; i++){
		phi[i]= rphi[i] * 2 * exp(k*curr_s) * 1/(dz*i);
	};
	phi[0] = 2 * dz_rphi[0] * exp(k*curr_s);
};

void initialize_CSV(const dvector &phi, const dvector &dz_rphi, const dvector &ddz_rphi, const dvector &rphi, const dvector &dz_phi, const dvector &ds_rphi){
	ofstream file("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_phi.csv");
	      for(int i = 0; i <= nz-1; i++){
	          file<< phi[i] <<",";
	      };
	      file << "\n";
	      file.close();

	ofstream file2("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_dz_rphi.csv");
	      	      for(int i = 0; i <= nz-1; i++){
	      	          file2<< dz_rphi[i] <<",";
	      	      };
	      	      file2 << "\n";
	      	      file2.close();
ofstream file2_1("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_ddz_rphi.csv");
								      	      for(int i = 0; i <= nz-1; i++){
								      	          file2_1<< ddz_rphi[i] <<",";
								      	      };
								      	      file2_1 << "\n";
								      	      file2_1.close();
	ofstream file3("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_dz_phi.csv");
	      	      for(int i = 0; i <= nz-1; i++){
	      	    	  file3<< dz_phi[i] <<",";
	      	      };
	      	      file3 << "\n";
	      	      file3.close();
    ofstream file4("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_ds_rphi.csv");
	      	      for(int i = 0; i <= nz-1; i++){
	      	    	   file4<< ds_rphi[i] <<",";
	      	      };
	      	      file4 << "\n";
	      	     file4.close();
	ofstream file5("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_rphi.csv");
						 							for(int i = 0; i <= nz-1; i++){
						 								 file5<< rphi[i] <<",";
						 							};
						 							file5 << "\n";
						 						 file5.close();
}

void update_CSV(const dvector &phi, const dvector &dz_rphi, const dvector &ddz_rphi, const dvector &rphi, const dvector &dz_phi, const dvector &ds_rphi){
	ofstream file("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_phi.csv", ios_base::app);
		  for(int i = 0; i <= nz-1; i++){
		       file<< phi[i] <<",";
		   };
		   file << "\n";
		   file.close();

	ofstream file2("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_dz_rphi.csv", ios_base::app);
		   		  for(int i = 0; i <= nz-1; i++){
		   		       file2<< dz_rphi[i] <<",";
		   		   };
		   		   file2 << "\n";
		   		   file2.close();

		ofstream file2_1("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_ddz_rphi.csv", ios_base::app);
					 		   		  for(int i = 0; i <= nz-1; i++){
					 		   		       file2_1<< ddz_rphi[i] <<",";
					 		   		   };
					 		   		   file2_1 << "\n";
					 		   		   file2_1.close();
    ofstream file3("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_dz_phi.csv", ios_base::app);
		  for(int i = 0; i <= nz-1; i++){
		       file3<< dz_phi[i] <<",";
		   };
		   file3 << "\n";
		   file3.close();
    ofstream file4("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_ds_rphi.csv", ios_base::app);
		   for(int i = 0; i <= nz-1; i++){
		   	    file4<< ds_rphi[i] <<",";
		    };
		    file4 << "\n";
		    file4.close();
		ofstream file5("/Users/Singh/eclipse-workspace/ss_wave_new/simulation_rphi.csv", ios_base::app);
					 for(int i = 0; i <= nz-1; i++){
								file5<< rphi[i] <<",";
						};
						file5 << "\n";
						file5.close();
}

void max_value(double &max, dvector &vec){
	int len = vec.size(); double temp = 0;
	for (int i =0; i<= len-1; i++){
		 if(abs(vec[i]) > temp) {
			 temp = abs(vec[i]);
		 }
	}
	max = temp;
	cout << max << "\n";
};
