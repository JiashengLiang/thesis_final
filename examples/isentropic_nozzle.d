/**
 * isentropic_nozzle.d
 * 
 * Module to calculate the exit velocity of an circular nozzle under isentropic and 
 * adiabatic operation. Refer to the final report for more details.
 * 
 * Contents:
 * 	 1. A(r) & r(x)
 *   2. ideal()
 *	 3.	iapws() 
 *	 4. main()   
 * Author: Jiasheng(Jason) Liang
 * Version:
 */
import std.stdio;
import std.math;
import IAPWS_local;

//-------------------------------------------------------------------------------
//PART 1. Functions to return the inner surface radius and cross-section area 
//		  at x [mm] for the nozzle 
//-------------------------------------------------------------------------------

//radius [mm]
double r(double x){
	assert(x>=0);
	double r; //radius in [mm]
	if(x<=5)
	{
		return 5;
	}
	else if(5<x && x<=6.65)
	{
		r = sqrt(abs(1.65^^2-(x-5)^^2))+3.35;
		return r;
	}
	else if(6.65<x && x<=8.45)
	{
		r = -sqrt(abs(1.8^^2-(x-8.45)^^2))+3.35;
		return r;
	}
	else
	{
		r = tan(PI*5/180)*(x-8.45)+1.55;
		return r;
	}
}

//Area [m^2]
double A(double x){ return PI * (1e-3*r(x))^^2;}

//-------------------------------------------------------------------------------
//PART 2. Function using ideal-gas laws and isentropic processing relations to 
//		  return a set of updated variables after a step in horizontal
//		  location (delta_x [mm]) is taken 
//-------------------------------------------------------------------------------
double[][] ideal(double[] phy_var, double[] thermo_var, double delta_x){
	//constants
	assert(IAPWS_local.R); 						/// specific gas constant[J/kg/K]
	double kappa = 1.308;     					/// isentropic exponent kappa 

	//initial physical variables 
	double x_old = phy_var[0];					/// horizontal location [mm]
	double A_old = phy_var[1];					/// cross-section area [m^2]
	double V_old = phy_var[2];					/// fluid velocity [m/s]
	//initial thermaldynamic variables
	double p_old = thermo_var[0];				/// pressure [Pa]
	double T_old = thermo_var[1];				/// temperature [K]

	//Intermedia thermo variables
	double rho = p_old/R/T_old;					/// Density [kg/m^3]
	double M = V_old/sqrt(kappa*R*T_old);		/// Mach number

	//update variables
	//--x
	double x_new = x_old + delta_x;
	//--Area
	double A_new = A(x_old+delta_x); 
	double delta_A = A_new- A_old;
	//--Velocity
	double delta_V = -V_old/(1-M^^2)*delta_A/A_old;
	double V_new = V_old + delta_V; 
	//--temperature
	double delta_T = (1-kappa)*(M^^2)*T_old*delta_V/V_old;
	double T_new = T_old + delta_T;
	//--pressure
	double p_new = p_old * (T_new/T_old)^^(kappa/(kappa-1));

	//update variable lists
	phy_var = [x_new, A_new, V_new]; thermo_var = [p_new, T_new];

	return [phy_var,thermo_var]; 
}

//-------------------------------------------------------------------------------
//PART 3. Function using local IAPWS database and isentropic processing relations 
//		  to return a set of updated variables after a step in horizontal
//		  location (delta_x [mm]) is taken 
//-------------------------------------------------------------------------------
double[][] iapws(double[] phy_var, double[] thermo_var, double delta_x){
	//constants
	assert(IAPWS_local.R); 						/// specific gas constant[J/kg/K]

	//initial physical variables 
	double x_old = phy_var[0];					/// horizontal location [mm]
	double A_old = phy_var[1];					/// cross-section area [m^2]
	double V_old = phy_var[2];					/// fluid velocity [m/s]
	//initial thermaldynamic variables
	double p_old = thermo_var[0];				/// pressure [Pa]
	double T_old = thermo_var[1];				/// temperature [K]

	//local IAPWS equations of state
	IAPWS iapws_old = new IAPWS(p_old, T_old, 1);/// an object of class IAPWS
	double rho = iapws_old.rho;					/// Density [kg/m^3]
	double M = V_old/iapws_old.a;				/// Mach number	
	double s_old = iapws_old.s;					/// Specific entropy [J/K/kg]

	//update variables
	//--x
	double x_new = x_old + delta_x;
	//--Area
	double A_new = A(x_old+delta_x); 
	double delta_A = A_new- A_old;
	//--Velocity
	double delta_V = -V_old/(1-M^^2)*delta_A/A_old;
	double V_new = V_old + delta_V;
	//--Pressure
	double delta_p = -rho*V_old*delta_V;
	double p_new = p_old + delta_p;
	//--temperature
	double Newton(){
		/*
		* Use 1D Newton's Method to find the T_new such that:
		* 	s(p_new, T_new) - s_old = 0	= f(p_new, T_new)  
		*/
		double f(double p, double T){
			IAPWS guess = new IAPWS(p, T, 1); 		///assume always be gas
			return guess.s-s_old;
		}
		double dfdT(double p, double T, double _h){
			/* use 2-point central finite difference to differentiate f(p,T)*/
			return (f(p, T+_h)-f(p, T-_h))/(2*_h);
		}
		//initial guess 
		double delta_T = 1; 
		double temp = T_old + delta_T; //local variable for temperature [K]
		double f_new = f(p_new, temp);
		//take the percentage uncertainty of isobaric enthalpy for IAPWS-region2 as the iteration tolerance
		//Reference: fig.2.32, Wanger, W., & Kretzschmar, H.(2008). 
		//           International Steam Tables. Berlin, Heidelberg: 
		//			 Springer Berlin Heidelberg.
		double tol;
		if(p_new<=1e6){tol=0.002*s_old;}
		else{tol = 0.003*s_old;}
		while(f_new>=tol){
			delta_T = -f_new/dfdT(p_new,temp,delta_T);
			//prepare for next iteration
			temp += delta_T;
			f_new = f(p_new,temp);
		}
		return temp;
	}
	double T_new = Newton();
	double delta_T = T_new - T_old;

	//update variable lists
	phy_var = [x_new, A_new, V_new]; thermo_var = [p_new, T_new];

	return [phy_var,thermo_var]; 
}

//-------------------------------------------------------------------------------
//PART 4. main() to execute functions and generate results 
//-------------------------------------------------------------------------------
void main(){
	double x_delta; 
	
	writeln("importing IAPWS_local.d...");	
	assert(IAPWS_local.R==461.526);
	writeln("local IAPWS database is imported...");
	write("step size of x in mm:");
	readf("%f", &x_delta);
	writeln(10*x_delta);
}