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
import std.file;
static import _IAPWS = IAPWS_local;


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
	else if(5<x && x<6.65)
	{
		r = sqrt(abs(1.65^^2-(x-5)^^2))+3.35;
		return r;
	}
	else if(6.65<x && x<8.45)
	{
		r = -sqrt(abs(1.8^^2-(x-8.45)^^2))+3.35;
		return r;
	}
	else if(x==6.65) //infinite gradient
	{
		return 3.35;
	}
	else if(x==8.45) //zero gradient
	{
		return 1.55;
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
	assert(_IAPWS.R); 						/// specific gas constant[J/kg/K]
	double kappa = 1.308;     					/// isentropic exponent kappa 

	//initial physical variables 
	double x_old = phy_var[0];					/// horizontal location [mm]
	double A_old = phy_var[1];					/// cross-section area [m^2]
	double V_old = phy_var[2];					/// fluid velocity [m/s]
	//initial thermaldynamic variables
	double p_old = thermo_var[0];				/// pressure [Pa]
	double T_old = thermo_var[1];				/// temperature [K]

	//Intermedia thermo variables
	double rho = p_old/_IAPWS.R/T_old;					/// Density [kg/m^3]
	double M = V_old/sqrt(kappa*_IAPWS.R*T_old);		/// Mach number

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
	assert(_IAPWS.R); 						/// specific gas constant[J/kg/K]
	File logs=File("log.txt", "a"); 
	//initial physical variables 
	double x_old = phy_var[0];					/// horizontal location [mm]
	double A_old = phy_var[1];					/// cross-section area [m^2]
	double V_old = phy_var[2];					/// fluid velocity [m/s]
	//initial thermaldynamic variables
	double p_old = thermo_var[0];				/// pressure [Pa]
	double T_old = thermo_var[1];				/// temperature [K]

	//local IAPWS equations of state
	_IAPWS.IAPWS iapws_old = new _IAPWS.IAPWS(p_old, T_old, 1);/// an object of class IAPWS
	double rho = iapws_old.rho;					/// Density [kg/m^3]
	double M = V_old/iapws_old.a;				/// Mach number	
	double h_old = 2.7201154e6;					/// Specific enthalpy at inlet [J/kg]

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
	double Bisection(){
		/*
		* Use 1D Bisection Method to find the T_new such that:
		* 	h(p_new, T_new) - h_old = 0	= f(p_new, T_new)  
		*/
		double f(double p, double T){
			_IAPWS.IAPWS guess = new _IAPWS.IAPWS(p, T, 1); 		///assume always be gas
			return guess.h-h_old;
		}
		double dfdT(double p, double T){
			/* use 2-point central finite difference to differentiate f(p,T)*/
			double _h = 1e-3;
			return (f(p, T+_h)-f(p, T-_h))/(2*_h);
		}
		
		//update Temperature
		double T_0 = T_old; //T remains the same if there is no Area change 
							  //and f_new = 0
		//initial guess 
		double h = 5e-4; //step in temperature
		double f_new = f(p_new, T_0);
		int i=0;			  //iteration number
		while(abs(f_new)>=1e3 && i<30){ //entropy tolerance 10 [J/kg K]
										 //maximum iterations  
			//calculate suitable step size for each iteration 
			if(dfdT(p_new,T_0)>0)
			{
				if(f(p_new,T_0)<0)
				{
					if(f(p_new,T_0+h)<0) {h*=2;}
					else				 {h/=2;}
				}
				else if(f(p_new,T_0)>0)
				{
					h*=-1;
					if(f(p_new,T_0+h)>0) {h*=2;}
					else				 {h/=2;}
				}
			}
			else if(dfdT(p_new,T_0)<0)
			{
				if(f(p_new,T_0)<0)
				{
					h*=-1;
					if(f(p_new,T_0+h)>0) {h/=2;}
					else				 {h*=2;}
				}
				else if(f(p_new,T_0)>0)
				{
					if(f(p_new,T_0+h)<0) {h/=2;}
					else				 {h*=2;}
				}
			}

			//prepare for next iteration
			T_0+=h;
			//--if the updated T has acrossed the saturation boundary
			if(T_0 < _IAPWS.get_Ts(p_new))
			{
				T_0 = _IAPWS.get_Ts(p_new);
			}
			f_new = f(p_new,T_0);
			i+=1;
			logs.writeln("x:",x_new,", V:", V_new, ", Pressure:", p_new,", iteration:",i,", T:", T_0, ", f_new:", f_new,", dfdT:",dfdT(p_new,T_0));
			writeln("processing... at x=", x_new);
		}
		logs.close();
		return T_0;
	}
	double T_new = Bisection();

	//update variable lists
	phy_var = [x_new, A_new, V_new]; thermo_var = [p_new, T_new];

	return [phy_var,thermo_var]; 
}

//-------------------------------------------------------------------------------
//PART 4. main() to execute functions and generate results 
//-------------------------------------------------------------------------------
void main(){
	//construct the txt files at which the simulation data will be store
	File ideal_data = File("isentropic_nozzle_ideal.txt", "w");
	ideal_data.writeln("[x[mm], Area[m], Velocity[m/s]],[Pressure[Pa], Temperature[K]]");
	File iapws_data = File("isentropic_nozzle_iapws.txt", "w");
	iapws_data.writeln("[x[mm], Area[m], Velocity[m/s]],[Pressure[Pa], Temperature[K]]");
	
	double delta_x; //step size [mm]
	
	assert(_IAPWS.R==461.526);
	writeln("local IAPWS database is imported...");
	
	//inputs
	write("step size of x in mm:");readf("%f", &delta_x);

	//variables at the nozzle inlet
	//-- [physical: (x[mm], area[m^2], velocity [m/s]),
	//-- thermal: (pressure[Pa], temperature [K]) ]
	double[][] ideal_var = [[5, A(5), 30.62983],[270e3, 403.15]];
	double[][] iapws_var = ideal_var;
	double[] ideal_phy,ideal_thermo,iapws_phy,iapws_thermo;
	writeln("Inlet conditions:", ideal_var);
	
	//stepping along x-axis
	while(ideal_var[0][0]<67.95)
	{
		iapws_phy = iapws_var[0]; iapws_thermo = iapws_var[1];
		iapws_var = iapws(iapws_phy,iapws_thermo, delta_x);
		iapws_data.writeln(iapws_var[0][0]," ", iapws_var[0][1], " ",iapws_var[0][2],
							" ",iapws_var[1][0]," ",iapws_var[1][1]);
		ideal_phy = ideal_var[0]; ideal_thermo = ideal_var[1];
		ideal_var = ideal(ideal_phy,ideal_thermo, delta_x);
		ideal_data.writeln(ideal_var[0][0]," ", ideal_var[0][1], " ",ideal_var[0][2],
							" ", ideal_var[1][0]," ",ideal_var[1][1]);
	} 
	//close txt files
	/*iapws_data.close(); ideal_data.close();
	_IAPWS.IAPWS guess = new _IAPWS.IAPWS(270e3, 403.15, 1);
	writef("%.8f",guess.h);*/
}