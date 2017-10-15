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
	if(x<=30){return (23./6000.)*x^^2+(-23./100.)*x+5;}
	else{return (3./2023.)*x^^2+(-180./2023.)*x+116713./40460.;}
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
	assert(_IAPWS.R); 							/// specific gas constant[J/kg/K]
	//double R = 287.058;						    /// R = 287.058 [J/kg/K] for air;
	double R = _IAPWS.R; 								/// R = _IAPWS.R for steam
	double kappa = 1.308;     					/// isentropic exponent kappa for steam
												/// 1.4 for air, 1.308 for steam 

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
double s_old = 7026.906765;					/// Specific entropy at inlet [J/kg K]

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
		//---local functions to assist the iteration 
	double f(double p, double T){
		_IAPWS.IAPWS guess = new _IAPWS.IAPWS(p, T, 1); ///assume always be gas
		return (guess.s - s_old);
	}
	double dfdT(double p, double T){
		/* use 2-point central finite difference to differentiate f(p,T)*/
		double _h = 1e-5;
		if(f(p, T-_h) != f(p, T-_h)) // when the left boundary has acrossed the saturation boundary
		{
			_h = (T-_IAPWS.get_Ts(p))/2;
		}
		return (f(p, T+_h)-f(p, T-_h))/(2*_h);
		
	}
	double Newton(){
		/*
		* Use 1D Newton Method to find the T_new such that:
		* 	s(p_new, T_new) - s_old = 0	= f(p_new, T_new)  
		*/
		//update Temperature
		double T_0 = T_old; //T remains the same if there is no Area change 
							//and f_new = 0
		//initial guess
		double h = 5e-5; 
		double f_old = f(p_new,T_0);
		int i=0;			  //iteration number
		while((f_old!=f_old || abs(f_old)>0.002*s_old) && i<30)
		{  //entropy tolerance .2% of initial entropy and maximum iterations  
			// check if the updated temperature has acrossed the saturation boundary
			if(f_old != f_old) // execute when f_new = nan
			{
				writeln("adjust T_new from ",T_0," [K]");
				T_0 = _IAPWS.get_Ts(p_new) + (_IAPWS.get_Ts(p_new)-T_0)/2;
				writeln("to ", T_0," [K]");
				f_old = f(p_new,T_0);
			}
			//calculate suitable step size for this iteration 
			double dfdt = (f(p_new,T_0+h)-f_old)/h;
			h = -(f_old-0.002*s_old)/dfdt;

			//prepare for next iteration
			T_0 += h;
			f_old = f(p_new,T_0);
			i+=1;

/**/		//if iapws is not enough to find T_0 when entropy is fixed, let it equal to something slightly bigger than T_sat
			if(f_old!=f_old && i==30)
			{
				writeln("Exceed IAPWS-Region 2 valid range, let update temperature be",
						 "saturation temperature at update pressure.");
				T_0 = _IAPWS.get_Ts(p_new)+5e-5;
				f_old = f(p_new,T_0);
			}
			logs.writeln("x:",x_new,", V:", V_new, ", Pressure:", p_new,", iteration:",i,", T:", T_0, ", f_old:", f_old,", dfdT:",dfdt,", h:", h);
		}
		return T_0;
	}//end Newton()
	double T_new = Newton(); 

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
	
	assert(_IAPWS.R==461.526);
	writeln("local IAPWS database is imported...");

	//variables at the nozzle inlet
	//-- [physical: (x[mm], area[m^2], velocity [m/s]),
	//-- thermal: (pressure[Pa], temperature [K]) ]
	double[][] init_var = [[0,A(0),30.62983],[270e3,403.15]];
	double[][] ideal_var = init_var;
	double[][] iapws_var = init_var;

	writeln("Inlet conditions:", init_var);

	double delta_x;
	write("Step size of x in mm:");readf("%f", &delta_x);
	
	//stepping along x-axis
	while(ideal_var[0][0]<89.5)
	{
		iapws_var = iapws(iapws_var[0],iapws_var[1], delta_x);
		iapws_data.writeln(iapws_var[0][0]," ", iapws_var[0][1], " ",iapws_var[0][2],
							" ",iapws_var[1][0]," ",iapws_var[1][1]);
		ideal_var = ideal(ideal_var[0],ideal_var[1], delta_x);
		ideal_data.writeln(ideal_var[0][0]," ",ideal_var[0][1], " ",ideal_var[0][2],
							" ", ideal_var[1][0]," ",ideal_var[1][1]);
		writeln("Processing...up to x = ", ideal_var[0][0]);
	}
	
	writeln("outlet conditions:ideal:", ideal_var, ", iapws:",iapws_var);
	//close txt files
	iapws_data.close(); ideal_data.close();
}