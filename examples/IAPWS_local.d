/**
 * steam.d
 * 
 * Contents:
 * 	 1.IAPWS formulations for state calculation.
 * 	 	1.1 Formulations to calculate saturation (p,T) & boundary equation
 *			for region 2 and 3.
 *		1.2 IAPWS-Region2 formulation struct.
 *		1.3 IAPWS base class whose object has all the calculated values
 *			of state properties.
 * Author: Jiasheng(Jason) Liang
 * Version:
 */

import std.stdio;
import std.math;
import std.string;
import std.conv;


//steam constants:
immutable double R=461.526; /// specific gas constant[J/kg/K]
immutable double rho_c=322; /// critical density [kg/m^^3]	
immutable double T_c=647.096; /// critical thermal temperature [K]
immutable double p_c=22.064e6; /// critical pressure [Pa]
immutable double k_c=1e-3; /// critical thermal conductivity [W/K/m]
immutable double mu_c=1e-6; // critical dynamic viscosity [Pa.s]

//-------------------------------------------------------------------------------
//PART1.1 Formulations to calculate saturation (p,T) & boundary equation
//						for region 2 and 3.
//-------------------------------------------------------------------------------

//table for calculating saturation pressure and Temperature 
double[10] table_2_19=
    [0.11670521452767e4,-0.72421316703206e6,-0.17073846940092e2, 0.1202082470247e5,
     -0.32325550322333e7,0.1491510861353e2,-0.48232657361591e4,0.40511340542057e6,
     -0.23855557567849,0.65017534844798e3];

double get_pb23(double T){
    /*
     *	calculate boundary value of pressre with a given temperature using eqn 2.1 
     *	and the coefficients inf table 2.1.
     */

    return 1e6*(0.34805185628969e3-0.11671859879975e1*T +0.10192970039326e-2*T^^2);
} 
double get_ps(double T){
    /*
     *	calculate saturation pressure with a given temperature using eqn 2.12b,
     * 	eqn 2.13 and the coefficients in table 2.19.
     */
    //initialise intermediate properties 
    double vartheta = T+table_2_19[8]/(T-table_2_19[9]);
    double A = vartheta^^2+table_2_19[0]*vartheta+table_2_19[1];
    double B = table_2_19[2]*vartheta^^2+table_2_19[3]*vartheta
	+table_2_19[4];
    double C = table_2_19[5]*vartheta^^2+table_2_19[6]*vartheta
	+table_2_19[7];
    return 1e6*(2*C/(-1*B+(B^^2-4*A*C)^^0.5))^^4;
}
double get_Ts(double p){
    /*
     *	calculate saturation Temperature with a given pressure using eqn 2.12a, 
     *	eqn 2.14 and the coefficients in table 2.19
     */	
    double _beta = (p/1e6)^^0.25;
    double E = _beta^^2+table_2_19[2]*_beta+table_2_19[5];
    double F = table_2_19[0]*_beta^^2+table_2_19[3]*_beta+table_2_19[6];
    double G = table_2_19[1]*_beta^^2+table_2_19[4]*_beta+table_2_19[7];
    double D = 2*G/(-1*F-(F^^2-4*E*G)^^0.5);
    return (table_2_19[9]+D-((table_2_19[9]+D)^^2-4*(table_2_19[8] +table_2_19[9]*D))^^0.5)/2;
}

//---------------------------------------------------------------------------------
//PART 1.2. IAPWS-Region2 formulation struct
//---------------------------------------------------------------------------------	
struct Region2{
private:
    //state properties
    double p,T;	//pressure & temperature
    double quality=1.0; //default value for vapor quality 
    //intermediate properties
    double g,gamma_o,gamma_r,pi,tau; 
    //derivatives of intermediate properties
    ///refer to table 2.9 for more details 
    double gamma_o_tau,gamma_o_pi,gamma_o_tautau,gamma_o_pipi,gamma_o_pitau;
    ///refer to table 2.10 for more details
    double gamma_r_tau,gamma_r_pi,gamma_r_tautau,gamma_r_pipi,gamma_r_pitau;

    //Implementing tables 
    double[2][9] table_2_6=
	[/*[J_i, n_i]*/
	 [0,-0.96927686500217e1], [1,0.10086655968018e2], [-5,-0.56087911283020e-2],
	 [-4,0.71452738081455e-1], [-3,-0.40710498223928],
	 /*i=5*/
	 [-2,0.14240819171444e1], [-1,-0.43839511319450e1], [2,-0.28408632460772],
	 [3,0.21268463753307e-1] ];
    double[3][43] table_2_7=
	[/*[I_i,J_i,n_i]*/
	 [1,0,-0.0017731742473213], [1,1,-0.017834862292358], [1,2,-0.045996013696365],
	 [1,3,-0.057581259083432], [1,6,-0.05032527872793], [2,1,-3.3032641670203e-05],
	 [2,2,-0.00018948987516315], [2,4,-0.0039392777243355], [2,7,-0.043797295650573],
	 [2,36,-2.6674547914087e-05], [3,0,2.0481737692309e-08], [3,1,4.3870667284435e-07],
	 [3,3,-3.227767723857e-05], [3,6,-0.0015033924542148], [3,35,-0.040668253562649],
	 /*16*/
	 [4,1,-7.8847309559367e-10], [4,2,1.2790717852285e-08], [4,3,4.8225372718507e-07],
	 [5,7,2.2922076337661e-06], [6,3,-1.6714766451061e-11], [6,16,-0.0021171472321355],
	 [6,35,-23.895741934104], [7,0,-5.905956432427e-18], [7,11,-1.2621808899101e-06],
	 [7,25,-0.038946842435739], [8,8,1.1256211360459e-11], [8,36,-8.2311340897998],
	 [9,13,1.9809712802088e-08], [10,4,1.0406965210174e-19], [10,10,-1.0234747095929e-13],
	 /*31*/
	 [10,14,-1.0018179379511e-09], [16,29,-8.0882908646985e-11], [16,50,0.10693031879409],
	 [18,57,-0.33662250574171], [20,20,8.9185845355421e-25], [20,35,3.0629316876232e-13],
	 [20,48,-4.2002467698208e-06], [21,21,-5.9056029685639e-26], [22,53,3.7826947613457e-06],
	 [23,39,-1.2768608934681e-15], [24,26,7.3087610595061e-29], [24,40,5.5414715350778e-17],
	 [24,58,-9.436970724121e-07] ];
    double[3][13] table_2_12=
	[/*[I_i,J_i,n_i]*/
	 [1,0,-0.73362260186506e-02], [1,2,-0.88223831943146e-01], [1,5,-0.72334555213245e-01],
	 [1,11,-0.40813178534455e-02], [2,1,0.20097803380207e-02], [2,7,-0.53045921898642e-01],
	 [2,16,-0.76190409086970e-02], [3,4,-0.63498037657313e-02], [3,16,-0.86043093028588e-01],
	 [4,7,0.75321581522770e-02], [4,10,-0.79238375446139e-02], [5,9,-0.22888160778447e-03],
	 [5,10,-0.26456501482810e-02] ];
	
    this(double _p, double _T, double _quality){
	p = _p;
	T = _T;
	quality = _quality;
	init;
    }

    void init(){
	/*
	 * contains all properties and coefficients required to proceed calculation in Region 2.
	 */
	pi=p*1e-6;
	tau=540/T;

	//pure steam
	if(quality==1.0)
	    {
		//table 2.9, derivatives of gamma_o
		///container for summing value during calculation
		double sum=0;
		double sum_tau=0;
		double sum_tautau=0;
		// [TODO] Jason, I prefer the foreach loop these days.
		// foreach (i; 0 .. 9) { ... }
		for (int i=0; i<9; i++){
		    sum += table_2_6[i][1]*(tau^^table_2_6[i][0]);
		    sum_tau += table_2_6[i][1]*table_2_6[i][0]*(tau)^^(table_2_6[i][0]-1);
		    sum_tautau += table_2_6[i][1]*table_2_6[i][0]*(table_2_6[i][0]-1)
			*(tau)^^(table_2_6[i][0]-2);
		}
		gamma_o = log(pi) + sum;
		gamma_o_tau = sum_tau;
		gamma_o_tautau = sum_tautau;
		gamma_o_pi = 1/pi; 
		gamma_o_pipi = -(pi)^^(-2); 
		gamma_o_pitau = 0;

		//table 2.10, derivatives of gamma_r
		sum=0;
		sum_tau=0;
		sum_tautau=0;
		double sum_pi=0;
		double sum_pipi=0;
		double sum_pitau=0;
		for (int i=0; i<43; i++){
		    sum += table_2_7[i][2]*((pi)^^(table_2_7[i][0]))
			*((tau-0.5)^^(table_2_7[i][1]));
		    sum_tau += table_2_7[i][2]*((pi)^^(table_2_7[i][0]))
			*((tau-0.5)^^(table_2_7[i][1]-1))*table_2_7[i][1];
		    sum_pi += table_2_7[i][2]*table_2_7[i][0]*(pi^^(table_2_7[i][0]-1))
			*((tau-0.5)^^(table_2_7[i][1]));
		    sum_tautau += table_2_7[i][2]*((pi)^^(table_2_7[i][0]))
			*((tau-0.5)^^(table_2_7[i][1]-2))*table_2_7[i][1]*(table_2_7[i][1]-1);
		    sum_pipi += table_2_7[i][2]*table_2_7[i][0]*(table_2_7[i][0]-1)
			*(pi^^(table_2_7[i][0]-2))*((tau-0.5)^^(table_2_7[i][1])); 
		    sum_pitau += table_2_7[i][2]*table_2_7[i][0]*(pi^^(table_2_7[i][0]-1))
			*table_2_7[i][1]*((tau-0.5)^^(table_2_7[i][1]-1));		
		}
		gamma_r = sum;
		gamma_r_tau = sum_tau;
		gamma_r_pi = sum_pi; 
		gamma_r_tautau = sum_tautau;
		gamma_r_pipi = sum_pipi;
		gamma_r_pitau = sum_pitau;
	    }

	//metastable-vapour region		
	if((quality<1.0)&&(quality>=0.95))
	    {
		//according to eqn 2.9, the gamma_o and its derivatives 
		//within supplementary equation for metastable-vapour region
		double sum=0;
		double sum_tau=0;
		double sum_tautau=0;
		table_2_6[0][1] = -0.96937268393049e1;
		table_2_6[1][1] = 0.10087275970006e2;
		for (int i=0; i<9; i++){
		    sum += table_2_6[i][1]*(tau^^table_2_6[i][0]);
		    sum_tau += table_2_6[i][1]*table_2_6[i][0]*(tau)^^(table_2_6[i][0]-1);
		    sum_tautau += table_2_6[i][1]*table_2_6[i][0]*(table_2_6[i][0]-1)
			*(tau)^^(table_2_6[i][0]-2);
		}
		gamma_o = log(pi) + sum;
		gamma_o_tau = sum_tau;
		gamma_o_tautau = sum_tautau;
		gamma_o_pi = 1/pi;
		gamma_o_pipi = -(pi)^^(-2);
		gamma_o_pitau = 0;

		//table 2.13, the gamma_r and its derivaties within supplementary 
		//equation for metastable-vapour region
		sum=0;
		sum_tau=0;
		sum_tautau=0;
		double sum_pi=0;
		double sum_pipi=0;
		double sum_pitau=0;
		for (int i=0; i<13; i++){
		    sum += table_2_12[i][2]*((pi)^^(table_2_12[i][0]))
			*((tau-0.5)^^(table_2_12[i][1]));
		    sum_tau += table_2_12[i][2]*((pi)^^(table_2_12[i][0]))
			*((tau-0.5)^^(table_2_12[i][1]-1))*table_2_12[i][1];
		    sum_pi += table_2_12[i][2]*table_2_12[i][0]*(pi^^(table_2_12[i][0]-1))
			*((tau-0.5)^^(table_2_12[i][1]));
		    sum_tautau += table_2_12[i][2]*((pi)^^(table_2_12[i][0]))
			*((tau-0.5)^^(table_2_12[i][1]-2))*table_2_12[i][1]*(table_2_12[i][1]-1);
		    sum_pipi += table_2_12[i][2]*table_2_12[i][0]*(table_2_12[i][0]-1)
			*(pi^^(table_2_12[i][0]-2))*((tau-0.5)^^(table_2_12[i][1])); 
		    sum_pitau += table_2_12[i][2]*table_2_12[i][0]*(pi^^(table_2_12[i][0]-1))
			*table_2_12[i][1]*((tau-0.5)^^(table_2_12[i][1]-1));		
		}
		gamma_r = sum;
		gamma_r_tau = sum_tau;
		gamma_r_pi = sum_pi;
		gamma_r_tautau = sum_tautau;
		gamma_r_pipi = sum_pipi;
		gamma_r_pitau = sum_pitau;
	    }

	//eqn 2.6 & eqn 2.9
	g = (gamma_r + gamma_o)*R*T;
    }	

public:
    //table 2.8, relations of thermodynamic properties to gamma_o and gamma_r 
    double SpecificVolume(){
	return (pi*(gamma_o_pi+gamma_r_pi))*R*T/p;
    } 
    double SpecificEnthalpy(){
	return tau*(gamma_o_tau + gamma_r_tau)*R*T;
    }  
    double SpecificInternalEnergy(){
	return SpecificEnthalpy - SpecificVolume*p;
    } 
    double SpecificEntropy(){
	return SpecificEnthalpy/T - (gamma_o+gamma_r)*R;
    }  
    double SpecificIsobaricHeatCapacity(){
	return -(tau^^2)*(gamma_o_tautau+gamma_r_tautau)*R;
    } 
    double SpecificIsochoricHeatCapacity(){
	return SpecificIsobaricHeatCapacity
	    -(((1+pi*gamma_r_pi-tau*pi*gamma_r_pitau)^^2)/(1-pi^^2*gamma_r_pipi))*R;
    }
    double SoundSpeed(){
	double _a,_b,_c,_d;
	_a = 1+2*pi*gamma_r_pi+pi^^2*gamma_r_pi^^2;
	_b = 1-pi^^2*gamma_r_pipi;
	_c = (1+pi*gamma_r_pi-tau*pi*gamma_r_pitau)^^2;
	_d = tau^^2*(gamma_o_tautau+gamma_r_tautau);
	return sqrt((_a/(_b+(_c/_d)))*R*T);
    } 	
    double IsobaricCubicExpansionCoefficient(){
	return (1+pi*gamma_r_pi - tau*pi*gamma_r_pitau)
	    /(1+pi*gamma_r_pi)/T;
    }
    double IsothermalCompressibility(){
	return (1-pi^^2*gamma_r_pipi)/(1+pi*gamma_r_pi)/p;
    }
} // end struct Region2


//-------------------------------------------------------------------------------
// PART 1.2. IAPWS base formulation class
//-------------------------------------------------------------------------------

class IAPWS{
    /*
     *	Unless specify, otherwise:
     *	Reference:
     *	 Wanger, W., & Kretzschmar, H.(2008). International Steam Tables.
     *	 Berlin, Heidelberg: Springer Berlin Heidelberg.
     */
private:
    // a char representing region of state within IAPWS
    string region; // e.g. '2': Formulation Region 2 of IAPWS
    // distance from the input p-T to the closest saturation point 
    double dis_satoff;
    //tolerance for the distance to the closest saturation point 
    immutable SATURATION_TOL = 0.01;

    //tables for calculating dynamic viscosity:
    double[4] table_3_1=[0.167752e-1,0.220462e-1,0.6366564e-2,-0.241605e-2];
    double[3][21] table_3_2=
	[/*[I_i,J_i,n_i]*/
	 [0,0,0.520094],[0,1,0.850895e-1],[0,2,-0.108374e1],[0,3,-0.289555],[1,0,0.222531],
	 /*5*/ [1,1,0.999115],[1,2,1.88797],[1,3,0.126613e1],[1,5,0.120573],[2,0,-0.281378],
	 [2,1,-0.906851],[2,2,-0.772479],[2,3,-0.489837],[2,4,-0.25704],[3,0,0.161913],
	 [3,1,0.257399],[4,0,-0.325372e-1],[4,3,0.698452e-1],[5,4,0.872102e-2],
	 [6,3,-0.435673e-2],[6,5,-0.593264e-3] ];
	
    //table for calculating thermal conductivity:
    double[5] table_1=[2.443221E-03,1.323095000E-02,6.770357000E-03,-3.454586000E-03,4.096266000E-04];
    double[6][5] table_2=
	[[1.60397357,-0.646013523,0.111443906,0.102997357,-0.0504123634,0.00609859258],
	 [2.33771842,-2.78843778,1.53616167,-0.463045512,0.0832827019,-0.00719201245],
	 [2.19650529,-4.54580785,3.55777244,-1.40944978,0.275418278,-0.0205938816],
	 [-1.21051378,1.60812989,-0.621178141,0.0716373224,0,0],
	 [-2.7203370,4.57586331,-3.18369245,1.1168348,-0.19268305,0.012913842]];
    double[5][6] table_6=
	[[6.53786807199516,6.52717759281799,5.35500529896124,1.55225959906681,1.11999926419994],
	 [-5.61149954923348,-6.30816983387575,-3.96415689925446,0.464621290821181,0.595748562571649],
	 [3.39624167361325,8.08379285492595,8.91990208918795,8.93237374861479,9.8895256507892],
	 [-2.27492629730878,-9.82240510197603,-12.033872950579,-11.0321960061126,-10.325505114704],
	 [10.2631854662709,12.1358413791395,9.19494865194302,6.1678099993336,4.66861294457414],
	 [1.97815050331519,-5.54349664571295,-2.16866274479712,-0.965458722086812,-0.503243546373828]];

    //determine which region the input (p,T) lies in
    void set_region(){

	//(p,T) lies in region 2
	if(((0<T)&&(T<623.15))&&(p<=get_ps(T)) ||
	   ((623.15<=T&&T<=863.15) &&(get_pb23(T)>p))||((863.15<=T)&&(T<=1073.15))) 
	{
		    region = "2";
	} 
	else
	{
		writeln("input state not in IAPWS-Region2");
		writeln("input p:",p, ", input T:", T, ", saturation p:", get_ps(T));
	}
    }

    void update_thermo(){
	set_region;

	//initialise an object of the corresponding region class 
	//and proceed the calculation within that class
	switch (region) {
	case "2":
	    Region2 _IAPWS = Region2(p,T,quality);
	    u = _IAPWS.SpecificInternalEnergy;
	    h = _IAPWS.SpecificEnthalpy;
	    s = _IAPWS.SpecificEntropy;
	    v = _IAPWS.SpecificVolume;
	    rho = 1/v;
	    Cp = _IAPWS.SpecificIsobaricHeatCapacity; 
	    Cv = _IAPWS.SpecificIsochoricHeatCapacity; 
	    a = _IAPWS.SoundSpeed; 
	    mu = DynamicViscosity;
	    mu = ThermalConductivity;
	    alpha_v = _IAPWS.IsobaricCubicExpansionCoefficient;
	    kappa_T = _IAPWS.IsothermalCompressibility;
	    break;	 

	default:
	    string msg;
	    msg ~= format("Warning in function: %s:\n", __FUNCTION__);
	    msg ~= format("    Input state is out of the valid range of IAPWS formulations of state.\n"); 
	    writeln(msg);
	    break;

	}//end switch 
    }// end update_thermo

public:
	//Thermodynamic properties:
    double T; /// thermal temperature [K]
    double p; /// pressure [Pa]
    double quality; /// vapour quality (1: pure gas) [-] 
    double u; /// specific internal energy [J/kg]
    double h; /// specfic enthalpy [J/kg]
    double s; /// specific entropy [J/K/kg]
    double v; /// specific volume [m^^3/kg]
    double rho; /// density [kg/m^^3] 
    double Cp; /// specific isobaric heat capacity [J/K/kg]
    double Cv; /// specific isochoric heat capacity [J/kg/K]
    double a; /// sound speed [m/s]
    double mu; /// dynamic viscosity [Pa.s]
    double k; /// thermal conductivity [W/m/K]
    double alpha_v; /// isobaric cubic expansion coefficient [1/K]
    double kappa_T; /// isothermal compressibility [1/Pa]
    this(double _p, double _T, double _quality){
	T = _T;p = _p; quality = _quality;
	update_thermo;
    }

    //function to compute dynamic viscosity but not in IAPWS-IF97
    //parameter: rho, T
    //valid in: 273.15K <= T <= 1173.15K and p <= 100 MPa
    double DynamicViscosity(){	
	//intermediate properties
	double delta=rho/rho_c;
	double theta=T/T_c;
	double psi_0,psi_1;
	double sum=0;
		
	//eqn 3.2
	for(int i=0;i<4;++i){
	    sum += table_3_1[i]*theta^^(-i);
	}
	psi_0 = theta^^0.5*(sum)^^-1;

	//eqn 3.3
	sum=0;
	for(int i=0;i<21;++i){
	    sum += table_3_2[i][2]*(delta-1)^^(table_3_2[i][0])
		*((theta^^-1-1)^^(table_3_2[i][1])); 
	}
	psi_1 = exp(delta*sum);

	//eqn 3.1
	return 1e-6*psi_0*psi_1;
    }
	
    double ThermalConductivity(){
	/*	
	 * contains everything implementing from IAPWS R15-11 for industrial use
	 * reference:
	 * 	IAPWS (2011). Release on the IAPWS Formulaiton 2011 for the Thermal 
	 *	Conductivity of Ordinary Water Substance,
	 * 	available at the IAPWS website http://www.iapws.org 
	 * 
	 * parameters: p, T, rho, mu
	 * valid in: 273.15K <= T <= 1173.15K and p <= 100 MPa;
	 */
	
	//intermediate properties
	double lambda_bar,lambda_0,lambda_1,lambda_2,dvdp_T,drhodp_T,zeta_T,
	    zeta_TR,xi,delta_x, y, Z;
	///eqn 7 ~ eqn 13
	double T_bar = T/T_c;
	double p_bar = p/p_c;
	double rho_bar = rho/rho_c;
	double mu_bar = mu/mu_c;
	double Cp_bar = Cp/R;
	double kappa = Cp/Cv; 
	///dummy sum container	
	double sum=0;
	double sum_1=0;

	//eqn 16
	for(int i=0;i<5;++i){sum += table_1[i]/(T_bar^^i);} 
	lambda_0 = sqrt(T_bar)/sum;
	//eqn 17
	sum=0;
	for(int i=0;i<5;++i){
	    for(int j=0;j<6;++j){sum_1 += table_2[i][j]*(rho_bar-1)^^j;}
	    sum += ((T_bar^^-1-1)^^i)*sum_1;sum_1=0;}
	lambda_1 = exp(rho_bar*sum);	

	/*
	 *lambda_2, critical enhancement of thermal conductivity
	 *	refer to http://www.twt.mpei.ac.ru/mcs/worksheets/iapws/wspTCPT.xmcd	
	 *	for more details
	 */
	if((Cp_bar<0)||(Cp_bar>1e13)){
	    Cp_bar = 1e13;Cp = Cp_bar*R;kappa = Cp/Cv;
	}	
	dvdp_T = -kappa_T*v;
	drhodp_T = rho^^2*-dvdp_T; 
	zeta_T = drhodp_T*p_c/rho_c;
	///eqn 26
	int _j=-1;
	do{
	    if(rho_bar<=0.310559006){_j=0;}
	    if((0.310559006<rho_bar)&&(rho_bar<=0.776397516)){_j=1;}
	    if((0.776397516<rho_bar)&&(rho_bar<=1.242236025)){_j=2;}
	    if((1.242236025<rho_bar)&&(rho_bar<=1.863354037)){_j=3;}
	    if(rho_bar>1.863354037){_j=4;}
	}while(_j==-1);
	sum=0;for(int i=0;i<6;++i){sum += table_6[i][_j]*rho_bar^^i;}
	zeta_TR = 1/sum; 
	///eqn 23
	delta_x = rho_bar*(zeta_T-zeta_TR*1.5/T_bar); 
	if(delta_x<0){delta_x=0;}
	///eqn 22
	xi = 0.13*(delta_x/0.06)^^(0.63/1.239); 
	if((xi<0)||(xi>1e13)){xi=1e13;}
	///eqn 20
	y = xi/0.4;
	///eqn 19
	if(y<1.2e-7){Z = 0;}
	else{Z = (2/PI/y)*((1-1/kappa)*atan(y)+y/kappa
			   -(1-exp(-1/(1/y+y^^2/3/(rho_bar^^2)))));}
	///eqn 18 
	lambda_2 = 177.8514*rho_bar*Cp_bar*T_bar*Z/mu_bar;
		
	//eqn 15
	lambda_bar = lambda_0*lambda_1+lambda_2;
	//eqn 10
	return lambda_bar * k_c;
    } // end ThermalConductivity   
} // end class IAPWS

