/**
 * sat_p(T).d
 * 
 * Local IAPWS functions to generate a set of (p,T) points in pT.txt
 * for the use of analysing uncertainties of steam.getpT_from_rhou
 */

import std.stdio, std.file, std.math;

//table for calculating saturation pressure and Temperature 
immutable double[10] table_2_19=
    [0.11670521452767e4,-0.72421316703206e6,-0.17073846940092e2, 0.1202082470247e5,
     -0.32325550322333e7,0.1491510861353e2,-0.48232657361591e4,0.40511340542057e6,
     -0.23855557567849,0.65017534844798e3];

double get_pb23(double T){
    /*
     *  calculate boundary value of pressre with a given temperature using eqn 2.1 
     *  and the coefficients inf table 2.1.
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

void main()
{
	File p_list = File("p.txt","w");
    File T_list = File("T.txt","w");
    //data points under saturation boundary
    for(double T = 300; T<623.15;T+=10){
        for(double p = 100; p<get_ps(T);p+=0.05*get_ps(T)){
            p_list.writeln(p);
            T_list.writeln(T);
        }
    } 
    //gas data points in IAPWS-region2 and T < 863.15
    for(double T = std.math.round(623.15)+1; T<863.15;T+=10){
        for(double p = 100; p<get_pb23(T); p+=0.05*get_pb23(T)){
            p_list.writeln(p);
            T_list.writeln(T);
        }
    }
    //gas data points in IAPWS-region2 and T > 863.15
    for(double T = std.math.round(863.15)+1; T<1073.15;T+=10){
        for(double p = 100; p<=100e6; p+=5e6){
            p_list.writeln(p);
            T_list.writeln(T);
        }
    }
    p_list.close(); T_list.close();
}