--GasModel Steam: Eilmer4 Input file for the Uncertainty Analysis 
--				   of updateThermoFromRHOU

gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}
Q.quality = 1 --> assume always in gas phase

file = io.open("updateThermoFromRHOU.txt","w")
p_data = io.open("p.txt","r")
T_data = io.open("T.txt","r")

-- lists of input (p,T) points
p_list = {}; T_list = {}  
for _p in p_data:lines() do
	table.insert(p_list,tonumber(_p))
end

for _T in T_data:lines() do
	table.insert(T_list,tonumber(_T))
end

-- lists of properties calculated by updateThermoFromPT
rho_list={}; u_list={}; h_list={}; s_list={}; Cp_list={} 
Cv_list={}; a_list={}; mu_list={}; k_list = {}
--lists of properties calculated by updateThermoFromRHOU
--based on the (rho, u) calculated previously
p_list_1={}; T_list_1={}; h_list_1={}; s_list_1={}; Cp_list_1={} 
Cv_list_1={}; a_list_1={}; mu_list_1={}; k_list_1 = {}

for i_T in pairs(T_list) do
	Q.T = T_list[i_T]
	Q.p = p_list[i_T]
	--calculate all the thermo properties 
	gmodel:updateThermoFromPT(Q)
	gmodel:updateSoundSpeed(Q)
	gmodel:updateTransCoeffs(Q)
	table.insert(rho_list, Q.rho)
	table.insert(u_list, Q.u)
	table.insert(h_list, gmodel:enthalpy(Q))
	table.insert(s_list, gmodel:entropy(Q))
	table.insert(Cv_list, gmodel:Cv(Q))
	table.insert(Cp_list, gmodel:Cp(Q))
	table.insert(a_list, Q.a)
	table.insert(mu_list, Q.mu)
	table.insert(k_list, Q.k)
end

--lists of properties calculated by updateThermoFromRHOU
--based on the (rho, u) calculated previously
p_list_1={}; T_list_1={}; h_list_1={}; s_list_1={}; Cp_list_1={} 
Cv_list_1={}; a_list_1={}; mu_list_1={}; k_list_1 = {}
for i_rho in pairs(rho_list) do
	Q.rho = rho_list[i_rho]
	Q.u = u_list[i_rho]
	gmodel:updateThermoFromRHOU(Q)
	gmodel:updateSoundSpeed(Q)
	gmodel:updateTransCoeffs(Q)
	table.insert(p_list_1, Q.rho)
	table.insert(T_list_1, Q.T)
	table.insert(h_list_1, gmodel:enthalpy(Q))
	table.insert(s_list_1, gmodel:entropy(Q))
	table.insert(Cv_list_1, gmodel:Cv(Q))
	table.insert(Cp_list_1, gmodel:Cp(Q))
	table.insert(a_list_1, Q.a)
	table.insert(mu_list_1, Q.mu)
	table.insert(k_list_1, Q.k)
end

print(rho_list[1000]," ", u_list[1000]," ",p_list_1[1000]," ", T_list_1[1000]," ", h_list_1[1000]," ", s_list_1[1000]," ", Cp_list_1[1000]," ", Cv_list_1[1000]," ", a_list_1[1000]," ", mu_list_1[1000]," ", k_list_1[1000])


Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--1. updateThermoFromPT:")
gmodel:updateThermoFromPT(Q)
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Density [kg/m3]=", Q.rho, "Interal Energy [J/kg]=", Q.u, "Sound Speed [m/s]=",
	 Q.a,"\n")
file:write("Density [kg/m3]=")
file:close()





