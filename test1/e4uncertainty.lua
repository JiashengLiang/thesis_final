--GasModel Steam: Eilmer4 Input file for the Uncertainty Analysis 
--				   of updateThermoFromRHOU

gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}
Q.quality = 1 --> assume always in gas phase

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
p_data:close(); T_data:close()

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
--based on the (rho, u) calculated from updateThermoFromPT
p_list_1={}; T_list_1={}; h_list_1={}; s_list_1={}; Cp_list_1={} 
Cv_list_1={}; a_list_1={}; mu_list_1={}; k_list_1 = {}
for i_rho in pairs(rho_list) do
	Q.rho = rho_list[i_rho]
	Q.u = u_list[i_rho]
	gmodel:updateThermoFromRHOU(Q)
	gmodel:updateSoundSpeed(Q)
	gmodel:updateTransCoeffs(Q)
	table.insert(p_list_1, Q.p)
	table.insert(T_list_1, Q.T)
	table.insert(h_list_1, gmodel:enthalpy(Q))
	table.insert(s_list_1, gmodel:entropy(Q))
	table.insert(Cv_list_1, gmodel:Cv(Q))
	table.insert(Cp_list_1, gmodel:Cp(Q))
	table.insert(a_list_1, Q.a)
	table.insert(mu_list_1, Q.mu)
	table.insert(k_list_1, Q.k)
end

--lists of properties calculated by updateThermoFromPS
--based on the s calculated from updateThermoFromPT
--and with the initial temperature is 20 K larger or smaller 
--than the actual
T_list_2={}; h_list_2={}; Cp_list_2={}; Cv_list_2={}; a_list_2={}
mu_list_2={}; k_list_2 = {}
for i_s in pairs(s_list) do
	s = s_list[i_s]
	Q.p = p_list[i_s]
	T = T_list[i_s]+20
	if T > 1073.15 then T = T_list[i_s]-20 end
	Q.T = T
	gmodel:updateThermoFromPS(Q,s)
	gmodel:updateSoundSpeed(Q)
	gmodel:updateTransCoeffs(Q)
	table.insert(T_list_2, Q.T)
	table.insert(h_list_2, gmodel:enthalpy(Q))
	table.insert(Cv_list_2, gmodel:Cv(Q))
	table.insert(Cp_list_2, gmodel:Cp(Q))
	table.insert(a_list_2, Q.a)
	table.insert(mu_list_2, Q.mu)
	table.insert(k_list_2, Q.k)
end


--write the data to txt file
p_txt = io.open("uncertainty_p.txt","w")
p_txt:write("updateThermoFromPT ", " updateThermoFromRHOU",
	"  difference","  percentage","\n")
for i in pairs(p_list) do
	p_txt:write(p_list[i],"	 ",p_list_1[i],"  ",
	       math.abs(p_list[i]-p_list_1[i]),"  ",
	       100*math.abs((p_list[i]-p_list_1[i])/p_list[i]),"\n")
end
p_txt:close()

T_txt = io.open("uncertainty_T.txt","w")
T_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(T_list) do
	T_txt:write(T_list[i],"	 ",T_list_1[i],"  ",
	       math.abs(T_list[i]-T_list_1[i]),"  ",
	       100*math.abs((T_list[i]-T_list_1[i])/T_list[i]),"  ",
	       T_list_2[i],"  ",math.abs(T_list[i]-T_list_2[i]),"  ",
	       100*math.abs((T_list[i]-T_list_2[i])/T_list[i]),"\n")
end
T_txt:close()

h_txt = io.open("uncertainty_h.txt","w")
h_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(h_list) do
	h_txt:write(h_list[i],"	 ",h_list_1[i],"  ",
	       math.abs(h_list[i]-h_list_1[i]),"  ",
	       100*math.abs((h_list[i]-h_list_1[i])/h_list[i]),"  ",
	       h_list_2[i],"  ",math.abs(h_list[i]-h_list_2[i]),"  ",
	       100*math.abs((h_list[i]-h_list_2[i])/h_list[i]),"\n")
end
h_txt:close()

s_txt = io.open("uncertainty_s.txt","w")
s_txt:write("updateThermoFromPT ", " updateThermoFromRHOU",
	"  difference","  percentage","\n")
for i in pairs(s_list) do
	s_txt:write(s_list[i],"	 ",s_list_1[i],"  ",
	       math.abs(s_list[i]-s_list_1[i]),"  ",
	       100*math.abs((s_list[i]-s_list_1[i])/s_list[i]),"\n")
end
s_txt:close()

Cp_txt = io.open("uncertainty_Cp.txt","w")
Cp_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(Cp_list) do
	Cp_txt:write(Cp_list[i],"  ",Cp_list_1[i],"  ",
	       math.abs(Cp_list[i]-Cp_list_1[i]),"  ",
	       100*math.abs((Cp_list[i]-Cp_list_1[i])/Cp_list[i]),"  ",
	       Cp_list_2[i],"  ",math.abs(Cp_list[i]-Cp_list_2[i]),"  ",
	       100*math.abs((Cp_list[i]-Cp_list_2[i])/Cp_list[i]),"\n")
end
Cp_txt:close()

Cv_txt = io.open("uncertainty_Cv.txt","w")
Cv_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(Cv_list) do
	Cv_txt:write(Cv_list[i],"  ",Cv_list_1[i],"  ",
	       math.abs(Cv_list[i]-Cv_list_1[i]),"  ",
	       100*math.abs((Cv_list[i]-Cv_list_1[i])/Cv_list[i]),"  ",
	       Cv_list_2[i],"  ",math.abs(Cv_list[i]-Cv_list_2[i]),"  ",
	       100*math.abs((Cv_list[i]-Cv_list_2[i])/Cv_list[i]),"\n")
end
Cv_txt:close()

a_txt = io.open("uncertainty_a.txt","w")
a_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(a_list) do
	a_txt:write(a_list[i],"	 ",a_list_1[i],"  ",
	       math.abs(a_list[i]-a_list_1[i]),"  ",
	       100*math.abs((a_list[i]-a_list_1[i])/a_list[i]),"  ",
	       a_list_2[i],"  ",math.abs(a_list[i]-a_list_2[i]),"  ",
	       100*math.abs((a_list[i]-a_list_2[i])/a_list[i]),"\n")
end
a_txt:close()

mu_txt = io.open("uncertainty_mu.txt","w")
mu_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(mu_list) do
	mu_txt:write(mu_list[i],"  ",mu_list_1[i],"  ",
	       math.abs(mu_list[i]-mu_list_1[i]),"  ",
	       100*math.abs((mu_list[i]-mu_list_1[i])/mu_list[i]),"  ",
	       mu_list_2[i],"  ",math.abs(mu_list[i]-mu_list_2[i]),"  ",
	       100*math.abs((mu_list[i]-mu_list_2[i])/mu_list[i]),"\n")
end
mu_txt:close()

k_txt = io.open("uncertainty_k.txt","w")
k_txt:write("updateThermoFromPT  ", "  updateThermoFromRHOU:",
	"  [difference","  percentage]","    updateThermoFromPS:",
	"  [difference","  percentage]","\n")
for i in pairs(k_list) do
	k_txt:write(k_list[i],"	 ",k_list_1[i],"  ",
	       math.abs(k_list[i]-k_list_1[i]),"  ",
	       100*math.abs((k_list[i]-k_list_1[i])/k_list[i]),"  ",
	       k_list_2[i],"  ",math.abs(k_list[i]-k_list_2[i]),"  ",
	       100*math.abs((k_list[i]-k_list_2[i])/k_list[i]),"\n")
end
k_txt:close()







