--GasModel Steam: Eilmer4 Input file for the Uncertainty Analysis 
--				   of thermo-update Functions 

--Considered thermo-update Functions:
-- 		-updateThermoFromRHOU
-- 		-updateThermoFromRHOT
--		-updateThermoFromRHOP
--		-updateThermoFromPS
--		-updateThermoFromHS


gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}

file = io.open("updateThermoFromRHOU.txt","w")
file:write("First run:")
Q.rho = 2.532197743e-2 
Q.T= 300 
Q.quality = 1
print("--1. updateThermoFromPT:")
gmodel:updateThermoFromRHOT(Q)
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Density [kg/m3]=", Q.rho, "Interal Energy [J/kg]=", Q.u, "Sound Speed [m/s]=",
	 Q.a,"\n")
file:write("Density [kg/m3]=")
file:close()





