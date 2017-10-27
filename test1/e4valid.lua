--Steam GasModel: Validation Eilmer4 Input file for Gasmodel Steam

--Validated Lua GasModel functions:
-- 		-updateThermoFromPT
-- 		-updateThermoFromRHOU
-- 		-updateSoundSpeed
-- 		-updateTransCoeffs
-- 		-intEnergy
-- 		-enthalpy
-- 		-entropy
-- 		-Cv
-- 		-Cp
-- 		-gasConstant


gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}

Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--1. updateThermoFromPT:")
gmodel:updateThermoFromPT(Q)
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Density [kg/m3]=", Q.rho, " Specific Interal Energy [J/kg]=", Q.u, " Sound Speed [m/s]=",
	 Q.a,"\n")

Q.rho = 5.0
Q.u = 2.8e6
print("--2. updateThermoFromRHOU:")
gmodel:updateThermoFromRHOU(Q)
print("Density [kg/m3]=",Q.rho," Specific Interal Energy [J/kg]=",Q.u," quality=",
		Q.quality)
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"\n")

Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--3. updateSoundSpeed:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
gmodel:updateSoundSpeed(Q)
print("Sound Speed [m/s]=", Q.a,"\n")
Q.p = 0.3e6 
Q.T= 650 
Q.quality = 1
print("--4. updateTransCoeffs:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
gmodel:updateTransCoeffs(Q)
print("Dynamic Viscosity [Pa.s]=", Q.mu, "Thermal Conductivity [W/m/K]=", Q.k, "\n")
Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--5. intEnergy:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("u [J/kg]=", gmodel:intEnergy(Q),"\n")
print("--6. enthalpy:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("h [J/kg]=", gmodel:enthalpy(Q),"\n")
print("--7. entropy:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("s [J/kg]=", gmodel:entropy(Q),"\n")
print("--8. Cv:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Cv [J/K/kg]=", gmodel:Cv(Q),"\n")
print("--9. Cp:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Cp [J/K/kg]=", gmodel:Cp(Q),"\n")
print("--10. gasConstant:")
print("Temperature [K]=",Q.T,"Pressure [Pa]=",Q.p,"quality=",Q.quality)
print("Specific Gas Constant [J/K/kg]=", gmodel:gasConstant(Q),"\n")



