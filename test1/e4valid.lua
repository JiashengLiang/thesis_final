--GasModel Steam: Validation Eilmer4 Input file for Gasmodel Steam

--Validated Lua GasModel functions:
-- 		-updateThermoFromPT
-- 		-updateThermoFromRHOU
--		-updateThermoFromPS
-- 		-updateSoundSpeed
-- 		-updateTransCoeffs
-- 		-intEnergy
-- 		-enthalpy
-- 		-entropy
-- 		-Cv
-- 		-Cp
-- 		-gasConstant

print("\n\n")
print("========================================================")
print("Validation Results: Implement Steam Gasmodel in Eilmer4")
print("========================================================")
print("\n")
gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}

Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--1. updateThermoFromPT:")
gmodel:updateThermoFromPT(Q)
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nDensity [kg/m3]=", Q.rho, "\nSpecific Interal Energy [J/kg]=", Q.u, 
	"\nSound Speed [m/s]=", Q.a,"\n\n")

Q.rho = 5.0
Q.u = 2.8e6
print("--2. updateThermoFromRHOU:")
gmodel:updateThermoFromRHOU(Q)
print("Density [kg/m3]=",Q.rho,"\nSpecific Interal Energy [J/kg]=",Q.u,
	"\nquality=", Q.quality)
print("\nTemperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\n\n")

Q.T = 380.0
Q.p = 0.0035e6
print("--3. updateThermoFromPS:")
gmodel:updateThermoFromPS(Q,0.852238967e4)
print("Pressure [Pa]=",Q.p,"\ns [J/kg]=",0.852238967e4)
print("\nExpected Temperature [K]= 300", "\nActual Calculated Temperature [K]=",
		 Q.T,"\n\n")

Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--4. updateSoundSpeed:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
gmodel:updateSoundSpeed(Q)
print("\nSound Speed [m/s]=", Q.a,"\n\n")
Q.p = 0.3e6 
Q.T= 650 
Q.quality = 1
print("--5. updateTransCoeffs:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
gmodel:updateTransCoeffs(Q)
print("\nDynamic Viscosity [Pa.s]=", Q.mu, "\nThermal Conductivity [W/m/K]=",
	 Q.k, "\n\n")
Q.p = 0.0035e6 
Q.T= 300.0 
Q.quality = 1
print("--6. intEnergy:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nu [J/kg]=", gmodel:intEnergy(Q),"\n\n")
print("--7. enthalpy:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nh [J/kg]=", gmodel:enthalpy(Q),"\n\n")
print("--8. entropy:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\ns [J/kg]=", gmodel:entropy(Q),"\n\n")
print("--9. Cv:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nCv [J/K/kg]=", gmodel:Cv(Q),"\n\n")
print("--10. Cp:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nCp [J/K/kg]=", gmodel:Cp(Q),"\n\n")
print("--11. gasConstant:")
print("Temperature [K]=",Q.T,"\nPressure [Pa]=",Q.p,"\nquality=",Q.quality)
print("\nSpecific Gas Constant [J/K/kg]=", gmodel:gasConstant(Q),"\n\n")



