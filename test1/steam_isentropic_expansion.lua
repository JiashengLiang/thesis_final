gmodel = GasModel:new{'steam.lua'}
Q = GasState:new{gmodel}
Q.p = 270e3 -- Pa
Q.T = 403.15 -- K
Q.quality = 1.0 
gmodel:updateThermoFromPT(Q)
-- Compute enthalpy and entropy at stagnation conditions
h0 = gmodel:enthalpy(Q) + 0.5*30.62983^2
s0 = gmodel:entropy(Q)
-- Set up for stepping process
dp = 1.0 -- Pa, use 1 Pa as pressure step size
Q.p = Q.p - dp
M_tgt = 1.0
-- Begin stepping until M = M_tgt
while true do
   gmodel:updateThermoFromPS(Q, s0)
   h1 = gmodel:enthalpy(Q)
   v1 = math.sqrt(2*(h0 - h1))
   if Q.p == 269856 then
      print("If use ideal gas model with water, at p = 269856 Pa:\n")
      print("Temperature = ",Q.T,"[K],"," Velocity = ", v1,"[m/s],",
         " Specific Enthalpy = ", h1,"[J/kg]\n\n")
   end
   gmodel:updateSoundSpeed(Q)
   M1 = v1/Q.a
   if M1 >= M_tgt then
      print("Stopping at M= ", M1)
      break
   end
   Q.p = Q.p - dp
end

print("Gas state at sonic conditions are:")
print("p= ", Q.p)
print("T= ", Q.T)