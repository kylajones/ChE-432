#Modeling of BFB reactor for plant design II hazelnut pyrolysis
#Capacity (C): Oregon produces 40,000 metric tons of hazelnuts/year
#We assume that Cascade foods is responsible for 15,000
C = 15e3 #metric tons
C = C*1000 #kg
xₛ = 0.40 #shell wt% of hazelnut
mₛ = xₛ*C #mass of hazelnut shells
Dₚ = 0.180 #[mm] avg particle diameter

##Reaction modeling
#IC's:
W₀ = 12 #[kg] initial weight of biomass
T₀ = 12 #[K] initial temperature of reactor
#time discretization:
t₀ = 0  #[s] initial time
tf = 12 #[s] final time
τ = 0.1 #[s] time step
t = [t₀:τ:tf] #[s] discretized time
#forward Euler approximation
T = T₀, W = W₀
for j = 1:length(t)-1
    dWdt = f(t,T,W)
    W(j+1) = W(j)+τ*dWdt
    #n = 1 #reaction order
    #X = ((W₀-W)/W₀)ⁿ #degree of pyrolysis
    #rₓ = k*X #rate of biomass reaction
end

function f(t,T,W)
    #initial conditions
    W₀ = 12 #[kg] initial weight of biomass
    T₀ = 12 #[K] initial temperature of reactor
    #Reation constants:
    A = 12 #pre-exponential factor
    R = 8.314 #[J/mol*K]
    E = 12 #[J/mol] activation energy
    b = 12 #[K/s] heating rate
    #time dependant paremeters:
    T = T₀+b*t #[K] temperature
    k = A*exp(-E/(R*T)) #rate constant
    dWdt = -W₀*k*exp(-E/(R*T))*(W₀-W)/W₀
    return k dWdt


