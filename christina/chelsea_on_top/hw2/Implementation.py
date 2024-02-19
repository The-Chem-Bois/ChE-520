#Christina Mohan
#20844467

from Shortcut_Methods import *

#Reactor
#Compute the linear mass balances for each component in the reactor based on fixed conversion model
#EL
eta = 0.07
K = 0.2
muEL1 = 96
muEL2 = (1-eta)*muEL1

#W
muW1 = 0.6*muEL1
muW2 = muW1 - eta*muEL1

#EA
muEA1 = 0
muEA2 = eta*muEL1 + muEA1

#DEE
muDEE2 = (muEA2)**2/(muW2*K) #Based on equilibrium equation

#print('EL',muEL1,muEL2, '\nW', muW1, muW2, '\nEA', muEA1, muEA2, '\nDEE', muDEE2)

#Flash
#Compute split fractions for flash unit
antoine_coeffs = np.array(
 [
        [8.10765, 1750.286, 235.0], # water
        [8.04494, 1554.3, 222.65], # ethanol
        [6.74756, 585, 255], # ethylene
        [7.4021, 1391.4, 273.16], # DEE
    ]
)

antoine_coeffs[:, :2] *= np.log(10) # Convert A, B from log base 10 to log base e.
antoine_coeffs[:, 2] -= 273.15 # Convert C coefficient from C to K.

A = antoine_coeffs[:,0]
B = antoine_coeffs[:,1]
C = antoine_coeffs[:,2]

eps_flash = shortcut_flash(51379.22, 310, A,B,C, 0.5, 3)

#print('Flash')
#print('W, EA, EL, DEE')
#print(eps_flash)

#Absorber
#Compute split fractions for absorber unit
eps_abs_V = shortcut_absorber(1,51004.2,310,A,B,C)[0]
eps_abs_L = shortcut_absorber(1,51004.2,310,A,B,C)[1]

#print('Absorber')
#print('W - vapor, EA - vapor, EL - vapor, DEE - vapor ; W - liquid, EA - liquid, EL - liquid, DEE - liquid')
#print(eps_abs_V, eps_abs_L)

#Tearing Algorithm
#Ethylene
mu01el = 96
mu1el = mu01el/(1-(eps_abs_V[2]*eps_flash[2]*(1-eta))-((1-eps_flash[2])*(1-eta))-((1-eps_abs_V[2])*eps_flash[2]*(1-eta)))
mu2el = (1-eta)*mu1el
mu31el = eps_flash[2]*mu2el
mu32el = (1-eps_flash[2])*mu2el
mu41el = (eps_abs_V[2])*mu31el
mu42el = (1-eps_abs_V[2])*mu31el
mu51el = mu41el
mu52el = 0
mu6el = mu32el + mu42el
mu71el = mu6el
mu72el = 0
mu81el = mu71el
mu82el = 0
mu91el = 0
mu92el = 0

mass_EL = np.array([mu01el,mu1el,mu2el,mu31el,mu32el,mu41el,mu42el,mu51el,mu52el,mu6el,mu71el, mu72el, mu81el, mu82el, mu91el, mu92el])
#print('Ethylene Linear Mass Balance')
#print('mu01,mu1,mu2,mu31,mu32,mu41,mu42,mu51,mu52, mu6, mu71, mu72, mu81, mu82, mu91, mu92')
#print(mass_EL)

#Water
mu02w = 0.6*mu01el
mu1w = mu02w/(1-(eps_abs_V[0]*eps_flash[0]*(1-eta)))
mu2w = (1-eta)*mu1w
mu31w = eps_flash[0]*mu2w
mu32w = (1-eps_flash[0])*mu2w
mu41w = (eps_abs_V[0])*mu31w
mu42w = (1-eps_abs_V[0])*mu31w
mu51w = mu41w
mu52w = 0
mu6w = mu32w + mu42w
mu71w = 0.1*mu6w
mu72w = 0.9*mu6w
mu81w = mu71w
mu82w = mu72w

mass_W = np.array([mu02w,mu1w,mu2w,mu31w,mu32w,mu41w,mu42w,mu51w,mu52w,mu6w,mu71w, mu72w, mu81w, mu82w])
#print('Water Linear Mass Balance')
#print('mu02,mu1,mu2,mu31,mu32,mu41,mu42,mu51,mu52, mu6, mu71, mu72, mu81, mu82')
#print(mass_W)

#Ethanol
mu1ea = 1/(1-(eps_abs_V[1]*eps_flash[1]*(1+eta))-((1-eps_flash[1])*(1+eta)+((1-eps_abs_V[1])*eps_flash[1]*(1+eta)))*0.005*0.995)
mu2ea = (1+eta)*mu1ea
mu31ea = eps_flash[1]*mu2ea
mu32ea = (1-eps_flash[1])*mu2ea
mu41ea = (eps_abs_V[1])*mu31ea
mu42ea = (1-eps_abs_V[1])*mu31ea
mu51ea = mu41ea
mu52ea = 0
mu6ea = mu32ea + mu42ea
mu71ea = 0.995*mu6ea
mu72ea = 0.005*mu6ea
mu81ea = 0.005*mu71ea
mu82ea = 0.995*mu71ea
mu91ea = 0.995*mu82ea
mu92ea = 0.005*mu82ea

mass_EA = np.array([mu1ea,mu2ea,mu31ea,mu32ea,mu41ea,mu42ea,mu51ea,mu52ea,mu6ea,mu71ea, mu72ea, mu81ea, mu82ea, mu91ea, mu92ea])
#print('Ethanol Linear Mass Balance')
#print('mu1,mu2,mu31,mu32,mu41,mu42,mu51,mu52, mu6, mu71, mu72, mu81, mu82, mu91, mu92')
#print(mass_EA)

#Diethyl Ether
mu1d = 1/(1-(eps_abs_V[3]*eps_flash[3]*K*(mu2ea**2)/mu2w)-((1-eps_flash[3])*K*(mu2ea**2)/mu2w +((1-eps_abs_V[3])*eps_flash[3]*K*(mu2ea**2)/mu2w)))
mu2d = K*(mu2ea**2)/mu2w
mu31d = eps_flash[3]*mu2d
mu32d = (1-eps_flash[3])*mu2d
mu41d = (eps_abs_V[3])*mu31d
mu42d = (1-eps_abs_V[3])*mu31d
mu51d = mu41d
mu52d = 0
mu6d = mu32d + mu42d
mu71d = mu6d
mu72d = 0
mu81d = 0
mu82d = mu71d
mu91d = mu82d
mu92d = 0

mass_DEE = np.array([mu1d,mu2d,mu31d,mu32d,mu41d,mu42d,mu51d,mu52d,mu6d,mu71d, mu72d, mu81d, mu82d, mu91d, mu92d])
#print('Diethyl Ether Linear Mass Balance')
#print('mu1,mu2,mu31,mu32,mu41,mu42,mu51,mu52, mu6, mu71, mu72, mu81, mu82, mu91, mu92')
#print(mass_DEE)

#Temperature Levels
#Stream 1
# mu1 = np.array([mu1w,mu1ea,mu1el,mu1d])
# T1 = bubble_point(750,300,mu1,A,B,C,Find_T=True)

#Stream 31 and 32
mu2 = np.array([mu2w,mu2ea,mu2el,mu2d])
#mu31_OLD= np.array([mu31w,mu31ea,mu31el,mu31d]) 
#mu32_OLD= np.array([mu32w,mu32ea,mu32el,mu32d])
# breakpoint()
flash1 = case1_flash(0.5,68.5*750,393, mu2, A,B,C, 2)
mu31 = flash1[2]*flash1[3]
mu32 = flash1[1]*flash1[4]
# print(mu31,mu32)

#Stream 03, 41, 42
absorber1 = absorber(mu31,68.5*750, 300,A,B,C,2)
mu03 = absorber1[1]
mu41 = absorber1[2]
mu42 = absorber1[3]
#print(mu03, mu41, mu42)

#Stream 6
mu6 = mu32+mu42
T6 = bubble_point(68.5*750, 350, mu6, A, B, C, Find_T=True)

#Stream 71 and 72
distillation1 = distillation(mu6,T6,20*750,(1,0.995),(0,0.005),A,B,C)
mu71 = distillation1[1]
mu72 = distillation1[0]
Tdist1 = distillation1[4]
#print(mu71,mu72)

#Stream 81 and 82
distillation2 = distillation(mu71,Tdist1,10*750,(3,0.995),(1,0.005),A,B,C)
mu81 = distillation2[1]
mu82 = distillation2[0]
Tbot2 = distillation2[5]
#print(mu81,mu82)

#Stream 91 and 92
distillation3 = distillation(mu82,Tbot2,1*750,(1,0.995),(2,0.005),A,B,C)
mu91 = distillation3[1]
mu92 = distillation3[0]
#print(mu91,mu92)
