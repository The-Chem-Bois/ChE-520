#Reactor

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
muDEE2 = (muEA2)**2/(muW2*K)

print('EL',muEL1,muEL2, '\nW', muW1, muW2, '\nEA', muEA1, muEA2, '\nDEE', muDEE2)