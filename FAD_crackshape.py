import numpy as np
import matplotlib.pyplot as plt
import math
import sys

def fad_diag(Lr):
    return (1+0.5*Lr**2)**(-1/2)*(0.3+0.7*np.exp(-0.6*Lr**2))

def fad_assess_line(Lr_assessment,Kr_assessment,K_ic,Ki_secondary,Lr):
    return (Kr_assessment-Ki_secondary/K_ic)/Lr_assessment * Lr + Ki_secondary/K_ic

def N_cycle(a_final,a_initial,geometryfactor_Y,sigma_range,C_paris,m):
    return (a_final**(1-m/2)-a_initial**(1-m/2))/(C_paris*(1-m/2)*(geometryfactor_Y*sigma_range*math.sqrt(math.pi))**m)

def cycle2length(a_initial,Number_of_cycle,geometry_factor,sigma_range):
    return (a_initial**(1-m/2)+Number_of_cycle*C_paris*(1-m/2)*(geometry_factor*sigma_range*math.sqrt(math.pi))**m)**(2/(2-m))

""" Draw FAD for Steel B """
""" Given parameters are written below """
""" Material parameters """
sigma_y = 350 # Yielding Stress [MPa]
sigma_u = 500 # Ultimate Tensile Strength [MPa]
K_ic = 130 # Fracture toughness [MPa*m^0.5]
sigma_res = 70 # Residual stress [MPa]
C_paris = 10**(-11) # Parameters of Paris's law
m = 3

a_i  = 0.003 # Initial Crack length [m]
a_detectable = 0.004 # Detectable Crack length [m]
a_repair_max = 0.015 # Length of crack that will be repaired 
D_pipe = 2.5 # Pipe diameter [m]
h_n = 568 # Static Water head [m]
h_max = 744 # Water hammer head [m] 
rho = 1000 # Wate density [kg/m^3]
g = 9.8 # Acceleration of gravity [m/s^2]

""" Geometric parameters """
a_to_c =0.3
Y_elliptic  = 1.12/math.sqrt(1+1.464*(a_to_c**1.65)); Y_penet = 1.0
Y_ave = (Y_elliptic+Y_penet)/2 # Geometry correction factor (penetrated crack)
print('Geometric Calibration Factor Elliptic %f, Penetrated %f, Averaged %f'%(Y_elliptic,Y_penet,Y_ave))
""" Cycles """
Cycle_per_year = 7500; Service_year = 50
Cycle_total = Cycle_per_year*Service_year

""" Water Pressure Inside """
p_normal = rho*g*h_n #[Pa]
p_max = rho*g*h_max #[Pa]

""" Safety Parameters """
stress_limit = 0.6; SoF = 5
Inspection_period = 10 # how many times do we inspect in 50 years

""" Determine the thickness of the wall """
wall_t = (p_max*D_pipe)/(2*stress_limit*sigma_y*10**6)

""" Calculate the nominal stress range """
sigma_0 = p_normal*D_pipe/(2*wall_t)*(10**-6) # [MPa]
sigma_max = p_max*D_pipe/(2*wall_t)*(10**-6) # [MPa]

""" Calculate the critical crack length """
a_crt = (1/math.pi)*(K_ic/Y_ave/sigma_0)**2 # [m]
print(a_crt)

""" Calculate the number of cycles to failure """
N_failure = N_cycle(a_crt,a_i,Y_ave,sigma_0,C_paris,m)
print('Initial Estimation of the Wall thickness %f' %wall_t)
print('Initial Estimation of the number of cycles to failure %d' %N_failure)

if Cycle_total/Inspection_period < N_failure/SoF:
    sys.exit()
else:
    """ Modify the wall thickness """
    while Cycle_total/Inspection_period >= N_failure/SoF:
        wall_t = wall_t*1.0001
        """ Calculate the nominal stress range """
        sigma_0 = p_normal*D_pipe/(2*wall_t)*(10**-6) # [MPa]

        """ Calculate the critical crack length """
        a_crt = (1/math.pi)*(K_ic/Y_ave/sigma_0)**2 # [m]
        
        """ Calculate the number of cycles to failure, repairment and detectable size """
        N_failure = N_cycle(a_crt,a_i,Y_ave,sigma_0,C_paris,m)

wall_t_mm = wall_t*1000
print('Modified Wall Thickness:%f [mm]' % wall_t_mm)
print('Modified Number of cycles to failure:%f [mm]' %N_failure)
print('Modified Tensile Stree:%f [MPa]' % sigma_0)
#sys.exit()

N_inspection = N_failure/SoF # Number of cycles to inspection

""" Calculate the crack length for inspection """
a_repair = cycle2length(a_i,N_inspection,Y_ave,sigma_0)
print('Number of cycles to be repaired at SoF of %f: %f '% (SoF,N_inspection))
print('Crack Length to be repaired at SoF of %f (N_failure base): %f'% (SoF,a_repair))

if a_repair < a_detectable:
    a_repair = a_detectable

""" Calculate the number of cycles that a detectable length appears """
N_detectable = N_cycle(a_detectable,a_i,Y_ave,sigma_0,C_paris,m)
Detectable_num = Cycle_total/N_detectable
Appear_freq = Service_year*12/Detectable_num
sof_at_detect = N_failure/N_detectable

print('Number of cycles to the detectable crack length %f '% N_detectable)
print('Detectable crack appears every %f mothns'% Appear_freq)
print('SoF when a crack becomes detectable %f '% sof_at_detect)

N_when_detected = Cycle_per_year*(Service_year/Inspection_period)*math.ceil(N_detectable/(Cycle_per_year*(Service_year/Inspection_period)))
a_found = cycle2length(a_i,N_when_detected,Y_ave, sigma_0)*1000
sof_found = N_failure/N_when_detected
print('Crack size when it is found in a regular inspection %f [mm]'% a_found)
print('SoF when found %f '% sof_found)

a_repair = a_found/1000


""" Draw FAD at a crack length to be repaired/detected """
sigma_0 = sigma_max
""" Determine the cut off value of Lr """
Lr_max = (sigma_y+sigma_u)/(2*sigma_y)
Lr = np.arange(0,Lr_max,0.0001,float)
Kr_cutoff = np.array([0,fad_diag(Lr_max)])
Lr_cutoff = np.array([Lr_max,Lr_max])

""" Select the diagram to be used """
Kr_diag = fad_diag(Lr)

""" Decide the assessment point """
Ki_prime,Ki_secondary = Y_ave*sigma_0*np.sqrt(np.pi*a_repair), Y_ave*sigma_res*np.sqrt(np.pi*a_repair)

Lr_assess = 2*sigma_0/(sigma_u+sigma_y)
Kr_assess = 1/K_ic * (Ki_prime + Ki_secondary)

Kr_base = fad_assess_line(Lr_assess,Kr_assess,K_ic,Ki_secondary,Lr)

""" Decide the crossing point """
idLr = np.argwhere(np.sign(np.round(Kr_diag - Kr_base, 4)) == 0)

############# Y sensitivity check ##################################################################################################
""" Decide the assessment point positive error case"""
Y_plus = 1.2
Ki_prime_plus,Ki_secondary_plus = Y_plus*sigma_0*np.sqrt(np.pi*a_repair), Y_plus*sigma_res*np.sqrt(np.pi*a_repair)

Lr_assess_plus = 2*sigma_0/(sigma_u+sigma_y)
Kr_assess_plus = 1/K_ic * (Ki_prime_plus + Ki_secondary_plus)

Kr_plus = fad_assess_line(Lr_assess_plus,Kr_assess_plus,K_ic,Ki_secondary_plus,Lr)

""" Decide the crossing point """
idLr_plus = np.argwhere(np.sign(np.round(Kr_diag - Kr_plus, 4)) == 0)

""" Decide the assessment point positive error case"""
Y_minus = 0.8
Ki_prime_minus,Ki_secondary_minus = Y_minus*sigma_0*np.sqrt(np.pi*a_repair), Y_minus*sigma_res*np.sqrt(np.pi*a_repair)

Lr_assess_minus = 2*sigma_0/(sigma_u+sigma_y)
Kr_assess_minus = 1/K_ic * (Ki_prime_minus + Ki_secondary_minus)

Kr_minus = fad_assess_line(Lr_assess_minus,Kr_assess_minus,K_ic,Ki_secondary_minus,Lr)

""" Decide the crossing point """
idLr_minus = np.argwhere(np.sign(np.round(Kr_diag - Kr_minus, 4)) == 0)

""" Plot the results """
figure_ndt = plt.figure()
plt.plot(Lr, Kr_diag, color='red', label='FAD')
plt.plot(Lr_cutoff,Kr_cutoff, color='red',label='')
plt.xlim(0,1.5); plt.ylim(0,1.2)
plt.xlabel("Lr")
plt.ylabel("Kr")

SoF_fad = 0
####### Base case
if len(idLr)==0:
    plt.plot(Lr_max,fad_assess_line(Lr_assess,Kr_assess,K_ic,Ki_secondary,Lr_max),'ms', ms=5, label='Failure Point', color='green')
    plt.text(Lr_max,fad_assess_line(Lr_assess,Kr_assess,K_ic,Ki_secondary,Lr_max),'({Lr},{Kr})'.format(Lr=np.round(Lr_max,decimals=3),Kr=np.round(fad_assess_line(Lr_assess,Kr_assess,K_ic,Ki_secondary,Lr_max),decimals=3)), fontsize=10)
    SoF_fad = Lr_max/Lr_assess
else:
    plt.plot(Lr[idLr], Kr_base[idLr], 'ms', ms=5, label='Failure Point', color='green')
    for i in idLr.ravel():
        plt.text(Lr[i], Kr_base[i], '({Lr}, {Kr})'.format(Lr=np.round (Lr[i],decimals = 3), Kr=np.round(Kr_base[i],decimals=3)), fontsize=10,horizontalalignment = 'left')
    SoF_fad = Lr[idLr]/Lr_assess
print(SoF_fad)
plt.plot(Lr_assess, Kr_assess, 'ms', ms=5, label='Assessment point', color='black')
plt.text(Lr_assess,Kr_assess,'({Lr},{Kr})'.format(Lr=np.round(Lr_assess,decimals=3),Kr=np.round(Kr_assess,decimals=3)), fontsize=10)
plt.plot(Lr, Kr_base, color='black', label='Base (Y = 1.01)')

###### Positive error
if len(idLr_plus)==0:
    plt.plot(Lr_max,fad_assess_line(Lr_assess_plus,Kr_assess_plus,K_ic,Ki_secondary_plus,Lr_max),'ms', ms=5, label='', color='green')
    plt.text(Lr_max,fad_assess_line(Lr_assess_plus,Kr_assess_plus,K_ic,Ki_secondary_plus,Lr_max),'({Lr},{Kr})'.format(Lr=np.round(Lr_max,decimals=3),Kr=np.round(fad_assess_line(Lr_assess_plus,Kr_assess_plus,K_ic,Ki_secondary_plus,Lr_max),decimals=3)), fontsize=10)
else:
    plt.plot(Lr[idLr_plus], Kr_plus[idLr_plus], 'ms', ms=5, label='', color='green')
    for i in idLr_plus.ravel():
        plt.text(Lr[i], Kr_plus[i], '({Lr}, {Kr})'.format(Lr=np.round (Lr[i],decimals = 3), Kr=np.round(Kr_plus[i],decimals=3)), fontsize=10,horizontalalignment = 'left')

plt.plot(Lr_assess_plus, Kr_assess_plus, 'ms', ms=5, label='', color='black')
plt.text(Lr_assess_plus,Kr_assess_plus,'({Lr},{Kr})'.format(Lr=np.round(Lr_assess_plus,decimals=3),Kr=np.round(Kr_assess_plus,decimals=3)), fontsize=10,ha='right',va='bottom')
plt.plot(Lr, Kr_plus, color='coral', label='Y = 1.2')

####### Negative error case
if len(idLr_minus)==0:
    plt.plot(Lr_max,fad_assess_line(Lr_assess_minus,Kr_assess_minus,K_ic,Ki_secondary_minus,Lr_max),'ms', ms=5, label='', color='green')
    plt.text(Lr_max,fad_assess_line(Lr_assess_minus,Kr_assess_minus,K_ic,Ki_secondary_minus,Lr_max),'({Lr},{Kr})'.format(Lr=np.round(Lr_max,decimals=3),Kr=np.round(fad_assess_line(Lr_assess_minus,Kr_assess_minus,K_ic,Ki_secondary_minus,Lr_max),decimals=3)), fontsize=10)
else:
    plt.plot(Lr[idLr_minus], Kr_minus[idLr_minus], 'ms', ms=5, label='', color='green')
    for i in idLr_minus.ravel():
        plt.text(Lr[i], Kr_minus[i], '({Lr}, {Kr})'.format(Lr=np.round (Lr[i],decimals = 3), Kr=np.round(Kr_minus[i],decimals=3)), fontsize=10,horizontalalignment = 'left')

plt.plot(Lr_assess_minus, Kr_assess_minus, 'ms', ms=5, label='', color='black')
plt.text(Lr_assess_minus,Kr_assess_minus,'({Lr},{Kr})'.format(Lr=np.round(Lr_assess_minus,decimals=3),Kr=np.round(Kr_assess_minus,decimals=3)), fontsize=10,ha='left',va='top')
plt.plot(Lr, Kr_minus, color='blue', label='Y = 0.8')


plt.legend()
plt.show()
figure_ndt.savefig("FAD_crackshape.png")
