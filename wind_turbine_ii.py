from math import cos, sin, atan, pi, exp, acos, tan
from numpy import zeros, linspace

def initial_values(n, lambda_d, Cld, R, Z):
    lambdai = zeros(n)      # Tip speed ratio at ith position
    phi_i = zeros(n)            # Angle of the relative wind velocity at ith position
    Ci = zeros(n)               # Chord length at ith position
    rR = linspace(0.20, 0.90, n)    # r/R coordinates

    for i in range(n):
        lambdai[i] += lambda_d*rR[i]
        phi_i[i] += (2/3)*atan(1/lambdai[i])
        Ci[i] += 8*pi*rR[i]*R/(Z*Cld)*(1-cos(phi_i[i]))
    
    return [lambdai, Ci, rR]

def blade_geometry(initial_blade, Cd_Cl, n, Cld, alpha_d, R, Z):
    
    tip_ratio = initial_blade[0]
    initial_Ci = initial_blade[1]
    alpha_d *= pi/180
    
    imax = 10
    a = zeros((imax, n))
    ap = zeros((imax, n))
    sigma = zeros(n)
    phi = zeros((imax, n))
    Cd = Cld*Cd_Cl
    Cx = zeros((imax, n))
    Cy = zeros((imax, n))
    Cp = 0
    # Ct = zeros((imax, n))
    theta = zeros((imax, n))
    Ft = zeros((imax, n))
    Chord = zeros(n)
    rR = linspace(0.2, 0.90, n)  # r/R coordinates
    TSR = (tip_ratio[1]-tip_ratio[0])/(rR[1]-rR[0])*(1-rR[0])+tip_ratio[0]


    for j in range(n):              #Bucle para calcular la geometria en cada posicion r/R
        a[0][j] = 0
        ap[0][j] = 0

        phi[0][j] = atan((1 - a[0][j])/((1 + ap[0][j])*tip_ratio[j]))
        Ft[0][j] = 2/pi*acos(exp(-Z/2*(1-rR[j])/(rR[j]*sin(phi[0][j]))))
        sigma[j] = Z*initial_Ci[j]/(2*pi*rR[j]*R)
        theta[0][j]  = phi[0][j] - alpha_d

        Cx[0][j] = Cld*cos(phi[0][j]) + Cd*sin(phi[0][j])
        Cy[0][j] = Cld*sin(phi[0][j]) - Cd*cos(phi[0][j])
        
        
        for i in range(1, imax):    #Bucle para convergencia de resultados
            
            a[i][j] = sigma[j]*Cx[i-1][j]/(4*Ft[i-1][j]*(sin(phi[i-1][j]))**2 + sigma[j]*Cx[i-1][j])
            ap[i][j] = sigma[j]*Cy[i-1][j]/(4*Ft[i-1][j]*sin(phi[i-1][j])*cos(phi[i-1][j]) - sigma[j]*Cy[i-1][j])

            phi[i][j] = atan((1 - a[i][j])/((1 + ap[i][j])*tip_ratio[j]))
            theta[i][j]  = phi[i][j] - alpha_d
            Ft[i][j] = 2/pi*acos(exp(-Z/2*(1-rR[j])/(rR[j]*sin(phi[i][j]))))
            
            Cx[i][j] = Cld*cos(phi[i][j]) + Cd*sin(phi[i][j])
            # Ct = sigma*((1-a[i][j])**2)*Cx[i-1][0]/((sin(phi[i][j]))**2)
            Cy[i][j] = Cld*sin(phi[i][j]) - Cd*cos(phi[i][j])
                
        
    for k in range(n):
        Chord[k] = 8*pi*rR[k]*R/(Z*Cld)*(1-cos(phi[-1][k]))
        initial_Ci[k] = Chord[k]
        Cp = Cp + 8/((tip_ratio[-1])*n)*(Ft[-1][k]*(sin(phi[-1][k])**2))*(cos(phi[-1][k])-tip_ratio[k]*sin(phi[-1][k]))*(sin(phi[-1][k])+tip_ratio[k]*cos(phi[-1][k]))*(1-drag_to_lift*(1/tan(phi[-1][k])))*(tip_ratio[k])**2

        
        
    
    return [rR, Chord, theta[-1], phi[-1], Cx[-1], Cy[-1], Cp]

def print_geometry_values(geometry):
    rR = geometry[0]
    Blade_Chord = geometry[1]
    Twist_Angle = geometry[2]*180/pi
    phi = geometry[3]*180/pi
    Fx_Coeff = geometry[4]
    Fy_Coeff = geometry[5]
    Cp = geometry[6]
    

    with open('blade_dimensions.txt', 'a') as file:
        file.write('\n\nr/R\tC[m]\tTwist Angle[deg]\tphi[deg]\tNormal Force coefficient\tAxial Force coefficient')
        for i in range(len(rR)):
            file.write('\n{}\t{}\t{}\t{}\t{}\t{}'.format(rR[i], Blade_Chord[i], Twist_Angle[i], phi[i], Fx_Coeff[i], Fy_Coeff[i]))
        file.write('\nPower Coefficient\t{}'.format(Cp))

    with open('average_cp.txt', 'a') as file:
        file.write('Power Coefficient\t{}\n'.format(Cp))

def blade_design(N, Rd, Cld, alpha_d, Z, tip_ratio_d, Lift_to_drag_d):
    Initial_geometry = initial_values(N, tip_ratio_d, Cld, Rd, Z)
    Final_geometry = blade_geometry(Initial_geometry, Lift_to_drag_d, N, Cld, alpha, Rd, Z)
    print_geometry_values(Final_geometry)

N = 20
R = 0.75
Cl = 1.08
alpha = 8
B = 3
tip = range(2,26)
drag_to_lift = 1/51.89

for i in tip:
    blade_design(N, R, Cl, alpha, B, i, drag_to_lift)

print('PROGRAMA EJECUTADO CON EXITO :D')


