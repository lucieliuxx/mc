import numpy as np 
import matplotlib.pyplot as plt 
plt.rcParams.update({'font.size': 13})
plt.rcParams['lines.linewidth'] = 2

# number constants
C = (4*np.pi/3) * (6-18/5*2**(2/3))
A = C**(1/3)*(6**(2/3)+12/6**(1/3))
D = (4*np.pi/3)*(9/2/np.pi)**(2/3) + 12*(2*np.pi/9)**(1/3)
Lx, Ly, Lz = 20,20,28
#print("C = {0:.2f}, A = {1:.2f}, D = {2:.2f}".format(C,A,D))

# functions
def f(sig):
    fval = (sig**2 + (1+sig)**(5/3) - 1)/5
    fval += (sig + 1)*((sig + 1) - 2 * (sig + 1)**(2/3) + 1) 
    return fval

def E_bilayer(qh,d,Lx,Ly,Lz,eps, sig):
    e_bilayer = 4*np.pi/3 * d**3 * qh**2 * Lx * Ly *(sig + 1)/sig
    e_bilayer += 4*eps*Lx*Ly - eps*3*Lx*Ly*Lz
    return e_bilayer

def E_soluble(qh,d,Lx,Ly,Lz,eps):
    return -2*np.sqrt(np.pi)*(d*qh**2)**(3/2) *(Lx*Ly*Lz)**(1/2) - (3*Lz - 24*d) * Lx*Ly * eps

def E_micelle(qh,d,Lx,Ly,Lz,eps,sig):
    #return A * qh**(2/3) * d * L**2 * eps**(2/3) - 3*L**3*eps
    return D * qh**(2/3) * d * Lx*Ly * f(sig)**(1/3) / sig * eps**(2/3) - 3*eps*Lx*Ly*Lz

def U_micelle_min(eps,Lx,Ly,Lz,sig,nh,z,rhom):
    return eps*Lx*Ly*Lz*(-3 + D*f(sig)**(1/3) * nh**(1/3) / 2 / sig * (z**2 * rhom**3 / eps)**(1/3) )

def U_lamellar_min(eps,Lx,Ly,Lz,sig,nh,z,rhom):
    return eps*Lx*Ly*Lz*(-3 +  (18*np.pi * (sig + 1)/sig * z**2 * nh * rhom**3 / eps)**(1/3) )

def r_micelle_inner(qh,eps,sig):
    #return (6*eps/(C*qh**2))**(1/3)
    return (9/(2*np.pi*f(sig)) * eps / qh**2)**(1/3)

def r_micelle_outer(qh,eps,sig):
    return r_micelle_inner(qh,eps,sig) *(1+sig)**(1/3)

def n_layers_opt(eps,Lz,sig,nh,z,rhom):
    return (np.pi/12 * z**2 * (nh + nh/sig) * rhom**3 * Lz**3 / eps )**(1/3)


# independent variables
z, eps = 0.6, 1
sig = 1
nh = 1
rhom = 0.25
qh = z/nh
nt = nh/sig
qt = - qh * sig

#print(C,4*np.pi/3*f(1)/1)
#print(6/C, 9/(2*np.pi*f(1)))
#print(A, D * f(sig)**(1/3)/sig)

# derived varibales
C1 = C**(1/3)*(6**(2/3)+12/6**(2/3))*eps**(2/3)/(8*np.pi)
#ds = (A/8/np.pi)**(1/2) * qh**(-2/3) * eps**(1/3) 
ds = (D/4/np.pi)**(1/2) * qh**(-2/3) * eps**(1/3) * f(sig)**(1/6)/(sig+1)**(1/2)
#print("constants = ", (A/8/np.pi)**(1/2), 4 - 2/3 * A**(3/2) / (8*np.pi)**(1/2))

"""
dE_over_L2_analytical = (4 - 2/3 * D**(3/2) / (4*np.pi)**(1/2) * ( f(sig)/(sig+1)/sig**2 )**(1/2) ) * eps
dE_over_L2_numerical = E_bilayer(qh,ds,Lx,Ly,Lz,eps,sig)/Lx/Ly-E_micelle(qh,ds,Lx,Ly,Lz,eps,sig)/Lx/Ly
print("inner radius of the micelle = {0:.2f}, outer radius = {1:.2f}".format(r_micelle_inner(qh,eps,sig),r_micelle_outer(qh,eps,sig)))
print("optimal head layer thickness dh = {0:.2f}".format(ds))
print("min dE/L2 at optimal d: analytical {0:.2f}, numerical {1:.2f}".format(dE_over_L2_analytical,dE_over_L2_numerical))
"""

print("n* = {0:.1f}".format(n_layers_opt(eps,Lz,sig,nh,z,rhom)))
print("Nh = {0:.0f}, Nt = {1:.0f}, Nw = {2:.0f}".format(rhom*Lx*Ly*Lz*nh, rhom*Lx*Ly*Lz*nt, Lx*Ly*Lz - rhom*Lx*Ly*Lz*(nh + nt)))
print("z = {0:.2f}, eps = {1:.2f}".format(z,eps))
dU_opt = U_lamellar_min(eps,Lx,Ly,Lz,sig,nh,z,rhom) - U_micelle_min(eps,Lx,Ly,Lz,sig,nh,z,rhom)
print("u_lam - u_mic = {0:.2f}".format(dU_opt/Lx/Ly/Lz))

"""
# plots
plt.figure(3)
s = np.linspace(0.02,2,100)
dU = U_lamellar_min(eps,Lx,Ly,Lz,s,nh,z,rhom) - U_micelle_min(eps,Lx,Ly,Lz,s,nh,z,rhom)
plt.plot(s,U_lamellar_min(eps,Lx,Ly,Lz,s,nh,z,rhom)/Lx/Ly/Lz,color='tab:red',label='lamellar')
plt.plot(s,U_micelle_min(eps,Lx,Ly,Lz,s,nh,z,rhom)/Lx/Ly/Lz,color='tab:blue',label='micellar')
plt.plot(s,dU/Lx/Ly/Lz,color='k',label='lamellar - micellar')
plt.legend(loc='best')
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\beta U \delta^3/L_xL_yL_z$')
plt.tight_layout()
plt.savefig('lamellar_vs_micellar.png',dpi=200)

plt.figure(2)
dpts = np.linspace(0,4,100)
ddots = np.linspace(1,4,4)
#plt.plot(dpts,E_bilayer(qh,dpts,L,eps)/L**2-E_soluble(qh,dpts,L,eps)/L**2)
#plt.plot(dpts,(E_soluble(qh,dpts,L,eps)+3*L**3*eps)/L**2,color='tab:red',label='soluble')
#plt.plot(ddots,(E_soluble(qh,ddots,L,eps)+3*L**3*eps)/L**2,'o',color='tab:red')
plt.plot(dpts,(E_bilayer(qh,dpts,L,eps,sig)+3*L**3*eps)/L**2,color='tab:blue',label='bilayer')
plt.plot(ddots,(E_bilayer(qh,ddots,L,eps,sig)+3*L**3*eps)/L**2,'o',color='tab:blue')
plt.plot(dpts,(E_micelle(qh,dpts,L,eps,sig)+3*L**3*eps)/L**2,color='orange',label='micelle')
plt.plot(ddots,(E_micelle(qh,ddots,L,eps,sig)+3*L**3*eps)/L**2,'o',color='orange')
plt.legend(loc='best')
plt.xlabel(r'layer thickness $d/\delta$')
plt.ylabel(r'$(E - E_{water box})/L^2$')
plt.tight_layout()
plt.savefig('phase.png',dpi=200)

plt.figure(1)
sigpts = np.linspace(0.02,2,100)
plt.plot(sigpts,f(sigpts)/sigpts**2 ,color='black')
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$f(\sigma)$')
plt.tight_layout()
plt.savefig('f(sig).png',dpi=200)
"""
