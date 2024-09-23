import math 
import ROOT,array
import numpy as np

def Fluxopener():
    file="/home/trevorh/rat/src/extgen/Marley/marley-1.1.1/rootfiles/Boron_HEP_spectrum.root"  # flux in Graham
    f = ROOT.TFile(file,"READ")
    solarflux = f.Get("SolarNeutrinos")
    solarflux.SetLineWidth(1)
    solarflux.SetLineColor(ROOT.kBlue)
    solarflux.Draw("AL")
    c1.SaveAs("flux.png")
    return solarflux

def Sigma(E_n): #insert energy of incident neutrio
    #From DUne F. Capozzi et al. (DUNE Collaboration), (2019),supplemental_material.pdf (aps.org)
    G0f=1.1663787e-11 ## MeV^-2 pg 313 griffiths
    vud=0.9738 ## also Grifiths
    hbar=6.582e-22 ## MeV*s wiki 
    c=3e10 #cm/s speed of light
    Gf=G0f*math.pow(hbar*c,3) ## MeV*cm^3 from griffiths 
    Qgs=1.504 ##MeV
    matrix=[1.64, 1.49, 0.06,
    0.16, 0.44, 4.0,
    0.86, 0.48, 0.59,
    0.21, 0.48, 0.71,
    0.06, 0.14, 0.97]
    deltaE=[2.333, 2.775, 3.204,
    3.503, 3.870, 4.384,
    4.421, 4.763, 5.162,
    5.681, 6.118, 6.790,
    7.468, 7.795, 7.952]
    #summ over all possible states
    sigma=0
    #sigma2=0
    for i in range(len(matrix)):
        K_e=E_n-(deltaE[i]+Qgs) ## % MeV kinetic energy
        if K_e<=0: continue
        E_e=K_e+0.511 ## MeV
        p_e=math.sqrt(E_e*E_e-0.511*0.511) #MeV/c
        sigma+= math.pow(G0f*hbar*c*vud,2)*(1/math.pi)*matrix[i]*E_e*p_e*Fermi(19,E_e)
                            ##MeV^-2*cm^2*                         MeV *MeV/c*c =cm
                            ##Gf^2/hbar^4*c^3 =MeV^-2*cm^6/(s*cm^3)=MeV^-2*cm^2*(c)
        #sigma2+= math.pow(Gf*vud,2)/math.pow(hbar,c,3)*(1/math.pi)*matrix(i)*E_e*p_e*Fermi(19,E_e)
    return sigma
 
def Fermi(Z,E):
    #approximation from  G. K. Schenter & P. Vogel, "A Simple Approximation of the Fermi Function in Nuclear Beta Decay"
    ## https://www-tandfonline-com/doi/abs/10.13182/NSE83-A17574
    T=E-0.511 # MeV
    if T>=1.2*0.511:    
        alpha=-8.46e-2+(2.48e-2)*Z+(2.37e-4)*Z*Z
        beta=1.15e-2+(3.58e-4)*Z-(6.17e-5)*Z*Z
    else:
          alpha=-0.811+(4.46e-2)*Z+(1.08e-4)*Z*Z
          beta=0.673-(1.82e-2)*Z-(6.38e-5)*Z*Z
    
    p=math.sqrt(E**2-0.511**2) #apparently units MeV since fermi function is unitless instead of MeV/c
    fermi=E/p*math.exp(alpha+beta*math.sqrt(E/0.511-1))
    return fermi #maybe unitless maybe units c

def multiplyTF1(cross,flux):
    min=3.5#MeV
    max=20#MeV
    x=np.linspace(min,max,100)
    y = array.array("d" ,[0]* len (x))

    for i,xi in enumerate (x):
        y[i]=cross.Eval(xi)*flux.Eval(xi)
    fluxXsigma = ROOT.TGraph(len(x),x,y)
    fluxXsigma.Draw("AL")
    c1.SaveAs("FxS.png")
    return fluxXsigma

#start a rooft canvas
c1=ROOT.TCanvas("c1")
c1.SetLogy(1)
# open flux TF1 object flux
flux=Fluxopener()

#create cross section TF1 object 
x=np.linspace(3.5,15,100)
y = array . array ("d" ,[0]* len (x))# array of zeros .
for i,xi in enumerate (x):
    y[i] = Sigma(xi)
cross= ROOT.TGraph(len(x),x,y)
cross.Draw("AL")
c1.SaveAs("crosssection.png")
# create a TF1 Object that is the Flux*crosssection
sigma_flux=multiplyTF1(cross,flux)

