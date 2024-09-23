import math 
import ROOT,array
import numpy as np
#add spectum files
dir="/project/6004969/trevorh/cross/flux_files/"
b8=dir+"Boron8spectrum.root"
hep=dir+"Hepspectrum.root"

c1=ROOT.TCanvas("c1")
c1.SetLogy(1)

#open Hep and Boron8 flux files 
def Fluxopener(file,graph):
    f = ROOT.TFile(file,"READ")
    solarflux = f.Get(graph)
    solarflux.SetLineWidth(1)
    solarflux.SetLineColor(ROOT.kBlue)
    solarflux.Draw("AL")
    return solarflux

#function computes cross section using 2009 model
def Sigma(E_n): #insert energy of incident neutrio
    #From DUne F. Capozzi et al. (DUNE Collaboration), (2019),supplemental_material.pdf (aps.org)
    G0f=1.1663787e-11 ## MeV^-2 pg 313 griffiths
    vud=0.9738 ## also Grifiths 0.97270+-0.00014https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf
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
        sigma+= math.pow(G0f*hbar*c*vud,2)*(1/math.pi)*matrix[i]*E_e*p_e*Fermi(18,E_e)
    return sigma

#function computes cross section using 1998 model
def Sigma_old(E_n): #insert energy of incident neutrio
    #From M. Bhattacharya (1998), https://journals.aps.org/prc/pdf/10.1103/PhysRevC.58.3677
    G0f=1.1663787e-11 ## MeV^-2 pg 313 griffiths
    vud=0.9738 ## also Grifiths 0.97270+-0.00014https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf
    hbar=6.582e-22 ## MeV*s wiki 
    c=3e10 #cm/s speed of light
    Gf=G0f*math.pow(hbar*c,3) ## MeV*cm^3 from griffiths 
    Qgs=1.504 ##MeV
    #need to add contributions of Sc atoms 
    deltaE = [2.28988 , 2.73038 , 2.95070,
    3.10975 , 3.14644 , 3.29300,
    3.73850 , 3.79758 , 3.84025,
    3.89800 , 3.99600 , 4.35200,
    4.38400 , 4.69700 , 4.76100,
    4.78865 , 4.84800 , 5.02700,
    5.223, 5.696,6.006]
  
    matrix= [0.90 , 1.50 , 0.11,
    0.06 , 0.04 , 0.01,
    0.16 , 0.26 , 0.01,
    0.05 , 0.11 , 0.29,
    3.84 , 0.31 , 0.38,
    0.47 , 0.36 , 0.23,
    0.03, 0.11,0.13]

    #summ over all possible states
    sigma=0
    #sigma2=0
    for i in range(len(matrix)):
        K_e=E_n-(deltaE[i]+Qgs) ## % MeV kinetic energy
        if K_e<=0: continue
        E_e=K_e+0.511 ## MeV
        p_e=math.sqrt(E_e*E_e-0.511*0.511) #MeV/c
        sigma+= math.pow(G0f*hbar*c*vud,2)*(1/math.pi)*matrix[i]*E_e*p_e*Fermi(18,E_e)
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
    min=2#MeV
    max=20#MeV
    x=np.linspace(min,max,1000)
    y = array.array("d" ,[0]* len (x))

    for i,xi in enumerate (x):
        y[i]=cross.Eval(xi)*flux.Eval(xi)
    fluxXsigma = ROOT.TGraph(len(x),x,y)
    fluxXsigma.GetXaxis().SetTitle("Energy (MeV)")
    fluxXsigma.GetYaxis().SetTitle("Cross section*flux (events*s^-1*MeV^-1")
    fluxXsigma.GetXaxis().CenterTitle()
    fluxXsigma.GetYaxis().CenterTitle()
    fluxXsigma.Draw("AL")
    c1.SaveAs("FxS.png")
    return fluxXsigma

def addTF1(graph1,graph2):
    min=2#MeV
    max=20#MeV
    x=np.linspace(min,max,1000)
    y = array.array("d" ,[0]* len (x))

    for i,xi in enumerate (x):
        y[i]=graph1.Eval(xi)+graph2.Eval(xi)
    
    sumgraph = ROOT.TGraph(len(x),x,y)
    sumgraph.GetXaxis().SetTitle("Energy (MeV)")
    sumgraph.GetYaxis().SetTitle("Flux neutrinos (s^-1*cm^-2)")
    sumgraph.GetXaxis().CenterTitle()
    sumgraph.GetYaxis().CenterTitle()
    c1.SaveAs("sumgraph.png")
    sumgraph.Draw("AL")
    return sumgraph

def IntegrateTF1(function,xmin,xmax):
    #simpsons COmposite rule WIKI
    n=2000
    h=(xmax-xmin)/n
    xj=[]
    for j in range(0,n+1):
        xj.append(xmin+j*h)
    sum=0
    for j in range(1,n/2+1):
        sum+=function.Eval(xj[2*j-2])+4*function.Eval(xj[2*j-1])+function.Eval(xj[2*j])
    integral=sum*h/3 #events per second per target atom
    print(integral)
    #convert to tons of LAr
    years=60*60*24*365 #number of seconds in a year
    avagodro=6.022e23
    molarmass=39.948 #g/mol# taken as natural abundance is this right?
    tonnes=1e6 #grams
    targets=(tonnes/molarmass)*avagodro #number of target nuclei in a tonne
    Rate=integral*years*targets
    print("Revised rate in units events per second per ton",Rate)
    return Rate

   
def make_cross_tgraph(model):
    #create cross section TF1 object 
    x=np.linspace(3.5,20,1000)
    y = array . array ("d" ,[0]* len (x))# array of zeros .
    for i,xi in enumerate (x):
        if model=="new":
            y[i] = Sigma(xi)
        else:
            y[i] = Sigma_old(xi)
    cross= ROOT.TGraph(len(x),x,y)
    #cross.GetXaxis().SetTitle("Energy (MeV)")
    #cross.GetYaxis().SetTitle("Cross section (cm^2)")
    #cross.GetXaxis().CenterTitle()
    #cross.GetYaxis().CenterTitle()
    #cross.Draw("AL")
    #c1.SaveAs("crosssection.png")
    return cross

#start a rooft canvas

def doAll(b8_spectrum,old_new):
    flux1=Fluxopener(b8,b8_spectrum)
    flux2=Fluxopener(hep,"Hep_flux")
    solarflux=addTF1(flux1,flux2)
    cross=make_cross_tgraph(old_new)
    sigma_flux=multiplyTF1(cross,solarflux)
    rate=IntegrateTF1(sigma_flux,3.5,20)
    return rate

#Draw 2 cross sections
mg = ROOT.TMultiGraph()
mg.SetTitle(" v_{e} + {}^{40}Ar -> e^{-} + {}^{40}K ; Energy (MeV); #sigma (cm^2) ")
old=make_cross_tgraph("old")
new=make_cross_tgraph("new")
old.SetLineColor(ROOT.kRed)
new.SetLineColor(ROOT.kBlue)
mg.Add(old)
mg.Add(new)
mg.Draw("AL")

leg = ROOT.TLegend(0.55, 0.35, 0.80, 0.50)
leg.SetTextSize(14)
leg.AddEntry(old, 'GT strengths (1998) model')
leg.AddEntry(new, ' GT strengths (2009) model')
#leg.Draw("same")
c1.SaveAs("cross_section_comparison.png")


