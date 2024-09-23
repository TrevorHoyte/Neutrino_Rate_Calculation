## this code reads in data for the boron 8 spectrum from http://www.sns.ias.edu/~jnb/SNdata/b8spectrum.html
# need to download the text file and put its location and name in dir and file locations
#SNO flux https://journals.aps.org/prc/pdf/10.1103/PhysRevC.88.025501
# The code then takes the data and normalizes the boron 8 flux and t eplus minus 3 sigma confidence interval to whatever value you would like and then saves the Tgraph object as a root file
#written by Trevor Hoyte

import ROOT,array
import numpy as np

dir="/project/6004969/trevorh/cross/flux_files/"
file1="b8spectrum.txt"
normalization=5.244e6  #  have to change the value from 5.25e6 so that the integrated flux using simpsons rule is correct# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.88.025501
 
fluxfile = open(dir+file1, 'r')
lines = fluxfile.readlines()
nbins=829 # number of data points+1
h1=ROOT.TH1D("h1","Nominal",nbins,0,16.56)
h2=ROOT.TH1D("h2","Overflow",nbins,0,16.56) #+3sigma
h3=ROOT.TH1D("h3","underflow",nbins,0,16.56) #-3sigma
energy=[]
#read through text files and extract info
for line in lines:
    data=line.strip().split()
    if len(data)==4:
        #fill hist and Tgraphs
        h1.Fill(float(data[0]),float(data[1]))
        h2.Fill(float(data[0]),float(data[2]))
        h3.Fill(float(data[0]),float(data[3]))
        energy.append(float(data[0]))
en=np.array(energy)

def Scalehist(histogram,scale,x1,x2):
    minbin=int(x1/0.02+1)
    maxbin=int(x2/0.02+1)
    factor=scale/histogram.Integral(minbin,maxbin,"width")
    histogram.Scale(factor)
    #make sure everythings right
    histogram.Draw()
    print(histogram.Integral(minbin,maxbin,"width"))
    return histogram

def histtograph(histogram):
    value=[]
    for k in range(1,nbins+1):
        value.append(histogram.Integral(k,k))
    vl=np.array(value)
    tgraph= ROOT.TGraph(len(en),en,vl)
    return tgraph

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
    print(integral,"from simpsons rule")

def fancygraph(g):
    g.GetXaxis().SetTitle("Energy (MeV)")
    g.GetYaxis().SetTitle("Flux neutrinos (s^-1*cm^-2)")
    g.GetXaxis().CenterTitle()
    g.GetYaxis().CenterTitle()
    g.Draw("ALP")
    return g

#save as root file
filesave = ROOT.TFile.Open("Boron8spectrum.root", "RECREATE")
def savefiles(graph,title):
    filesave.cd()
    graph.Write(title)

c1=ROOT.TCanvas("c1")
c1.SetLogy(1)

#chose the nominal or plus,muinus for hist, choose normalization, chose saved name in rootfile
def doAll(hist,norm,savename):
    r1=Scalehist(hist,norm,0,16.56)
    g1=histtograph(r1)
    IntegrateTF1(g1,0,16.56)
    g1=fancygraph(g1)
    savefiles(g1,savename)

doAll(h1,5.5138e6,"Nominal_low_B8") #for lower limit
doAll(h1,4.984e6,"Nominal_high_B8") #for upper limit
doAll(h1,5.244e6,"Nominal_B8") #this is nominal flux

doAll(h3,4.984e6,"spectrumlow_fluxlow_B8") #for spectrum shape -3sigma limit and flux lower
doAll(h2,5.5138e6,"spectrumhigh_fluxhigh_B8") #for spectrum shape + 3sigma and flux upper limit

doAll(h3,5.244e6,"spectrumlow_B8") #this is -3sigma flux
doAll(h2,5.244e6,"spectrumhigh_B8") #this is +3sigma flux 

