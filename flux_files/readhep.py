## this code reads in data for the Hep spectrum from http://www.sns.ias.edu/~jnb/SNdata/Export/Hepspectrum/hepspectrum.dat
# need to download the text file and put its location and name in dir and file locations
#The code then takes the data and normalizes the Hep flux to whatever value you would like and then saves the Tgraph object as a root file
#written by Trevor Hoyte

import ROOT,array
import numpy as np

normalization=9.3e3 # from SM200 from Bahacall https://iopscience.iop.org/article/10.1086/321493/pdf
dir="/project/6004969/trevorh/cross/flux_files/"
file1="hepspectrum.txt"

fluxfile = open(dir+file1, 'r')
lines = fluxfile.readlines()
nbins=1001 # number of data points+1
h1=ROOT.TH1D("h1","Hep",nbins,0.018,18.784) #change range to the range in the text file
c1=ROOT.TCanvas("c1")
c1.SetLogy(1)

energy=[]
#read through text files and extract info
for line in lines:
    data=line.strip().split()
    if len(data)==2:
        #fill hist and Tgraphs
        h1.Fill(float(data[0]),float(data[1]))
        energy.append(float(data[0]))
en=np.array(energy)

def Scalehist(histogram,scale,x1,x2):
    minbin=int(x1/0.019+1)
    maxbin=int(x2/0.019+1)
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
    n=1000
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
    c1.SaveAs("Hep.png")
    return g

#save as root file
filesave = ROOT.TFile.Open("Hepspectrum.root", "RECREATE")
def savefiles(graph,title):
    filesave.cd()
    graph.Write(title)


#scale histogram make it a Tgraph and draw + 
print("unnormalized flux: ", h1.Integral("weight"))
r1=Scalehist(h1,normalization,0.018,18.78)
#change hist to Tgraph object 
g1=histtograph(r1)
# check normalizations quick things using simpsons method
IntegrateTF1(g1,0.018,18.78)
#make Graphs pretty
g1=fancygraph(g1)
#save graphs in a root file
savefiles(g1,"Hep_flux")

