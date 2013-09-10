# parse zhaoxia's 

base = 'http://zmeng.web.cern.ch/zmeng/SUSY2012/EightTeV2012/'

plots = [ 'B' + i + '_' + j + '_DataSetDataJetHT2012' 
          for i in [str(k) for k in range(0,3)] 
          for j in ['0To-1j','2To3j','4To15j'] ]

import ROOT 

mg = ROOT.TMultiGraph()

region = ['C1','C2','C3','C4','D']
alphat = [0.51,0.52,0.53,0.54,0.55,0.60]

for njet in ['0To-1j','2To3j','4To15j'] :

    x = []
    y = []
    xe = []
    ye = []

    for i in len(region) :
        plot = str(region[i])+'_ReverseMHToverMHT_'+njet+'_DataSetDataJetHT2012'

        url = base+'Json21Sep2012_11p7ifb_data533p2/MuonIso012/PeriodABC/QCDk_fitAllHTbinns/'+plot+'.txt'
        #print url

        import urllib
        f = urllib.urlopen(url)

        import StringIO
        txt = StringIO.StringIO(f.read())
        
        s = txt.read().split()
        print plot,"->",s[7],s[10],s[12],s[15],s[17],s[19],s[21]
        
        width = alphat[i+1] - alphat[i]
        centre = alphat[i] + width/2.

        x.append(centre)
        xe.append(width/2.)
        y.append(s[12])
        ye.append(s[15])

    

    g = ROOT.TGraph(len(x),x,y)

#g.SetTitle("cosine in x=[%.1f, %.1f]" % (x[0], x[-1]))
#g.GetXaxis().SetTitle("x")
#g.GetYaxis().SetTitle("y")
#g.Draw("AL")

