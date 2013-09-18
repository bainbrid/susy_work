from uncertainties import unumpy
import numpy as np

class Table :

    def __init__(self,verbose=False) :
        self._verbose = verbose
        self._alphat_bins = []
        self._mhtmet_bins = []
        self._file = open('qcd.tex','w')
        if self._verbose :
            if self._file.closed is False :
                print "Opened file '%s'"%self._file.name
            else :
                print "Problem opening file '%s'!"%self._file.name
        self._called_preamble = False
        self._called_postamble = False
        self._preamble()

    def __del__(self) :
        self._postamble()
        if not self._file.closed and self._verbose :
            self._file.close() 
            print "Closed file '%s'"%self._file.name

    def _preamble(self) :
        if self._called_preamble is True : return 
        font_size=""#\scriptstyle"
        s  = "\documentclass[portrait,a4paper]{article}\n"
        s += "\usepackage[total={19.0cm,27.7cm},top=1.0cm,left=1.0cm,includefoot]{geometry}\n"
        s += "\\begin{document}\n"
        s += "\\renewcommand{\\arraystretch}{1.3}\n"
        s += "\\newcommand\AlphaT{\ensuremath{\\alpha_{"+font_size+"\\textrm{T}}}}\n"
        s += "\\newcommand\HT{\ensuremath{H_{"+font_size+"\\textrm{T}}}}\n"
        s += "\\newcommand\ET{\ensuremath{E_{"+font_size+"\\textrm{T}}}}\n"
        s += "\\newcommand\MHT{\ensuremath{\HT^{"+font_size+"\\textrm{miss}}}}\n"
        s += "\\newcommand\MET{\ensuremath{\ET^{"+font_size+"\\textrm{miss}}}}\n"
        s += "\\newcommand\stat{\ensuremath{{"+font_size+"\,\\textrm{(stat)}\,}}}\n"
        s += "\\newcommand\syst{\ensuremath{{"+font_size+"\,\\textrm{(syst)}\,}}}\n"
        if False : s += "\\newcommand{\general}[5]{\ensuremath{#1^{+#2}_{-#3}\stat^{+#3}_{-#4}\syst}}\n"
        else : s += "\\newcommand{\\fixedpoint}[5]{\ensuremath{#1^{+#2}_{-#3}{}^{+#4}_{-#5}}}\n"
        s += "\\newcommand{\scientific}[6]{\ensuremath{(#2^{+#3}_{-#4}{}^{+#5}_{-#6})\cdot10^{#1}}}\n"
        s += "\\newcommand{\zero}{\entry{0}{0}{0}{0}{0}}\n"
        if not self._file.closed :
            self._file.write(s)
            self._called_preamble = True
            if self._verbose :
                print "Added preamble to file '%s'"%self._file.name

    def _postamble(self) :
        if self._called_postamble is True : return  
        s  = "\n"
        s += "\end{document}\n"
        if not self._file.closed :
            self._file.write(s)
            self._called_postamble = True
            if self._verbose :
                print "Added postamble to file '%s'"%self._file.name

    def _scientific(self,data,i,j,decimal_places=2) :
        num = ("%e"%data[2][j,i]).split('e')
        numbers = [float(num[0])]
        expo = int(num[1])
        for k in range(4) : 
            num = ("%e"%data[k+3][j,i]).split('e')
            numbers.append(float(num[0])*(10.**(float(num[1])-float(expo))) )
        s = "\scientific{%i}"%expo
        for k in numbers : s += ("{%."+str(decimal_places)+"f}")%k
        return s
        
    def _fixed_point(self,data,i,j,decimal_places=2) :
        s = "\\fixedpoint"
        for k in range(5) : s += ("{%."+str(decimal_places)+"f}")%data[k+2][j,i]
        return s

    def _general_format(self,data,i,j) :
        expo = int(("%e"%data[2][j,i]).split('e')[1])
        if expo < -2 or expo > 1 :
            return self._scientific(data,i,j) 
        else :
            return self._fixed_point(data,i,j) 

    def misc(self,text,flush=True) :
        s = "\n\centerline{\LARGE\\bf "+text+"}\n"
        if flush and not self._file.closed :
            self._file.write(s)
            if self._verbose :
                print "Added header to file '%s'"%self._file.name
        return s

    def header(self,text,flush=True) :
        s = "\n\centerline{\LARGE\\bf "+text+"}\n"
        if flush and not self._file.closed :
            self._file.write(s)
            if self._verbose :
                print "Added header to file '%s'"%self._file.name
        return s
        
    def newpage(self,header_text=None) :
        s  = "\n\\newpage\n"
        if header_text is not None : 
            s += self.header(header_text,False) 
        if not self._file.closed :
            self._file.write(s)
            if self._verbose :
                print "Added newpage to file '%s'"%self._file.name

    def alphat_bins(self,binning) :
        self._alphat_bins = binning

    def mhtmet_bins(self,binning) :
        self._mhtmet_bins = binning

    def add_table(self,data,caption="",format=None,cumu=False) :
        if len(self._mhtmet_bins) == 0 or len(self._alphat_bins) == 0 : 
            print "ERROR: No size to MHT/MET or AlphaT binning!"
            return
        s  = "\n"
        s += "\\begin{table}[h!]\n"
        if caption == "" : caption = "A test caption."
        s += "\centering\n"
        s += "\\scriptsize\n"
        #s += "\\footnotesize\n"
        s += "\caption{"+caption+"}\n"
        s += "\label{tab:test}\n"
        s += "\\begin{tabular}{" + "c"*(len(self._mhtmet_bins)+1) + "}\n"
        s += "\hline\n"
        s += "& \multicolumn{" + str(len(self._mhtmet_bins)) + "}{c}{\MHT/\MET} \\\\[0.1cm]\n"
        s += "\cline{2-" + str(len(self._mhtmet_bins)+1) + "}\n"
        s += "\AlphaT"
        for i in range(len(self._mhtmet_bins)) :
            if i == len(self._mhtmet_bins)-1 : s += " & $>$%.2f"%self._mhtmet_bins[i]
            else : s += " & %.2f--%.2f"%(self._mhtmet_bins[i],self._mhtmet_bins[i+1]) 
        s += " \\\\\n"
        s += "\hline\n"
        for ii in range(len(self._alphat_bins)) :
            i = len(self._alphat_bins) - ii - 1
            if cumu is True : s += "$>$%.2f"%self._alphat_bins[i] 
            elif i == len(self._alphat_bins)-1 : s += "%.2f"%self._alphat_bins[i] 
            else : s += "%.2f--%.2f"%(self._alphat_bins[i],self._alphat_bins[i+1]) 
            if format == None :
                for j in range(len(self._mhtmet_bins)) : 
                    if np.isnan(data[j,i].n) : s += " & nan "
                    elif np.isinf(data[j,i].n) : s += " & inf "
                    else : s += " & " + "${:.3L}$".format(data[j,i])
#            elif format == "scientific" :
#                for j in range(len(self._mhtmet_bins)) : s += " & " + self._scientific(data,i,j)
#            elif format == "fixed" :
#                for j in range(len(self._mhtmet_bins)) : s += " & " + self._fixed_point(data,i,j)
#            else :
#                for j in range(len(self._mhtmet_bins)) : s += " & " + self._general_format(data,i,j)
            s += " \\\\\n"
        s += "\hline\n"
        s += "\end{tabular}\n"
        s += "\end{table}\n"
        if not self._file.closed :
            self._file.write(s)
            if self._verbose :
                print "Added table to file '%s'"%self._file.name

    def add_template(self) :
        s  = "\n"
        s += "\\begin{table}[h]\n"
        s += "\caption{A test caption.}\n"
        s += "\label{tab:test}\n"
        s += "\centering\n"
        s += "\\footnotesize\n"
        s += "\\begin{tabular}{ccccc}\n"
        s += "\hline\n"
        s += "& \multicolumn{4}{c}{\MHT/\MET} \\\\[0.1cm]\n"
        s += "\cline{2-5}\n"
        s += "\AlphaT & 0--1.25 & 1.25--2.5 & 2.5--3.75 & 3.75--5 \\\\\n"
        s += "\hline\n"
        s += "$>$0.55 & \zero & \zero & \zero & \zero \\\\\n"
        s += "0.54--0.55 & \zero & \zero & \zero & \zero \\\\\n"
        s += "0.53--0.54 & \zero & \zero & \zero & \zero \\\\\n"
        s += "0.52--0.53 & \zero & \zero & \zero & \zero \\\\\n"
        s += "0.51--0.52 & \zero & \zero & \zero & \zero \\\\[0.1cm]\n"
        s += "\hline\n"
        s += "\end{tabular}\n"
        s += "\end{table}\n"
        if not self._file.closed :
            self._file.write(s)
            if self._verbose :
                print "Added template to file '%s'"%self._file.name
        
if __name__=="__main__":
    tab = Table()
    tab.add_table()
    del tab
