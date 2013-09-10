ht_bin_options = {
    "200":"200_275_73_73_36",
    "275":"275_325_73_73_36",
    "325":"325_375_86_86_43",
    "375":"375_475_100_100_50",
    "475":"475_575_100_100_50",
    "575":"575_675_100_100_50",
    "675":"675_775_100_100_50",
    "775":"775_875_100_100_50",
    "875":"875_100_100_50",
    }

mht_met_options = {
    "with":(0.00,"_mht"),
    "without":(1.25,""),
    }

alpha_t_options = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.70]

njet_options = {
    "ge2j":"ge2j",
    "le3j":"le3j",
    "ge4j":"ge4j",
    }

def filename( option ) :
    return "text_HT"+ht_bin_options[option[0]]+\
        mht_met_options[option[1]][1]+\
        "_AlphaT_"+njet_options[option[2]]+".txt"

def alphat(val) :
    val = float(val)
    for i in range(len(alpha_t_options)-1) :
        if ( val >= alpha_t_options[i] ) and ( val < alpha_t_options[i+1] ) :
            return (str(alpha_t_options[i]),alpha_t_options[i])
    return (None,None)

def parse_files( path, files_options ) :
    effs = {}
    for option in files_options :
        name = filename(option) 
        print name
        file = open(path+name)
        for line in file.readlines() :
            entries = line.split()
            alpha_t = alphat(entries[2])
            mht_met = mht_met_options[option[1]][0]
            if ( alpha_t[0] == "0.60" and entries[9] == "Cumu" ) or entries[9] == "Differential" :
                effs[(str(mht_met),alpha_t[0])] = (entries[3],entries[5],entries[7])
    return effs

################################################################################
################################################################################
################################################################################

base="/Users/bainbrid/Repositories/susy_work/trigger/files/"
options = [
    ("200","with","ge2j"),
    ("200","without","ge2j"),
    ]
print parse_files(base,options)
