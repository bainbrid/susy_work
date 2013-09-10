import json
import ordereddict
import pandas as pd

masses = ['160','170','180','190']
cache = {}

for i in range( len(masses) ) :
    
    f = open('T2cc_200_'+masses[i]+'.log')
    data = json.load( f, object_pairs_hook = ordereddict.OrderedDict )
    tmp = data['Cut Flow'].items()[1]

    effs = []
    cuts = []
    
    try: 
        for j in range(100) :
            cut = tmp[0] 
            eff = tmp[1]['efficiency']
            tmp = tmp[1].items()[5]
            cuts.append( cut.strip() ) 
            effs.append( eff )
    except IndexError:
        pass
    
    cache[masses[i]] = pd.Series( effs, cuts )
    f.close()

#print pd.DataFrame( cache )

