import urllib2
c = urllib2.urlopen('file:///Users/bainbrid/Documents/WORK/susy/inventory/susy_work/utils/status4.html')
r = c.read()

from bs4 import BeautifulSoup
page = BeautifulSoup(r)

#print page.prettify()
#page.decompose()

import ordereddict
import pandas as pd

headers = ["jobid","release","tag","dataset","globaltag","path","status","user"]
jobids=[]
releases=[]
tags=[]
datasets=[] 
globaltags=[] 
paths=[]
statuses=[]
users=[]

alltags = page.body.findAllNext(name="div",attrs={"class","tagwrapper"})
for ii,i in enumerate(alltags) :

    release = i.findPrevious(name="a").contents[1].string.strip().rstrip(":") 
    tag = i.findPrevious(name="a").contents[2].string.strip() 

    alldatasets = i.findAllNext(name="div",attrs={"class","dsetwrapper"})
    for jj,j in enumerate(alldatasets) :

        releases.append(release)
        tags.append(tag)
        
        datasets.append( j.findPrevious(name="b").contents[0].strip() )
        globaltags.append( j.findPrevious(name="a").contents[2].strip().split()[0] )
        paths.append( j.contents[5].strip().replace("//","/") )
        jobids.append( j.findNext(name="div",attrs={"class","jobwrapper"}).attrs["id"].lstrip("job") )
        status = j.findNext(name="a").attrs["class"][0].strip()
        statuses.append(status)
        if ( status != "Unclaimed" ) :
            users.append( j.findNext(name="div",attrs={"class","jobwrapper"}).contents[2].strip().split("@")[0] )
        else :
            users.append(None)

s = {}
s[headers[0]] = pd.Series(jobids)
s[headers[1]] = pd.Series(releases)
s[headers[2]] = pd.Series(tags)
s[headers[3]] = pd.Series(datasets)
s[headers[4]] = pd.Series(globaltags)
s[headers[5]] = pd.Series(paths)
s[headers[6]] = pd.Series(statuses)
s[headers[7]] = pd.Series(users)

df = pd.DataFrame(s,columns=headers)
print df

