import urllib2
from bs4 import BeautifulSoup

def parse( page ) :
    lists = []
    cntr = 0 
    import re
    all_jobs = page.body.find_all(name="div",attrs={"class","jobwrapper"})
    for i in all_jobs :
        release = i.find_parent("div","tagwrapper").find_previous_sibling("a").contents[1].string.strip().rstrip(":")
        tag_set = i.find_parent("div","tagwrapper").find_previous_sibling("a").contents[2].string.strip()
        global_tag = i.find_parent("div","dsetwrapper").find_previous_sibling("a").contents[-1].split()[0].strip()
        job_id = i.attrs["id"].lstrip("job").strip()
        #status = i.find_parent("div","dsetwrapper").find("a").attrs["class"][0].strip()
        status = i.find_previous_sibling().attrs["class"][0].strip()
        user = None
        path = None
        if ( status != "Unclaimed" ) :
            user = i.find("br").contents[0].strip().split("@")[0] 
            path = i.previous_element.strip().replace("//","/").rstrip('/')
        data_sets = []
        for j in i.find_parent("div","dsetwrapper").find_previous_sibling("a").find("b").contents :
            data_sets = data_sets + filter( None, str(j).replace('</br>','<br>').split('<br>') )
        #print job_id,release,tag_set,global_tag,status,user,path,data_sets

        for k in data_sets :
            list = [job_id,release,tag_set,global_tag,status,user,path,k]
            lists.append(list)

    import numpy as np
    data = np.array(lists)
    return data

c = urllib2.urlopen('file:///Users/bainbrid/Repositories/susy_work/parse_icf/status4.html')
r = c.read()
page = BeautifulSoup(r)
data = parse(page)

#unique_tag_sets = set(zip(data[:,1],data[:,2]))

#print "Number of unique tag sets:",len(unique_tag_sets)
#print "List of unique tag sets:"
#print unique_tag_sets

#print set(data[:,1])

#for i in set(zip(data[:,1],data[:,2])) :
#    indices = data[:,1]==i[0] & data[:,2]==i[1]
#    print indices
#    print data[indices]
