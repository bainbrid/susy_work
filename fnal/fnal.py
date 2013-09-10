def users( path = "" ) :
    data = []
    file = open(path+"users.txt")
    for line in file.readlines() :
        data.append(line.strip())
    return data

def user_files( path = "", user_list = [] ) :
    data = users(path)
    data = [ entry+".txt" for entry in data ]
    return data

def dirs( path = "", dir_list = [] ) :
    data = []
    if len(dir_list) == 0 : 
        dir_list = user_files(path,users(path))
    for entry in dir_list :
        file = open(path+entry)
        for line in file.readlines() :
            temp = line.strip().replace('//','/').rstrip('/')
            if len(temp) > 0 and temp[0] != "#" : 
                # ignore blank or commented lines
                data.append(temp) 
    return data

def wget( url = "http://www.hep.ph.ic.ac.uk/~rjb3/FNAL/", users_file = "users.txt", keep_files = [] ) :
    local = "./download/"
    def get(local,file) :
        from subprocess import Popen, PIPE
        p1 = Popen(['wget','-P',local,file], stdout=PIPE, stderr=PIPE)
        out1,err1 = p1.communicate()
        if ( p1.returncode != 0 ) :
            print "Return code: ", p1.returncode, " stderr: ", err1.strip()
            return 
    import os
    if os.path.isdir(local) : 
        import shutil
        import time
        shutil.move(local,str(local.rstrip('/')+'.'+time.strftime('%Y_%m_%d_%H_%M_%S')+'/'))
    os.mkdir(local) 
    get(local,url+users_file)
    users_files = user_files( local, users(local) )
    all_files = users_files + keep_files
    print all_files
    for i in all_files :
        get(local,url+i)
    return local

def print_all_dirs() :
    path = "/Users/bainbrid/Desktop/FNAL/"
    data = dirs(path)
    print data
    print "Number of directories in files:",len(data)
    print "Number of unique dirs in files:",len(set(data))

def list_duplicates( dir_list ) :
    temp = dict((i,dir_list.count(i)) for i in dir_list)
    import numpy as np
    keys = np.array(temp.keys())
    values = np.array(temp.values())
    indices = values > 1
    if ( len(keys[indices]) == 0 ) :
        print "No duplicates"
    else :
        for i in zip(keys[indices],values[indices]) :
            print i

def check_directories( path, dir_list = [] ) :
    if len(dir_list) == 0 :
        if os.path.isdir(path) :
            print "Path exists: ",path
        else :
            print "Path does not exist: ",path
    else :
        number = 0
        for file in dir_list :
            if not os.path.isdir(path+file) :
                number = number + 1
                print "Directory does not exist: "
        print "Number of missing directories:",number

def du( dirs ) :
    from subprocess import Popen, PIPE
    total = 0
    totals = {}
    #print "Total number of directories:",len(dirs)
    for idx,dir in enumerate(dirs) :
        p = Popen(['du','-s',dir], stdout=PIPE, stderr=PIPE)
        out,err = p.communicate()
	if ( p.returncode != 0 ) :
            print "Return code: ", p.returncode, " stderr: ", err.strip()
        else :
            num = int(out.strip().split()[0])
	    total = total + num
            totals[dir] = num
    totals["Total"] = total
    return total,totals

def write_dirs( dirs_list = [], name = "temp" ) :
    dirs_list.sort()
    f = open(name+'.txt','w')
    for i in dirs_list :
        f.write(i+'\n')
    f.close()

#def write_dirs( dirs_list = {}, name = "temp" ) :
#    f = open(name+'.txt','w')
#    for key,value in dirs_list.iteritems() :
#        f.write(key+'\t'+value+'\n')
#    f.close()

keep_list = ['KEEP_RA1_HCP.txt',
             'KEEP_RA1_PARKED.txt',
             'KEEP_RA1_EXTRA.txt',
             'KEEP_RA4.txt',
             'KEEP_RA6.txt',
             ]

keep_list = ['KEEP_RA1_HCP.txt','KEEP_RA1_EXTRA.txt','KEEP_RA4.txt','KEEP_RA6.txt',]

drop_list = ['KEEP_RA1_PARKED.txt']

path = wget(keep_files=keep_list)
#path = "./"
#path = "./download/"
#path = "/Users/bainbrid/Desktop/FNAL/"

#print_all_dirs()

#all = set(dirs(path))
#keep = set(dirs(path,keep_list))

#print "all",len(all)
#print "keep",len(keep)
#print "intersect",len(all&keep)
#print "union",len(all|keep)
#print "diff",len(all-keep)
#print "diff",len(keep-all)
#print keep-all

#quit()

#list_duplicates(dirs(path,keep_list))
    
#write_dirs(list(keep),'KEEP')
#write_dirs(list(all-keep),'DELETE')
#yossof = [i for i in list(all-keep) if i.count('yeshaq') > 0]
#write_dirs(yossof,'yossof')

#print "KEEP:"
#k,kk = du(keep)
#print "Total:",k
#print "DELETE:"
#d,dd = du(list(all-keep))
#print "Total:",d
#print "YOSSOF:"
#y = du(yossof)
#print "Total:",y
