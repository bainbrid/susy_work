#!/usr/bin/env python

import sys
import subprocess
from optparse import OptionParser

# command line options

parser = OptionParser()
parser.add_option("-b",
                  "--base",
                  action="store",
                  dest="_base",
                  default=False,
                  help="specify base path")
(options, args) = parser.parse_args()
if len(options._base) == 0 : sys.exit() 

# human readable file sizes

def sizeof_fmt(num):
    for x in [' bytes',' KB',' MB',' GB']:
        if num < 1024.0 and num > -1024.0:
            return "%3.2f%s" % (num, x)
        num /= 1024.0
    return "%3.2f%s" % (num, ' TB')

# method to allow "du" on dcache

def srmdu( base,path,depth=-1 ) :
    skip = 7
    depth = depth + 1
    pieces = []

    try:
        proc = subprocess.Popen(["srmls",base,path],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        out, err = proc.communicate()
    except:
        print "Exception! ", base, path
        return 0

    pieces = out.strip().split()

    if pieces.index("512") != skip-2 :
        print "Error: ", pieces[:pieces.index("512")]
        return 0
        
    for piece in pieces :
        if piece.count("Exception") > 0 :
            print "Error: ", piece
            return 0

    if len(pieces) <= skip :
        print "Error: ", pieces
        return 0
        
    du = 0
    for index in range(skip,len(pieces),2) :
        try:
            tmp = int(pieces[index])
        except ValueError:
            "ValueError!"
            print pieces
            continue
        if tmp != 512 :
            du = du + tmp
        else :
            temp = srmdu( base,pieces[index+1],depth )
            du = du + temp
            if depth < 2 :
                print sizeof_fmt( temp ), pieces[index+1]

    return du

# execute code

base = options._base.split("SFN=")[0] + "SFN="
path = options._base.split("SFN=")[1]
size = srmdu( base,path )
print "Total:", sizeof_fmt( size )
