#kinit -l 1d bainbrid@CERN.CH
#cd ../../
#svn update notes/AN-13-000
#cd notes/AN-13-000/trunk/
eval `../../tdr runtime -sh` 
tdr --style an b AN-13-000
echo "tdr --style an b AN-13-000"
