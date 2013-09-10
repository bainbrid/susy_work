# Text contains output from command "du -h /pnfs/cms/WAX/11/store/user/lpcsusyra1/"

path = "/Users/bainbrid/Documents/WORK/susy/inventory/susy_work/utils/inventories/130429/"

input = open(path+'input.txt','r')
lines = input.readlines()

icf = [ line for line in lines if line.split()[1].count('ICF/automated') == 1 and line.split()[1][-19:].count('_') == 5 ]
output_icf = open('output_icf.txt','w')
output_icf.writelines(icf)
output_icf.close()

other = [ line for line in lines if \
              ( line.split()[1].count('ICF') == 0 and line.split()[1][39:].count('/') == 1 ) or \
              ( line.split()[1].count('automated') == 0 and line.split()[1].count('ICF') == 1 and line.split()[1][39:].count('/') == 2 ) ]
output_other = open('output_other.txt','w')
output_other.writelines(other)
output_other.close()

total = [ line for line in lines if len(line.split()[1][39:]) > 0 and line.split()[1][39:].count('/') == 0 ]
output_total = open('output_total.txt','w')
output_total.writelines(total)
output_total.close()

ra1 = ["henning","clucas","samr","yezhaq","zmeng","dburton"]
total_ra1 = [ line for line in total for user in ra1 if line.count(user) > 0 ]
output_ra1 = open('output_ra1.txt','w')
output_ra1.writelines(total_ra1)
output_ra1.close()

ra4 = ["agapitos","gouskas","nsaoulid","mstoye","karage","georgia"]
total_ra4 = [ line for line in total for user in ra4 if line.count(user) > 0 ]
output_ra4 = open('output_ra4.txt','w')
output_ra4.writelines(total_ra4)
output_ra4.close()


