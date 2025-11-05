import os

NUMTREES = 1000
TREESPERFILE = 10
old = ["DD0???"]
new = [""]

submit_str = "#!/bin/bash\n\n"

for DD in range(NUMTREES):

    py_file = open("catalogs/catalog_base.txt", "r")
    replacement = ""

    for line in py_file:
        line = line.strip()
        changes = line.replace("/scratch/08288/tg875874/trees/tree_NUM", "/scratch/08288/tg875874/trees/tree_"+str(DD))
        replacement = replacement + changes + "\n"

    py_file.close()
    # opening the file in write mode
    fout = open("catalogs/catalog_"+str(DD)+".txt", "w")
    fout.write(replacement)
    fout.close()

filename = "sbatch_scripts/trees_"
for DD in range(int(NUMTREES/TREESPERFILE)):

    replacement = ""

    sbatch_file = open(filename+"base.sh", "r")
    for line in sbatch_file:
        line = line.strip()
        changes = line.replace("myjob.oNUM", "myjobs/myjob.o"+str(DD))
        changes = changes.replace("myjob.eNUM", "myjobs/myjob.e"+str(DD))
        for x in range(TREESPERFILE):
            num = str(DD*TREESPERFILE + x)
            changes = changes.replace("catalogNUM"+str(x)+".txt", "catalog_"+num+".txt")
            changes = changes.replace("seed"+str(x), num)
        replacement = replacement + changes + "\n"


    sbatch_file.close()

    fout = open(filename+str(DD)+".sh", "w")
    fout.write(replacement)
    fout.close()
    submit_str += "sbatch "+ filename+str(DD)+".sh\n\n"



submit_name = "submit_trees.sh"
sbatch_submit = open(submit_name, "w")
sbatch_submit.write(submit_str)
sbatch_submit.close()
os.system("chmod +x "+submit_name)



submit_str = "#!/bin/bash\n\n"

filename = "sbatch_scripts/find_dcbhs_"
for DD in range(TREESPERFILE):

    sbatch_file = open(filename+"base.sh", "r")
    replacement = ""
    for line in sbatch_file:
        line = line.strip()
        changes = line.replace("myjob.oNUM", "myjobs/myjob_find.o"+str(DD))
        changes = changes.replace("myjob.eNUM", "myjobs/myjob_find.e"+str(DD))
        changes = changes.replace("BLOCK", str(DD))
        replacement = replacement + changes + "\n"

    sbatch_file.close()

    fout = open(filename+str(DD)+".sh", "w")
    fout.write(replacement)
    fout.close()

    submit_str += "sbatch "+filename+str(DD)+".sh\n\n"




submit_name = "submit_find_dcbh_jobs.sh"
sbatch_submit = open(submit_name, "w")
sbatch_submit.write(submit_str)
sbatch_submit.close()
os.system("chmod +x "+submit_name)

