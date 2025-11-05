import os

TIME = "24:00:00"
NUMTREES = 1000
TREESPJOB = 20
NUMJOBS = int(NUMTREES/TREESPJOB)
print(NUMJOBS)

submit_str = "#!/bin/bash\n\n"

fname = "sbatch_scripts/find_leo_candidates_"
for DD in range(NUMJOBS):

    replacement = ""

    sbatch_file = open(fname+"base.sh", "r")
    START = TREESPJOB*DD
    STOP = START+TREESPJOB
    for line in sbatch_file:
        line = line.strip()
        changes = line.replace("myjob.oNUM", "myjobs/myleojob.o"+str(DD))
        changes = changes.replace("myjob.eNUM", "myjobs/myleojob.e"+str(DD))
        changes = changes.replace("outputNUM.txt", "outputs/output"+str(DD)+".txt")
        changes = changes.replace("START", str(START))
        changes = changes.replace("STOP",  str(STOP))
        changes = changes.replace("24:00:00", TIME)
        replacement = replacement + changes + "\n"

    sbatch_file.close()

    fout = open(fname+str(DD)+".sh", "w")
    fout.write(replacement)
    fout.close()
    submit_str += "sbatch " +fname+str(DD)+".sh\n\n"







submit_name = "submit_find_leos.sh"
sbatch_submit = open(submit_name, "w")
sbatch_submit.write(submit_str)
sbatch_submit.close()
os.system("chmod +x "+submit_name)

