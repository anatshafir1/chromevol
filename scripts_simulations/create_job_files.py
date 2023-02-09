def create_job_format_for_all_jobs(path, job_name, memory, queue, ncpu, cmd):
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q "+ queue+ "\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N "+ job_name+"\n"
    text += "#PBS -e " + path + "/"+job_name+".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name +".OU"+ "\n"
    text += "#PBS -l select=ncpus="+ str(ncpu)+ ":mem="+ memory+"\n"
    text += "cd "+ path+"\n"
    text += cmd
    text+= "\n"

    return text
