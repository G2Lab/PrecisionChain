import os
import pandas as pd
import numpy as np

# Read envirental variable
NTASKS = int(os.environ.get('NTASKS', 1))
JOB_ID = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))


dirname = os.path.dirname(__file__)
filename = os.path.abspath(os.path.join(dirname, '..', 'data/clinical/person.csv'))

def getPersons(filename=filename):
    persons = pd.read_csv(filename)
    return persons.person_id.astype(str).sort_values().tolist()

# slice python list into chunks and return index of particular chunk
def getChunk(jobs_num=NTASKS, job_id=JOB_ID):
    array = getPersons()
    print("NTASKS: ", jobs_num, "JOB_ID: ", job_id, "Length of array: ", len(array))
    return np.array_split(array, jobs_num)[job_id].tolist()