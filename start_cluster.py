import os
from dask_jobqueue import PBSCluster
cluster = PBSCluster(env_extra=[f"export PYTHONPATH={os.path.abspath('Gillespie')}"])
cluster.scale(4)

