Please install enterprise, ptmcmcsampler, nestle, and corner if you haven't done so already.
For example, you can install them in conda like this:
$ conda create --name pta-gwa python=3.9
$ conda activate pta-gwa
$ conda install -c conda-forge enterprise-pulsar ptmcmcsampler nestle corner mpi4py 
$ pip install git+https://github.com/dfm/acor.git@master


