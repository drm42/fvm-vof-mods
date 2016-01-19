if ${?LD_LIBRARY_PATH} then
  setenv LD_LIBRARY_PATH /home/drmoser/fvm/build-ubuntu/lib:${LD_LIBRARY_PATH};
else
  setenv LD_LIBRARY_PATH /home/drmoser/fvm/build-ubuntu/lib;
endif
setenv PYTHONPATH /home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin:/home/drmoser/fvm/build-ubuntu/lib:/home/drmoser/fvm/build-ubuntu/lib64/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/lib/python2.7/site-packages:/home/drmoser/fvm/build-ubuntu/bin
# Add include and lib directories to compiler paths
if ${?C_INCLUDE_PATH} then
  setenv C_INCLUDE_PATH /home/drmoser/fvm/build-ubuntu/include:${C_INCLUDE_PATH};
else
  setenv C_INCLUDE_PATH /home/drmoser/fvm/build-ubuntu/include;
endif
if ${?CPLUS_INCLUDE_PATH} then
  setenv CPLUS_INCLUDE_PATH /home/drmoser/fvm/build-ubuntu/include:${CPLUS_INCLUDE_PATH};
else
  setenv CPLUS_INCLUDE_PATH /home/drmoser/fvm/build-ubuntu/include;
endif
setenv PATH /home/drmoser/fvm/build-ubuntu/bin:$PATH

# Need this to recompile MPM in its directory.
setenv MEMOSA_CONFNAME ubuntu
