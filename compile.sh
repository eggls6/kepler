#!/bin/bash

#F2PY
cd ./src/

f2py -c kind_m.f90 const_m.f90 trafo3_m.f90 orbdyn_m.f90 kepler_m.f90  -m kepler

wait

cp kepler.*.so ../lib/kepler.so
mv kepler.*.so ../notebooks/kepler.so
#in case that throws an error try this:

#LDFLAGS="$LDFLAGS -shared" python -m numpy.f2py -c kind_m.f90 const_m.f90 trafo3_m.f90 orbdyn2_m.f90 kepler_m.f90 -m kepler
