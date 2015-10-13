import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import cython_gsl



"""
"""

machine=os.environ['MACHINE']

if machine=='mylaptop':
    module = [
    Extension("CGAL-IO", ["cgalio.pyx", "../c/cgalio.c"], 
    libraries=cython_gsl.get_libraries()+["cuba"], 
    include_dirs=[ np.get_include(), cython_gsl.get_cython_include_dir(), 
		  "../c/"], 
    library_dirs=[cython_gsl.get_library_dir(), './', "../c/"]    ) , 
    ]

    

elif machine=='gwln':
    module = [
    Extension("cgalio", ["cgalio.pyx", "../c/cgalio.c"], 
    libraries=cython_gsl.get_libraries()+["cuba"], 
    include_dirs=[#"/home/wangxin/workspace/lib/mylib/script/clib", 
                  np.get_include(), cython_gsl.get_cython_include_dir(), 
		  "../c/"], 
    library_dirs=[cython_gsl.get_library_dir(), './', "../c/"]    ) , 
    ]


elif machine=='cita':
    module = [
    Extension("cgalio", ["cgalio.pyx", "../c/cgalio.c"], 
    libraries=cython_gsl.get_libraries(),
    include_dirs=[#"/home/wangxin/workspace/lib/mylib/script/clib", 
                  np.get_include(), cython_gsl.get_cython_include_dir(), 
		  "../c/"], 
    library_dirs=[cython_gsl.get_library_dir(), './', "../c/"]    ) , 
    ]


setup( name = 'Cgal Wrapper', cmdclass = {'build_ext': build_ext},
       ext_modules = module)

#-------
