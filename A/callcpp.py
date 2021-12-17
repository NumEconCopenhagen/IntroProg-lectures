import numpy as np
import os
import ctypes as ct

def compile_cpp(filename,
                vs_path = 'C:/Program Files (x86)/Microsoft Visual Studio/2017/Community/VC/Auxiliary/Build/'):      
    """compile cpp file to dll
    
    Args:
    
        filename (str): path to .cpp file (no .cpp extensions!)
        vs_path (str,optional): path to vs compiler
    
    """

    # a. compile string
    pwd_str = 'cd "' + os.getcwd() + '"\n'    
    path_str = f'cd "{vs_path}"\n'
    version_str = 'call vcvarsall.bat x64\n'
    compile_str = f'cl /LD /EHsc /Ox /openmp {filename}.cpp\n'
    lines = [path_str,version_str,pwd_str,compile_str]

    # b. write .bat
    with open('compile.bat', 'w') as txtfile:
        txtfile.writelines(lines)
                               
    # c. compile
    result = os.system('compile.bat')
    if result == 0:
        print('cpp files compiled')
    else: 
        raise ValueError('cpp files can not be compiled')
        
    # d. clean
    os.remove(f'{filename}.obj')
    os.remove(f'{filename}.lib')
    os.remove(f'{filename}.exp') 
        
def set_argtypes(cppfile,funcs):
    """ set argument types
    
    Args:
        cppfile (ctypes.CDLL): c++ library (result of ct.cdll.LoadLibrary('cppfile.dll'))
        funcs (list): list of functions with elements (functionname,[argtype1,argtype2,etc.])
        
    """

    for func in funcs:
        name = func[0]
        argtypes = func[1]        
        funcnow = getattr(cppfile,name)
        funcnow.restype = None
        funcnow.argtypes = argtypes

def link_cpp(filename,funcs): 
    """ link cpp library
        
    Args:
        filename (str): path to .dll file (no .dll extension!)
        funcs (list): list of functions with elements (functionname,[argtype1,argtype2,etc.])
        
    Return:
        cppfile (ctypes.CDLL): c++ library (result of ct.cdll.LoadLibrary('cppfile.dll'))
    
    """

    # a. link
    cppfile = ct.cdll.LoadLibrary(filename + '.dll')
    print('cpp files loaded')
    
    # b. functions
    set_argtypes(cppfile,funcs)
    cppfile.setup_omp() # must exist
    delink_cpp(cppfile,filename,do_print=False,do_remove=False)
    cppfile = ct.cdll.LoadLibrary(filename)
    set_argtypes(cppfile,funcs)

    return cppfile

def delink_cpp(cppfile,filename,do_print=True,do_remove=True):
    """ delink cpp library
        
    Args:
    
        cppfile (ctypes.CDLL): c++ library (result of ct.cdll.LoadLibrary('cppfile.dll'))
        filename (str): path to .dll file (no .dll extension!).
        do_print (bool,optional): print if successfull    
        do_remove (bool,optional): remove dll file after delinking
        
    """

    # a. get handle
    handle = cppfile._handle

    # b. delete linking variable
    del cppfile

    # c. free handle
    ct.windll.kernel32.FreeLibrary.argtypes = [ct.wintypes.HMODULE]
    ct.windll.kernel32.FreeLibrary(handle)
    if do_print:
        print('cpp files delinked')

    # d. remove dll file
    if do_remove:
        os.remove(filename + '.dll')
