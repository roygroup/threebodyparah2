import os
import sys
import pathlib
import subprocess

class KernelGeneratorInfo:
    def __init__(self, power1, power2, csvfile, kerneldir):
        # check if power1 and power2 are valid
        assert 0 <= power1 <= 6
        assert 0 <= power2 <= 6
        
        # setting up the input parameters
        self.executable = os.path.join("src", "makekernel.x")
        self.power1     = power1
        self.power2     = power2
        self.csvfile    = csvfile
        
        # creating the template for the file name 
        prefix              = pathlib.Path(csvfile).stem
        _kernelfiletemplate = "{0}.kernel".format(prefix)
        self.kerneltemplate = os.path.join(kerneldir, _kernelfiletemplate)
    
    def call_makekernel(self):
        # uses the command line to call the ./makekernel.x fortran code
        kernelfile = self.kerneltemplate
        cmd = "{} {} {} {} {}".format(self.executable, self.csvfile, kernelfile,
                                      self.power1, self.power2)
        
        print("Creating {}".format(kernelfile))
        subprocess.call(cmd, shell = True)

def make_kernel(prefix):
    csvdir    = os.path.join("..", "csvfiles")
    kerneldir = os.path.join("..", "kernels")
    os.mkdir(kerneldir)
    csvfile   = os.path.join(csvdir, prefix + ".csv")
    
    # first power doesn't actually matter, energies are empirically determined in that direction
    # second power should be 5
    kgen = KernelGeneratorInfo(3, 5, csvfile, kerneldir)
    kgen.call_makekernel()

if __name__ == '__main__':
    if len(sys.argv) == 1:
        make_kernel('threebodypes')
    else:
        kernelname = sys.argv[1]
        make_kernel(kernelname)
