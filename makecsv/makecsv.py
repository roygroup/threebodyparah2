import os
from subprocess import call

if __name__ == '__main__':
    print("making `threebodypes.csv`: it will be saved in this directory.")
    pwd = os.getcwd()

    # copy `threebodyparah2` here
    call(['cp', os.path.join('..', 'src', 'threebodyparah2'), '.'])
    # cd to RKHS
    os.chdir(os.path.join('..', 'RKHS'))
    # create the AVTZdata.kernel file
    call(['python', 'generate_kernel.py', 'AVTZdata'])
    # cd back to here
    os.chdir(pwd)

    # run `make_smodAVTZ_csvfile.py`
    print("applying long range s adjustments")
    call(['python', 'make_smodAVTZ_csvfile.py'])

    # run `make_RsmodAVTZ_csvfile.py`
    print("applying long range R adjustments")
    call(['python', 'make_RsmodAVTZ_csvfile.py'])

    # run `make_threebodypes_csvfile.py`
    print("applying short range R adjustments")
    call(['python', 'make_threebodypes_csvfile.py'])

    # do cleanup
    print("DONE!")
    call(['rm', 'threebodyparah2'])
    call(['rm', os.path.join('..', 'kernels', 'AVTZdata.kernel')])
    call(['rm', os.path.join('..', 'csvfiles', 'smodAVTZ.csv')])
    call(['rm', os.path.join('..', 'csvfiles', 'RsmodAVTZ.csv')])
