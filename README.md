# threebodyparah2
The threebodyparah2 program will use the RKHS method with the grid of aug-cc-pVTZ data points (with empirical adjustments) to calculate the three-body correction energy to the parahydrogen triplet. Only the isotropic term in the potential expansion is used. The program accepts as parameters the three side lengths of the triangle (in Angstroms), and returns the three-body correction energy in cm^{-1}.

The input file takes the form of three columns of whitespace-delimited floating point numbers.
Each line corresponds to a different triangular configuration for the three hydrogen molecules.
The order of the three pair distances isn't important (i.e. the labels {1, 2, 3} for the molecules is arbitrary); the distances are sorted within the program.

The 1st column represents the pair distances between molecules 1 and 2.
The 2nd column represents the pair distances between molecules 2 and 3.
The 3rd column represents the pair distances between molecules 1 and 3.

-----------------------------------------------------------------------------------

Usage:
    Part 1: create the PES executable

        cd src
        make                                            # creates the `threebodyparah2` executable
        cd ..

    Part 2: create the .kernel file
    
        cd RKHS
        make
        python generate_kernel.py                       # creates the threebodypes.kernel file in /kernels
        cd ..

    Part 3: test the executable
        
        cp kernels/threebodypes.kernel test
        cp src/threebodyparah2 test
        cd test
        ./threebodyparah2 threebodypes.kernel tsl_equilateral.dat    # test using equilateral triangles
        ./threebodyparah2 threebodypes.kernel tsl_random.dat         # test using random triangles

        python plot_examples.py                                      # several examples of curves using the three-body potential
                                                                     # uncomment lines at the bottom of the script for different plots
                                                                     # more instructions are inside the `plot_examples.py` file


    Part 4 (OPTIONAL): create the threebodypes.csv file from the AVTZ data
        > this step assumes
            > the `threebodyparah2` executable in Part 1 has been created
            > `make` has already been run in Part 2
        > descriptions of the adjustments are given at the top of the files:
            > `make_smodAVTZ_csvfile.py`
            > `make_RsmodAVTZ_csvfile.py`
            > `make_threebodypes_csvfile.py`
        > the resulting `threebodypes.csv` file is saved in the `makecsv` directory, not the `csvfiles` directory

        cd makecsv
        python makecsv.py

-----------------------------------------------------------------------------------

TESTED COMPILERS:
    GNU Fortran, version 9.3.0        # gfortran
    g++, version 9.3.0                # g++

PYTHON REQUIREMENTS
    Python 3.6.7 or above
    numpy
    scipy

-----------------------------------------------------------------------------------

NOTES:
	(1) The	`threebodypes.kernel` file must be in the same directory as the `threebodyparah2` executable.
	(2) The program does not check if the triangles inserted are 'impossible', i.e. if you include a "triangle" where one side length is greater than the sum of the other two, a junk value comes out, or a segmentation fault occurs.
