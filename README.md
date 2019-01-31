# Advanced HPC Coursework: LBM

This is the base coursework for the Advanced High Performance Computing class.
In this repository you will find:

* Source code in the `d2q9-bgk.c` file
* Results checking scripts are in the `check/` directory

## Compiling and running

To compile type `make`.
Editing the values for `CC` and `CFLAGS` in the Makefile can be used to enable different compiler options or use a different compiler.
These can also be passed on the command line:

    $ make CFLAGS="-O3 -fopenmp -DDEBUG"

Input parameters and obstacle files are all specified on the command line of the `d2q9-bgk` executable:

    $ ./d2q9-bgk <paramfile> <obstaclefile>

For example:

    $ ./d2q9-bgk input_256x256.params obstacles_256x256.dat

## Checking results

An automated result checking script is provided to validate your results.
The script is written in Python, so you will need to load the `languages/python-2.7.6` module before using it.
Running `make check` will check the output file (average velocities and final state) against some reference results.
By default, it should look something like this:

    $ make check
    python check/check.py --ref-av-vels-file=check/128x128.av_vels.dat --ref-final-state-file=check/128x128.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    Total difference in av_vels : 5.270812566515E-11
    Biggest difference (at step 1219) : 1.000241556248E-14
      1.595203170657E-02 vs. 1.595203170658E-02 = 6.3e-11%

    Total difference in final_state : 5.962977334129E-11
    Biggest difference (at coord (6,2)) : 1.000588500943E-14
      3.329122639178E-02 vs. 3.329122639179E-02 = 3e-11%

    Both tests passed!

This script takes both the reference results and the results to check (both average velocities and final state).
This is also specified in the makefile and can be changed like the other options:

    $ make check REF_AV_VELS_FILE=check/128x256.av_vels.dat REF_FINAL_STATE_FILE=check/128x256.final_state.dat
    python check/check.py --ref-av-vels-file=check/128x256.av_vels.dat --ref-final-state-file=check/128x256.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    ...

All the options for this script can be examined by passing the `--help` flag to it.

    $ python check/check.py --help
    usage: check.py [-h] [--tolerance TOLERANCE] --ref-av-vels-file
                    REF_AV_VELS_FILE --ref-final-state-file REF_FINAL_STATE_FILE
    ...


## Running on BlueCrystal Phase 3

You can use the provided job script to run an LBM job through the BCp3 queueing system:

    $ qsub job_submit_d2q9-bgk

You can monitor the progress of your jobs using:

    $ qstat -u $USER

There are [more instruction on how to use the queueing system in the getting started tutorial](https://github.com/UoB-HPC/hpc-course-getting-started/blob/master/3_Queueing_Systems.md).

When finished, the output from your job will be in a file called `d2q9-bgk.out`:

    $ less d2q9-bgk.out

If you wish to run a different set of input parameters, you should modify `job_submit_d2q9-bgk` with your chosen options.

## Checking submission content

Before handing in the coursework, you can use the `check_submission.sh` script to make sure that your code builds in a clean environment.
This will reduce the chances of the automarker failing to build or run your code.

To use the script, simply run it from the directory containing the files you intend to submit:

    $ /path/to/check_submission.sh

The script will:

1. Unload all the modules currently loaded.
2. Load your modules and environment variables specified in `env.sh`.
3. Use `make` to build your code and verify that an executable with the expected name is produced.

If the submission checking script prints any errors, you should try to address those before you hand in. 

Note that `check_submission.sh` does _not_ run your code, and so you _cannot_ verify that the results produced by your application validate just by running this script. You should check the correctness of your results separately, e.g. using `make check`.


## Serial output for sample inputs

This section shows running times for the provided code on a Phase 3 node.

- 128x128:
    ```
    $ ./d2q9-bgk input_128x128.params obstacles_128x128.dat
    ==done==
    Reynolds number:                9.751927375793E+00
    Elapsed time:                   58.832851 (s)
    Elapsed user CPU time:          58.837055 (s)
    Elapsed system CPU time:        0.004999 (s)
    ```

- 128x256:
    ```
    $ ./d2q9-bgk input_128x256.params obstacles_128x256.dat
    ==done==
    Reynolds number:                3.715003967285E+01
    Elapsed time:                   118.999340 (s)
    Elapsed user CPU time:          119.013907 (s)
    Elapsed system CPU time:        0.002999 (s)
    ```

- 256x256:
    ```
    $ ./d2q9-bgk input_256x256.params obstacles_256x256.dat
    ==done==
    Reynolds number:                1.005141162872E+01
    Elapsed time:                   477.089262 (s)
    Elapsed user CPU time:          477.164459 (s)
    Elapsed system CPU time:        0.005999 (s)
    ```

- 1024x1024:
    ```
    $ ./d2q9-bgk input_1024x1024.params obstacles_1024x1024.dat
    ==done==
    Reynolds number:                3.375851392746E+00
    Elapsed time:                   1965.561403 (s)
    Elapsed user CPU time:          1965.957129 (s)
    Elapsed system CPU time:        0.015997 (s)
    ```

## Visualisation

You can view the final state of the simulation by creating a .png image file using a provided Gnuplot script:

    $ gnuplot final_state.plt

