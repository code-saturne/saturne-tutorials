#!/usr/bin/python

# syrthes.py V5

import os
import os.path
import sys
import re
import shutil
import datetime
import fnmatch
import subprocess

from optparse import OptionParser
from time import gmtime, strftime

# Global environnement variables
s4home = None
arch = None
mpidir = None
#
syrthes_pp = None
syrthes_ppfonc = None
syrthes_ppfonc_tmc = None
syrthes_post = None
syrthes_ensight = None
syrthes_grouphist = None
syrthes_med = None
#
# *************************
# BEGIN OF CLASS DEFINITION
# *************************

#-------------------------------------------------------------------------------
# Definition of a class for managing SYRTHES case
#-------------------------------------------------------------------------------

class SyrthesCase:

    def __init__(self, name = None,
                 data_file = "syrthes_data.syd",
		 case_dir = None,
                 data_dir = None,
                 exec_dir = None,
                 src_dir = None,
                 n_procs = 1,
                 n_procs_ray = -1,
                 part_tool_name ="metis",
                 prepro = True,
                 post_mode = None,
                 debug = False,
                 logfile = None):
        """
        Initialize the structure for a SYRTHES computation
        """

        self.n_procs = n_procs           # Number of processes for syrthes computation
        self.n_procs_ray = n_procs_ray   # Number of processes for computing radiative transfer
        self.part_tool_name = part_tool_name # name of the partitionning tool
        self.prepro = prepro             # Do preprocessing
        self.parall = False              # Parallel run ?
        self.param = None                # Parameters defined in data_file (SyrParam object)
        self.post_mode = post_mode       # King of post-processing format (None, Ensight, MED)
        self.debug = debug

        self.c_mesh_dir = None         # Will be set while reading param file
        self.c_mesh_name = None        # Will be set while reading param file
        self.r_mesh_dir = None         # Will be set while reading param file
        self.r_mesh_name = None        # Will be set while reading param file

        self.resu_dir = None           # Will be set while reading param file

        #fluid1d
        self.f1d_mesh_dir = None         # Will be set while reading param file
        self.f1d_mesh_name = None        # Will be set while reading param file

        # Case directory
        if case_dir != None:
            self.case_dir = case_dir
        else:
            case_dir = os.path.dirname(data_file)
            if os.path.isabs(case_dir):
                self.case_dir = case_dir
            else:
                self.case_dir = os.path.abspath(case_dir)

        self.data_file = os.path.basename(data_file)

        # Execution directory
        if exec_dir == None or len(exec_dir) == 0:
            self.exec_dir = self.case_dir
        else:
            if os.path.isabs(exec_dir):
                self.exec_dir = exec_dir
            else:
                self.exec_dir = os.path.abspath(exec_dir)

        # Data directory (data file, input user data files, ...)
        if data_dir != None:
            if os.path.isabs(data_dir):
                self.data_dir = data_dir
            else:
                self.data_dir = os.path.abspath(data_dir)
        else:
            self.data_dir = self.case_dir

        # Sources directory (user subroutines, Makefile)
        if src_dir == None or len(src_dir) == 0:
            self.src_dir = self.case_dir
        else:
            if os.path.isabs(src_dir):
                self.src_dir = src_dir
            else:
                self.src_dir = os.path.abspath(src_dir)

        # Post-processing directory
        self.post_dir = os.path.join(self.exec_dir, 'POST')

        # Partionning directory
        if self.n_procs > 1:
            self.part_dir = os.path.join(self.exec_dir, 'PART')
        else:
            self.part_dir = None

        # SYRTHES Case name
        if name == None or len(name) == 0:
            name = os.path.basename(self.case_dir)
        self.name = name

        # Coupling parameters
        if self.n_procs > 1:
            self.parall = True

        # Update data_dile
        if self.n_procs > 1:
            self.update_datafile = True
        else:
            self.update_datafile = False

        # Output
        self.logfile = logfile     # Output filename
        self.echo = False          # Redirection or not of the output into a file
        if logfile != None:
            self.echo = True


#***************************************************************************************
    def set_exec_dir(self, exec_dir):
#***************************************************************************************

        self.exec_dir = exec_dir

        # Post-processing directory
        self.post_dir = os.path.join(self.exec_dir, 'POST')

        # Partionning directory
        if self.n_procs > 1:
            self.part_dir = os.path.join(self.exec_dir, 'PART')
        else:
            self.part_dir = None

        # Data directory (data file, input user data files, ...)
        if self.data_dir == None or len(self.data_dir) == 0:
            self.data_dir = self.exec_dir

        # Sources directory (user subroutines, Makefile)
        if self.src_dir == None or len(self.src_dir) == 0:
            self.src_dir = self.exec_dir


#***************************************************************************************
    def set_src_dir(self, src_dir, force = True):
#***************************************************************************************

        # Sources directory (user subroutines, Makefile)
        if self.src_dir == None or len(self.src_dir) == 0:
            self.src_dir = src_dir
        else:
            if force == True:
                self.src_dir = src_dir

#***************************************************************************************
    def set_data_dir(self, data_dir, force = True):
#***************************************************************************************

        # Data directory (data file, input user data files, ...)
        if self.data_dir == None or len(self.data_dir) == 0:
            self.data_dir = data_dir
        else:
            if force == True:
                self.data_dir = data_dir

#***************************************************************************************
    def set_n_procs(self, n_procs):
#***************************************************************************************

        self.n_procs = int(n_procs)

        if self.n_procs > 1:
            self.parall = True # Parallel run to do
            if self.part_dir == None or len(self.part_dir) == 0:
                self.part_dir = os.path.join(self.exec_dir, 'PART')

#***************************************************************************************
    def set_logfile(self, logfile_name):
#***************************************************************************************

        self.logfile = str(logfile_name)
        self.echo = True

#***************************************************************************************
    def read_data_file(self):
#***************************************************************************************

        self.param = read_syrthes_param(os.path.join(self.case_dir,
                                                     self.data_file))

        # Define self.c_mesh_dir and self.c_mesh_name
        self.c_mesh_name = os.path.basename(self.param.c_mesh_name)


        if os.path.dirname(self.param.c_mesh_name)[0:2] == "C:": # convertir format de chemin Windows en cygwin
            dirname_formatcygwin = os.path.dirname(self.param.c_mesh_name).replace('C:', '/cygdrive/c', 1)
            c_mesh_dir = os.path.join(self.case_dir, dirname_formatcygwin)
        else :
            c_mesh_dir = os.path.join(self.case_dir, os.path.dirname(self.param.c_mesh_name))

        if os.path.isabs(c_mesh_dir):
            self.c_mesh_dir = c_mesh_dir
        else:
            self.c_mesh_dir = os.path.abspath(c_mesh_dir)

        if self.c_mesh_dir != self.exec_dir:
            self.update_datafile = True
        if self.c_mesh_dir != c_mesh_dir:
            self.update_datafile = True

        # Define self.r_mesh_dir and self.r_mesh_name if needed
        if self.param.r_mesh_name != None and len(self.param.r_mesh_name) > 0:

            self.r_mesh_name = os.path.basename(self.param.r_mesh_name)



            if os.path.dirname(self.param.r_mesh_name)[0:2] == "C:" : # convertir format de chemin Windows en cygwin
                dirname_formatcygwin = os.path.dirname(self.param.r_mesh_name).replace('C:', '/cygdrive/c', 1)
                r_mesh_dir = os.path.join(self.case_dir, dirname_formatcygwin)
            else :
                r_mesh_dir = os.path.join(self.case_dir, os.path.dirname(self.param.r_mesh_name))


            if os.path.isabs(r_mesh_dir):
                self.r_mesh_dir = r_mesh_dir
            else:
                self.r_mesh_dir = os.path.abspath(r_mesh_dir)

            if self.r_mesh_dir != self.exec_dir:
                self.update_datafile = True
            if self.r_mesh_dir != r_mesh_dir:
                self.update_datafile = True

        # Define self.f1d_mesh_dir and self.f1d_mesh_name if needed
        if self.param.f1d_mesh_name != None and len(self.param.f1d_mesh_name) > 0:

            self.f1d_mesh_name = os.path.basename(self.param.f1d_mesh_name)



            if os.path.dirname(self.param.f1d_mesh_name)[0:2] == "C:" : # convertir format de chemin Windows en cygwin
                dirname_formatcygwin = os.path.dirname(self.param.f1d_mesh_name).replace('C:', '/cygdrive/c', 1)
                f1d_mesh_dir = os.path.join(self.case_dir, dirname_formatcygwin)
            else :
                f1d_mesh_dir = os.path.join(self.case_dir, os.path.dirname(self.param.f1d_mesh_name))


            if os.path.isabs(f1d_mesh_dir):
                self.f1d_mesh_dir = f1d_mesh_dir
            else:
                self.f1d_mesh_dir = os.path.abspath(f1d_mesh_dir)

            if self.f1d_mesh_dir != self.exec_dir:
                self.update_datafile = True
            if self.f1d_mesh_dir != f1d_mesh_dir:
                self.update_datafile = True
# isabelle
       # Define self.result_dir and self.result_name
        self.result_name = os.path.basename(self.param.result_name)


        if os.path.dirname(self.param.result_name)[0:2] == "C:": # convertir format de chemin Windows en cygwin
            dirname_formatcygwin = os.path.dirname(self.param.result_name).replace('C:', '/cygdrive/c', 1)
            result_dir = os.path.join(self.case_dir, dirname_formatcygwin)
        else :
            result_dir = os.path.join(self.case_dir, os.path.dirname(self.param.result_name))

        if os.path.isabs(result_dir):
            self.result_dir = result_dir
        else:
            self.result_dir = os.path.abspath(result_dir)

        if self.result_dir != self.exec_dir:
            self.update_datafile = True
        if self.result_dir != result_dir:
            self.update_datafile = True

        # In case of code coupling -> update data file
# isabelle
        if self.param.coupling == True:
            self.result_dir = self.exec_dir
            self.update_datafile = True

        # In case of coupling with fluid 1D
        if self.param.coupling_1D:
            self.part_dir = os.path.join(self.exec_dir, 'PART')


#***************************************************************************************
    def logfile_init(self):
#***************************************************************************************
        if self.echo == True:
            date=strftime("%a, %d %b %Y %H:%M:%S", gmtime())
            self.logfile = os.path.join(self.exec_dir, os.path.basename(self.logfile))
            os.system('echo ' + date + ' > ' + self.logfile)


#***************************************************************************************
    def clean(self):
#***************************************************************************************

        dir_files = os.listdir(os.getcwd())

        tmp1 = ['tmp.data', 'compile.log', 'listsyr', 'listing',
               'syrthes.log', 'syrthes', 'syr_user_fct.c', 'tmc_user_fct.c']
        tmp2 = fnmatch.filter(dir_files, '*~')
        tmp = tmp1 + tmp2

        for t in tmp:
            if os.access(os.path.join(self.exec_dir, t), 6) == 1:
                os.remove(os.path.join(self.exec_dir, t))


#***************************************************************************************
    def interpret_func(self, dest = None):
#***************************************************************************************

        # Check syrthes environment and set paths
        if s4home == None:
            retval = check_and_load_env()
            if retval != 0:
                err = '\n   Error during checking and loading environment.\n'
                sys.stderr.write(err)
                return 1

        cur_dir = os.getcwd()
        os.chdir(self.case_dir)

        if self.param == None:
            err = '\n  Try to use data from SYRTHES data file\n' \
                  '  but this file was not read yet !\n    -> Error !\n'
            sys.stderr.write(err)
            sys.exit('Stop.')

        if self.param.interpreted_func == True:
            sys.stdout.write(' \nBuilding user functions..\n')

            cmd = syrthes_ppfonc + ' -d ' + self.data_file
            if self.echo == True:
                ret = os.system(cmd + ' >> ' + self.logfile)
            else:
                ret = os.system(cmd)

            if dest == None:
                dest = self.src_dir

            if dest != self.case_dir:
                shutil.move(os.path.join(self.case_dir, 'syr_user_fct.c'),
                            os.path.join(dest, 'syr_user_fct.c'))

            sys.stdout.write(' \nBuilding user functions for the Thermal Monte Carlo solver..\n')
            cmd = syrthes_ppfonc_tmc + ' -d ' + self.data_file
            if self.echo == True:
                ret = os.system(cmd + ' >> ' + self.logfile)
            else:
                ret = os.system(cmd)

            if dest == None:
                dest = self.src_dir

            if dest != self.case_dir:
                shutil.move(os.path.join(self.case_dir, 'tmc_user_fct.c'),
                            os.path.join(dest, 'tmc_user_fct.c'))

        # Comme back to initial directory
        os.chdir(cur_dir)


#***************************************************************************************
    def prepare_run(self,
                    exec_srcdir = None,
                    compile_logname = None):
#***************************************************************************************
        """
        Create execution directory if needed
        Copy data files, source files into execution directory if needed
        exec_srcdir: source files directory in execution directory
        by default it is the same place (not true in cfd coupling mode)
        """

        retval = 0

        # Check syrthes environment and set paths
        if s4home == None:
            retval = check_and_load_env()
            if retval != 0:
                err = '\n   Error during checking and loading environment.\n'
                sys.stderr.write(err)
                return 1

        # Sanity checks
        retval = check_and_create_dir(self.exec_dir)
        if retval != 0:
            err = '\n   Error during check/create directory: ' + self.exec_dir + '\n'
            sys.stderr.write(err)
            return 1

        # Clean execution directory before computation
        self.clean()

        # Treament of interpreted functions
        self.interpret_func(dest = exec_srcdir)

        # Copy source files into exec_srcdir
        # ----------------------------------

        src_files = []
        if self.src_dir != exec_srcdir and exec_srcdir != None:

            srcdir_files = os.listdir(self.src_dir)
            src_files = (fnmatch.filter(srcdir_files, '*.c')
                         + fnmatch.filter(srcdir_files, 'Makefile'))

            if len(src_files) > 0:

                # Add header files to list so as not to forget to copy them
                src_files = src_files + fnmatch.filter(srcdir_files, '*.h')

                # Copy source files to execution directory
                if (self.exec_dir != exec_srcdir):
                    retval = check_and_create_dir(exec_srcdir)
                    if retval != 0:
                        err = '\n   Error during check/create directory: '
                        err += self.exec_dir + '\n'
                        sys.stderr.write(err)
                        return 1

                for f in src_files:
                    src_file = os.path.join(self.src_dir, f)
                    dest_file = os.path.join(exec_srcdir, f)
                    shutil.copy2(src_file, dest_file)

            self.src_dir = exec_srcdir

        # Compile and build syrthes
        # -------------------------

        retval = build_syrthes(parall = self.parall,
                               cfd_coupling = self.param.coupling,
                               fname = compile_logname,
                               debug = self.debug,
                               srcdir = self.src_dir,
                               destdir = self.exec_dir)
        if retval != 0:
            err = '\n   Error during building syrthes executable.\n'
            sys.stderr.write(err)
            return 2

        # Copy syrthes data file if needed
        # --------------------------------

        if self.case_dir != self.exec_dir:
            src = os.path.join(self.case_dir, self.data_file)
            dst = os.path.join(self.exec_dir, self.data_file)
            if not os.path.isfile(src) and not os.path.islink(src):
                err = src + ' is not a regular file or link.\n     -> Error !\n'
                sys.stderr.write(err)
                return 1
            shutil.copy2(src, dst)

        # Copy all user data files from data_dir to exec_dir
        # --------------------------------------------------

        if self.data_dir != self.exec_dir:

            # In case: data_dir = src_dir but data_dir != exec_dir
            lst = os.listdir(self.data_dir)
            exclude = ['tmp.data', 'syrthes', 'syrthes.py', self.c_mesh_name,
                       'syrthes_data.syd_example', 'listing', 'listsyr',
                       'syrthes.log', self.data_file]
            exclude.append(src_files)

            if self.r_mesh_name != None and len(self.r_mesh_name) > 0:
                exclude.append(self.r_mesh_name)

            if self.f1d_mesh_name != None and len(self.f1d_mesh_name) > 0:
                exclude.append(self.f1d_mesh_name)

            for e in exclude:
                if e in lst:
                    lst.remove(e)

            for f in lst:
                src = os.path.join(self.data_dir, f)
                dst = os.path.join(self.exec_dir, f)
                if not os.path.isdir(src):
                    shutil.copy2(src, dst)
                else:
                    shutil.copytree(src, dst)
                    
            self.data_dir = self.exec_dir

        return retval


#***************************************************************************************
    def _format(self,name,np):
#***************************************************************************************
        # insertion du numero de part00000 avant une extension .truc
        ext = name.split(".")[-1]
        rname = name.replace("."+ext,"")
        return rname + "_" + str(np).zfill(5) + "part00000"+ "." + ext

#***************************************************************************************
    def _format2(self,name,np,i):
#***************************************************************************************
        # insertion du numero de part0000n avant une extension .truc
        ext = name.split(".")[-1]
        rname = name.replace("."+ext,"")
        return rname + "_" + str(np).zfill(5) + "part" + str(i).zfill(5) + "." + ext

#***************************************************************************************
    def _format3(self,name,np):
#***************************************************************************************
        # insertion du numero de part00000 avant une extension _toto.truc
        extglob = name.split("_")[-1]
        rname = name.split("_")[-2]
        return rname + "_" + str(np).zfill(5) + "part00000"+ "_" + extglob

#***************************************************************************************
    def preprocessing(self):
#***************************************************************************************

        print("\n  ---------------------------")
        print("  Start SYRTHES preprocessing")
        print("  ---------------------------\n")

        os.chdir(self.exec_dir)

        # Check syrthes environment and set paths
        if s4home == None:
            retval = check_and_load_env()
            if retval != 0:
                err = '\n   Error during checking and loading environment.\n'
                sys.stderr.write(err)
                return 1

        retval = 0
        if self.n_procs > 1 or self.param.coupling_1D == True: # Parallel run => first, mesh partitionning

            retval = check_and_create_dir(self.part_dir)
            if retval != 0:
                return retval

            if self.param == None:
                err = '\n  Try to use data from SYRTHES data file\n' \
                      '  but this file was not read yet !\n    -> Error !\n'
                sys.stderr.write(err)
                sys.exit('Stop.')

            # Partionning
            # -----------

            # Modify mesh prefix to put partionned mesh into PART directory
            new_c_mesh_prefix = os.path.basename(self.param.c_mesh_prefix)
            new_c_mesh_prefix = os.path.join('PART', new_c_mesh_prefix)

            if self.prepro == True:

                sys.stdout.write('Pre-processing SYRTHES files.. \n')

                #isabelle cmd = syrthes_pp + ' -v '+' -n '+ str(self.n_procs)
                cmd = syrthes_pp +' -n '+ str(self.n_procs)
                cmd += ' --toolpart ' + self.part_tool_name
                cmd += ' -m '+ '"' + os.path.join(self.c_mesh_dir, self.c_mesh_name) + '"'
                cmd += ' -d ' + '"' + self.data_file + '"' + ' -o ' + '"'+new_c_mesh_prefix+'"'

                #Damien
                if self.param.coupling_1D == True:
                    cmd += ' -f ' + '"' + os.path.join(self.f1d_mesh_dir, self.f1d_mesh_name) + '"' + ' -F ' + '"'+os.path.join(self.part_dir,self.f1d_mesh_name)+'"'
                print(cmd)
                sys.stdout.write('  --> ' + cmd +'\n')

                if self.echo:
                    retval = os.system(cmd + ' >> '+ self.logfile)
                else:
                    retval = os.system(cmd)

                if retval != 0:
                    return retval
                else:
                    sys.stdout.write('   -> OK\n\n')

        # Update mesh name in data file
        if self.update_datafile == True:
            sys.stdout.write('Updating the mesh file name.. \n')
            try:
                if self.data_file != None or len(self.data_file) > 0:

                    if (self.n_procs > 1 and self.prepro) or (self.param.coupling_1D and self.prepro) == True:
                        new_abs_c_mesh_name = os.path.join(self.part_dir, self._format(self.c_mesh_name,self.n_procs))
                    else:
                        new_abs_c_mesh_name = os.path.join(self.c_mesh_dir, self.c_mesh_name)

                    if self.r_mesh_dir != None and self.r_mesh_name != None:
                        if self.n_procs_ray > 1 and self.prepro == True:
                            new_abs_r_mesh_name = os.path.join(self.part_dir, self._format(self.r_mesh_name,self.n_procs_ray))
                        else:
                            new_abs_r_mesh_name = os.path.join(self.r_mesh_dir, self.r_mesh_name)
                    else:
                        new_abs_r_mesh_name = None

                    if self.f1d_mesh_dir != None and self.f1d_mesh_name != None:
                        if self.param.coupling_1D and self.prepro:
                            new_abs_f1d_mesh_name = os.path.join(self.part_dir, self._format(self.f1d_mesh_name,self.n_procs))
                        else:
                            new_abs_f1d_mesh_name = os.path.join(self.f1d_mesh_dir, self.f1d_mesh_name)
                    else:
                        new_abs_f1d_mesh_name = None

                    if self.param.restart == True:
                        new_abs_restart = os.path.join(self.case_dir, self.param.restart_file)
                    else:
                        new_abs_restart = None

                    retval = self.param.update(self.n_procs,
                                               new_abs_c_mesh_name,
                                               new_abs_r_mesh_name,
                                               new_abs_f1d_mesh_name,
                                               new_abs_restart,
                                               new_filename = "tmp.data")
                    self.data_file = "tmp.data"

                else:
                    sys.stderr.write('\n  Data file not defined.\n')
                    sys.stderr.write('    Error while updating the mesh file name !\n')
                    retval = 1
            except:
                sys.stderr.write('\n  Error while updating the mesh file name !\n')
                retval = 1
            else:
                sys.stdout.write('   -> OK\n\n')

        return retval


#***************************************************************************************
    def build_cmdline(self):
#***************************************************************************************

        cmd = ''
        mpi_env=[None, None, None] # Start script, Run script, Finalize script

        if self.n_procs > 1: # Parallel run

            # Test if mpi library is MPICH2 or OPENMPI (otherwise exit with error)
            # MPI launcher
            mpi_bindir=os.path.join(mpidir, 'bin')

            # Update mpi_env with start script
            mpistart=os.path.join(mpi_bindir,'mpdboot')
            if os.path.isfile(mpistart) or os.path.islink(mpistart):
                mpi_env[0] = mpistart

            # Update mpi_env with run script
            mpiexec=os.path.join(mpi_bindir,'mpiexec')
            if os.path.isfile(mpiexec) or os.path.islink(mpiexec):
                #mpi_env[1] = mpiexec
                mpi_env[1] = 'mpiexec'   # isa : on met en dur mpiexec a la place de /usr/bin/mpiexec
            else:
                err='  ' + str(mpiexec)+ ' is not a valid file\n'
                sys.stderr.write(err)
                sys.exit('Stop SYRTHES execution')

            # Update mpi_env with finalize script
            mpiexit=os.path.join(mpi_bindir,'mpdallexit')
            if mpi_env[0] != None:
                if os.path.isfile(mpiexit) or os.path.islink(mpiexit):
                    mpi_env[2] = mpiexit
                else:
                    err= 'MPICH2 implementation found but not mpdallexit\n'
                    err+= 'mpipath used: ' + mpi_bindir
                    sys.stderr.write(err)
                    sys.exit('Stop SYRTHES execution')
            else:
                mpi_env[2] = None

            if mpi_env[0] != None:
                cmd += mpi_env[0] + ' && '

            cmd += mpi_env[1] + ' -n ' + str(self.n_procs)

        cmd += ' ./syrthes -d ' + self.data_file

        if self.n_procs_ray > 0:
            cmd += ' -r ' + str(self.n_procs_ray)

        if self.logfile != None:
            cmd += ' --log ' + str(self.logfile)

        if mpi_env[2] != None:
            cmd += ' && ' + mpi_env[2]

        return cmd

#***************************************************************************************
    def run(self):
#***************************************************************************************

        if not os.path.isdir(self.exec_dir):
            err = '\n  Execution directory is not a directory!\n' \
                  '     -> Error !\n'
            sys.stdout.write(err)
            return 1

        os.chdir(self.exec_dir)

        if self.param.coupling == True:
            err = '\n Cannot run syrthes in code coupling mode with this function.\n'
            sys.stderr.write(err)
            return 1

        # Build command line
        msg = 'Execution of SYRTHES.. \n'
        msg += '    -> number of processors for conduction = ' + str(self.n_procs) + ' \n'
        if (self.n_procs_ray > 0):
            msg += '    -> number of processors for radiation  = ' + str(self.n_procs_ray)
            msg += ' (if needed only: see SYRTHES data file) \n'
        sys.stdout.write(msg)

        if (self.n_procs_ray > 0 and self.n_procs_ray > self.n_procs):
            err='   Number of proc. for radiation > Number of proc. for conduction\n'
            err+='  This is not possible.\n'
            err+='  Stop SYRTHES execution'
            sys.stderr.write(err)
            return 1

        cmd = self.build_cmdline()

        # Execution of syrthes
        try:
            if self.echo == True:
                ret = os.system(cmd +' >> ' + self.logfile)
            else:
                ret = os.system(cmd)
        except:
            ret = 1

        return ret

#***************************************************************************************
    def postprocessing_merge(self, mesh_path, resu_path0, merge_path):
#***************************************************************************************
        retval=0
        
        cmd = syrthes_post + ' -n ' + str(self.n_procs) + ' -m '+ '"'+ mesh_path +'"'
        cmd += ' -r ' + '"'+resu_path0+'"' + ' -o ' + '"'+merge_path+'"'

        if (os.path.isfile(mesh_path) and os.path.isfile(resu_path0)):
            try:
                sys.stdout.write('\n  .merging results '+ resu_path0 + ' --> _all.* \n')
                sys.stdout.write('     --> ' + cmd +'\n')
                if self.echo:
                    retval = os.system(cmd + ' >> ' + self.logfile)
                else:
                    retval = os.system(cmd)
            except:
                err = '   Cannot execute tool to merge results files '+ resu_path0 + ' for post-processing\n'
                err += '   -> Error !\n'
                sys.stderr.write(err)
                retval=1

        return retval



#****************************************************************************************
    def postprocessing_ens_med(self, mode, mesh_path, resu_path, output):
#****************************************************************************************

        retval=0

        if mode== 'ens':
            # Execute syrthes -> ensight tool
            sys.stdout.write('  Generation du fichier Ensight/Paraview... \n')
            cmd = syrthes_ensight + ' -m ' + '"'+mesh_path+'"' + ' -r ' + '"'+resu_path+'"' + ' -o ' + '"'+output+'"'
        elif mode == 'med':
            # Execute syrthes -> med tool
            sys.stdout.write('  Generation du fichier MED... \n')
            cmd = syrthes_med + ' -m ' + '"'+mesh_path+'"' + ' -r ' + '"'+resu_path+'"' + ' -o ' + '"'+output + '.med"'
            
        sys.stdout.write('     --> ' + cmd +'\n')
            
        if (os.path.isfile(mesh_path) and os.path.isfile(resu_path)):
            try:
                os.chdir(self.post_dir)
                if self.echo:
                    retval = os.system(cmd + ' >> ' + self.logfile)
                else:
                    retval = os.system(cmd)
                    os.chdir(self.exec_dir)
            except:
                err = '   Cannot execute tool to convert boundary results ' + resu_path + ' into post-processing format\n'
                err += '   ->' + resu_path + '\n'
                err += '   -> Error !\n'
                sys.stdout.write(err)
                retval=1

        return retval

            
        

    
#***************************************************************************************
    def postprocessing(self, mode = 'ens'):
#***************************************************************************************

        sys.stdout.write('Post-processing.. \n')

        # Check syrthes environment and set paths
        if s4home == None:
            retval = check_and_load_env()
            if retval != 0:
                err = '\n   Error during checking and loading environment.\n'
                sys.stderr.write(err)
                return 1

        retval = check_and_create_dir(self.post_dir)
        if retval != 0:
            return retval

        retval1=0      # merge des resu en parallele
        retval2=0      # merge des resu transitoires (rdt)  en parallele
        retval3=0      # merge des resu bord en parallele
        retval4=0      # merge des resu bord transitoire  en parallele
        retval5=0      # merge des resu cplcfd  en parallele
        retval6=0      # merge des resu cplcfd transitoire  en parallele

        # Parallel run --------------------------------------------------------------------
        if self.n_procs > 1:

            if not os.path.isdir(self.part_dir):
                err = '\n  Partition directory is not a directory!\n' \
                      '     -> Error !\n'
                sys.stdout.write(err)
                return 1


            
            # Maillage volumique
            # ------------------
            if self.prepro:
                mesh_path = os.path.join(self.part_dir,  self.c_mesh_name)
            else:
                mesh_path = os.path.join(self.c_mesh_dir,  self.c_mesh_name)
 
            # resultat final
            resu_path  = os.path.join(self.part_dir, self.result_name + '.res')
            resu_path0 = os.path.join(self.part_dir, self._format(self.result_name+'.res',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_all')
            retval1=self.postprocessing_merge(mesh_path, resu_path0, merge_path)
            
            # resultats transitoires
            resu_path = os.path.join(self.part_dir, self.result_name + '.rdt')
            resu_path0 = os.path.join(self.part_dir, self._format(self.result_name+'.rdt',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_all')
            retval2=self.postprocessing_merge(mesh_path, resu_path0, merge_path)


            
            # maillage pour les flux de bord
            # ------------------------------
            mesh_path = os.path.join(self.part_dir, self._format3(self.result_name+'_bord.syr',self.n_procs))
           
            # resultat final
            resu_path = os.path.join(self.part_dir, self.result_name + '_bord.res')
            resu_path0 = os.path.join(self.part_dir, self._format3(self.result_name+'_bord.res',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_bord_all')
            retval3=self.postprocessing_merge(mesh_path, resu_path0, merge_path)

            # resultats transitoires
            resu_path = os.path.join(self.part_dir, self.result_name + '_bord.rdt')
            resu_path0 = os.path.join(self.part_dir, self._format3(self.result_name+'_bord.rdt',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_bord_all')
            retval4=self.postprocessing_merge(mesh_path, resu_path0, merge_path)

                
            # maillage pour le couplage CFD
            # ------------------------------
            mesh_path = os.path.join(self.part_dir, self._format3(self.result_name+'_cplcfd.syr',self.n_procs))
           
            # resultat final
            resu_path = os.path.join(self.part_dir, self.result_name + '_cplcfd.res')
            resu_path0 = os.path.join(self.part_dir, self._format3(self.result_name+'_cplcfd.res',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_cplcfd_all')
            retval5=self.postprocessing_merge(mesh_path, resu_path0, merge_path)

            # resultats transitoires
            resu_path = os.path.join(self.part_dir, self.result_name + '_cplcfd.rdt')
            resu_path0 = os.path.join(self.part_dir, self._format3(self.result_name+'_cplcfd.rdt',self.n_procs))
            merge_path = os.path.join(self.result_dir, self.result_name + '_cplcfd_all')
            retval6=self.postprocessing_merge(mesh_path, resu_path0, merge_path)

                   
            # si au moins les fichiers resultats ont ete merges on continue le postprocessing
            if (retval1==0):
                mesh_path = os.path.join(self.result_dir, self.result_name + '_all.syr')
                bord_mesh_path = os.path.join(self.result_dir, self.result_name + '_bord_all.syr')
                cplcfd_mesh_path = os.path.join(self.result_dir, self.result_name + '_cplcfd_all.syr')
                
                resu_path = os.path.join(self.result_dir, self.result_name + '_all.res')
                rdt_path = os.path.join(self.result_dir, self.result_name + '_all.rdt')
                bord_resu_path = os.path.join(self.result_dir, self.result_name + '_bord_all.res')
                bord_rdt_path = os.path.join(self.result_dir, self.result_name + '_bord_all.rdt')
                cplcfd_resu_path = os.path.join(self.result_dir, self.result_name + '_cplcfd_all.res')
                cplcfd_rdt_path = os.path.join(self.result_dir, self.result_name + '_cplcfd_all.rdt')
                
                output = self.result_name + '_all'
                outputrdt = self.result_name + '_rdt_all'
                bord_output = self.result_name + '_bord_all'
                bord_outputrdt = self.result_name + '_bord_rdt_all'
                cplcfd_output = self.result_name + '_cplcfd_all'
                cplcfd_outputrdt = self.result_name + '_cplcfd_rdt_all'
                
                if self.param.coupling_1D :
                    f1d_mesh_path = os.path.join(self.part_dir, self._format(self.f1d_mesh_name,self.n_procs))
                    f1d_resu_path = os.path.join(self.result_dir, self.result_name + '_f1d.res')
                    f1d_output = self.result_name + '_f1d'
            else:
                # erreur de merge sur les fichiers .res --> stop, impossible de passer vers ensight
                return retval1

        # Serial run --------------------------------------------------------------------
        else:
            mesh_path = os.path.join(self.c_mesh_dir, self.c_mesh_name)
            bord_mesh_path = os.path.join(self.result_dir, self.result_name + '_bord.syr')
            cplcfd_mesh_path = os.path.join(self.result_dir, self.result_name + '_cplcfd.syr')

            resu_path = os.path.join(self.result_dir, self.result_name + '.res')
            rdt_path = os.path.join(self.result_dir, self.result_name + '.rdt')
            bord_resu_path = os.path.join(self.result_dir, self.result_name + '_bord.res')
            bord_rdt_path = os.path.join(self.result_dir, self.result_name + '_bord.rdt')
            cplcfd_resu_path = os.path.join(self.result_dir, self.result_name + '_cplcfd.res')
            cplcfd_rdt_path = os.path.join(self.result_dir, self.result_name + '_cplcfd.rdt')
            
            output = self.result_name
            outputrdt = self.result_name + '_rdt'
            bord_output = self.result_name + '_bord'
            bord_outputrdt = self.result_name + '_bord_rdt'
            cplcfd_output = self.result_name + '_cplcfd'
            cplcfd_outputrdt = self.result_name + '_cplcfd_rdt'
            
            if self.param.coupling_1D :
                f1d_mesh_path = os.path.join(self.f1d_mesh_dir, self.f1d_mesh_name)
                f1d_resu_path = os.path.join(self.result_dir, self.result_name + '_f1d.res')
                f1d_output = self.result_name + '_f1d'

        # --------------------------------------------------------------------

        # global result file to postprocessing file
        self.postprocessing_ens_med( mode, mesh_path, resu_path, output)
        
        if self.param.coupling_1D:
            self.postprocessing_ens_med( mode, f1d_mesh_path, f1d_resu_path, f1doutput)


        # rdt file
        if (retval2==0 and os.path.isfile(rdt_path)):
            self.postprocessing_ens_med( mode, mesh_path, rdt_path, outputrdt)

        # bord_resu file
        if (retval3==0 and os.path.isfile(bord_resu_path)):
            self.postprocessing_ens_med( mode, bord_mesh_path, bord_resu_path, bord_output)

        # bord_rdt file
        if (retval4==0 and os.path.isfile(bord_rdt_path)):
             self.postprocessing_ens_med( mode, bord_mesh_path, bord_rdt_path, bord_outputrdt)

        # cplcfd_resu file
        if (retval5==0 and os.path.isfile(cplcfd_resu_path)):
            self.postprocessing_ens_med( mode, cplcfd_mesh_path, cplcfd_resu_path, cplcfd_output)

        # cplcfd_rdt file
        if (retval6==0 and os.path.isfile(cplcfd_rdt_path)):
             self.postprocessing_ens_med( mode, cplcfd_mesh_path, cplcfd_rdt_path, cplcfd_outputrdt)


        return retval

#***************************************************************************************
    def save_results(self, save_dir = None, horodat = True, overwrite = False):
#***************************************************************************************

        if save_dir == None or len(save_dir) == 0:
            save_dir = os.path.join(os.getcwd(), "SAVE.SYR")
        if horodat == True:
            now = datetime.datetime.now()
            save_dir += "." + now.strftime('%m%d%H%M')

        if os.path.isdir(save_dir):
            if overwrite == False:
                err = "   Directory to store results is already existing."
                err+=save_dir
                err+="   Do not overwrite existing data !\n"
                sys.stderr.write(err)
                return 1
            else:
                shutil.rmtree(save_dir)

        try:
            os.mkdir(save_dir)
        except:
            err = '   Cannot create directory to save SYRTHES result.\n'
            err += '   -> Stop this step !\n'
            sys.stderr.write(err)
            return 1

        # Copy files
        ls = os.listdir(self.exec_dir)
        exclude = ['tmp.data', 'syrthes', 'syrthes.py', 'Makefile', 'POST', 'PART']
        exclude += fnmatch.filter(ls, '*~')
        for e in exclude:
            if e in ls:
                ls.remove(e)

        for l in ls:
            src = os.path.join(self.exec_dir, l)
            if os.path.isfile(src) and not os.path.islink(src):
                dst = os.path.join(save_dir, l)
                shutil.copy2(src, dst)
            if os.path.isdir(src) and not os.path.islink(src):
                dst = os.path.join(save_dir, l)
                shutil.copytree(src, dst, symlinks = False)

        if self.post_dir != None:
            try:
                shutil.copytree(self.post_dir, os.path.join(save_dir, 'POST'),
                                symlinks = False)
            except:
                pass

        if self.part_dir != None:
            try:
                shutil.copytree(self.part_dir, os.path.join(save_dir, 'PART'),
                                symlinks = False)
            except:
                pass

        # Sucessful execution
        return 0

#***************************************************************************************
    def dump(self):
#***************************************************************************************

        print( "\n  SYRTHES Case summary:\n")
        print( "    Name =                        ", self.name)
        print( "    Data file =                   ", self.data_file)
        print( "    Update Data file =            ", self.update_datafile)
        print( "    Do preprocessing =            ", str(self.prepro))
        print( "    Debug =                       ", str(self.debug))
        print( "    Case dir. =                   ", self.case_dir)
        print( "    Execution dir. =              ", self.exec_dir)
        print( "    Data dir. =                   ", self.data_dir)
        print( "    Source dir. =                 ", self.src_dir)
        print( "    Post dir. =                   ", self.post_dir)

        if self.part_dir != None:
            print( "    Part dir. =                   ", self.part_dir)

        print( "\n    Conduction mesh dir. =        ", self.c_mesh_dir)
        print( "    Conduction mesh name =        ", self.c_mesh_name)

        if self.r_mesh_name != None and len(self.r_mesh_name) > 0:
            print( "    Radiative mesh dir. =         ", self.r_mesh_dir)
            print( "    Radiative mesh name =         ", self.r_mesh_name)

        if self.f1d_mesh_name != None and len(self.f1d_mesh_name) > 0:
            print( "    1D fluid mesh dir. =         ", self.f1d_mesh_dir)
            print( "    1D fluid mesh name =         ", self.f1d_mesh_name)

        print( "\n    Total num. of processes =     ", str(self.n_procs))

        if self.n_procs_ray > 0:
            print( "    Num. of processes for radiative transfer =     ", str(self.n_procs_ray))

        if self.echo == True:
            print( "    Logfile name            =     ", self.logfile)

        print( "    Echo =                        ", str(self.echo))
        print( "    Parallel run =                ", str(self.parall))
        print( "    Do preprocessing =            ", str(self.prepro))

        if self.param != None:
            self.param.dump()

#-------------------------------------------------------------------------------
# Definition of a class for managing SYRTHES parameters
#-------------------------------------------------------------------------------

class SyrthesParam:

    def __init__(self, name = "syrthes.data"):
        """
        Initialize the structure for a SYRTHES computation
        """

        self.name = os.path.basename(name) # parameter filename
        self.c_mesh_name = None            # conduction mesh name
                                           # (only *.syr files are accepted)
        self.c_mesh_prefix = None
        self.c_mesh_suffix = None

        self.r_mesh_name = None
        self.r_mesh_prefix = None
        self.r_mesh_suffix = None
# isabelle

        #fluid 1D
        self.f1d_mesh_name = None
        self.f1d_mesh_prefix = None
        self.f1d_mesh_suffix = None

        self.result_name = None
        self.result_file = None

        self.restart_file = ''
        self.restart = False
        self.coupling = False
        self.coupling_1D = False
        self.coupling_0D = False
        self.interpreted_func = False

        # Read param file
        try:
            fdata = open(name, 'r', encoding='utf-8')
        except:
            err = '\n  Unable to open the SYRTHES data file ' + name + '\n    -> Error !\n'
            sys.stderr.write(err)
            sys.exit('Stop SYRTHES execution.')

        coupl_apps = []

        # Delete comment lines
        kwds = []
        [kwds.append(ch.strip()) for ch in fdata if ch[0:1]!='/']

        # Scan param file
        for kw in kwds:

            kws = kw.split("=")
            k = kws[0].strip()

            if k == "MAILLAGE CONDUCTION":
                self.c_mesh_name = kws[1].strip().replace("\\", "/")
            elif k == "MAILLAGE RAYONNEMENT":
                self.r_mesh_name = kws[1].strip().replace("\\", "/")
            elif k == "MAILLAGE FL1D":
                self.f1d_mesh_name = kws[1].strip().replace("\\", "/")
            elif k == "PREFIXE DU RESULTAT PRECEDENT POUR SUITE DE CALCUL":
                self.restart_file = kws[1].strip().replace("\\", "/")
            elif k == "PREFIXE DES FICHIERS RESULTATS":
                self.result_name = kws[1].strip().replace("\\", "/")
            elif k == "CLIM":
                kk = kws[1].strip()
                kk = kk.strip('0123456789_')
                kk = kk.split()
                if kk[0] == "COUPLAGE_SURF_FLUIDE":
                    self.coupling = True
                elif kk[0] == "COUPLAGE_VOL_FLUIDE":
                    self.coupling = True
            elif k == "SUITE DE CALCUL":
                if kws[1].strip() == 'OUI':
                    self.restart = True
            #fluid 1D
            elif k == "MODELE 1D FLUIDE":
                if kws[1].strip() == 'OUI':
                    self.coupling_1D= True
            #fluid 0D
            elif k == "MODELE 0D FLUIDE":
                if kws[1].strip() == 'OUI':
                    self.coupling_0D= True
            elif k.count('_FCT') > 0:
                self.interpreted_func = True
            elif k.count('_INDOOR') > 0:
                self.interpreted_func = True
            elif k.count('_OUTDOOR') > 0:
                self.interpreted_func = True

        # Treatment of the conductive mesh name
        if self.c_mesh_name == None or len(self.c_mesh_name) == 0:
            err = '\n  Conduction mesh name is not specified.\n' \
                  '     -> Error !\n'
            sys.stderr.write(err)
            sys.exit('Stop SYRTHES execution.')

        self.c_mesh_prefix = os.path.basename(self.c_mesh_name)
        self.c_mesh_suffix = self.c_mesh_name[len(self.c_mesh_name)-4:]

        if self.c_mesh_suffix != ".syr":
            err = '\n  Invivalid mesh file format.\n' \
                  '  Only SYRTHES format ".syr" is allowed' \
                  '  If not, convert your format with SYRTHES tool "convert2syrthes"\n'
            sys.stderr.write(err)
            sys.exit('Stop SYRTHES execution.')

        # Treatment of the radiative mesh name
        if self.r_mesh_name != None and len(self.r_mesh_name) > 0:

            self.r_mesh_prefix = os.path.basename(self.r_mesh_name)
            self.r_mesh_suffix = self.r_mesh_name[len(self.r_mesh_name)-4:]

            if self.r_mesh_suffix != ".syr":
                err = '\n  Invalid mesh file format.\n' \
                      '  Only SYRTHES format ".syr" is allowed' \
                      '  If not, convert your format with SYRTHES tool "convert2syrthes"\n'
                sys.stderr.write(err)
                sys.exit('Stop SYRTHES execution.')

        # Treatment of the 1D fluid mesh name
        if self.f1d_mesh_name != None and len(self.f1d_mesh_name) > 0:

            self.f1d_mesh_prefix = os.path.basename(self.f1d_mesh_name)
            self.f1d_mesh_suffix = self.f1d_mesh_name[len(self.f1d_mesh_name)-4:]

            if self.f1d_mesh_suffix != ".syr":
                err = '\n  Invalid mesh file format.\n' \
                      '  Only SYRTHES format ".syr" is allowed' \
                      '  If not, convert your format with SYRTHES tool "convert2syrthes"\n'
                sys.stderr.write(err)
                sys.exit('Stop SYRTHES execution.')

#isabelle
        # Treatment of the result name
        if self.result_name == None or len(self.result_name) == 0:
            err = '\n  Results prefix is not specified.\n' \
                  '     -> Error !\n'
            sys.stderr.write(err)
            sys.exit('Stop SYRTHES execution.')


        fdata.close()

#***************************************************************************************
    def update(self, n_procs, c_mesh_abspath, r_mesh_abspath, f1d_mesh_abspath, restart_abspath,new_filename = "tmp.data"):
#***************************************************************************************

# Read param file
        try:
            fdata = open(self.name, 'r', encoding='utf-8')
        except:
            err = '\n  Unable to open the SYRTHES data file ' + self.name
            err += '\n   -> Error !\n'
            sys.stderr.write(err)
            return 1

        new_fdata = open(new_filename, 'w', encoding='utf-8')

        # Delete comment lines
        kwds = []
        [kwds.append(ch.strip()) for ch in fdata if (ch[0:1]!='/' or (ch[0:1]=='/' and ch[1:2]=='&')) ]

        # Scan param file
        for i in range(len(kwds)):

            kw = kwds[i]
            kws = kw.split("=")
            k = kws[0].strip()

            if k == "MAILLAGE CONDUCTION":
                kws[1] = c_mesh_abspath
                kws[0] += '='
                kwds[i] = " ".join(kws)

            elif k == "MAILLAGE RAYONNEMENT" and r_mesh_abspath != None:
                kws[1] = r_mesh_abspath
                kws[0] += '='
                kwds[i] = " ".join(kws)

            elif k == "MAILLAGE FL1D" and f1d_mesh_abspath != None:
                kws[1] = f1d_mesh_abspath
                kws[0] += '='
                kwds[i] = " ".join(kws)

            elif self.restart == True and k == "PREFIXE DU RESULTAT PRECEDENT POUR SUITE DE CALCUL":
                kws[1] = restart_abspath
                kws[0] += '='
                kwds[i] = " ".join(kws)

            elif k == "PREFIXE DES FICHIERS RESULTATS":
                init_name = kws[1].strip()
# isabelle
                if n_procs > 1 or self.coupling == True:
                    init_name = os.path.basename(init_name)

                if n_procs > 1:
                    kws[1] = os.path.join('PART', init_name)
                else:
                    kws[1] = init_name
                kws[0] += '='
                kwds[i] = " ".join(kws)

        for kw in kwds:
            new_fdata.write(kw)
            new_fdata.write('\n')

        fdata.close()
        new_fdata.close()

        return 0

#***************************************************************************************
    def dump(self):
#***************************************************************************************
        print("\n   SYRTHES Param summary")
        print("    Param file name =           ", self.name)
        print("    Conduction mesh name =      ", self.c_mesh_name)
        print("    Radiation mesh name =       ", self.r_mesh_name)
        print("    1D fluid mesh name =        ", self.f1d_mesh_name)
        print("    Result prefix. =            ", self.result_name)
        if self.restart == True:
            print( "    Restart file =              ", self.restart_file)

        print("    Restart =                   ", str(self.restart))
        print("    Coupling =                  ", str(self.coupling))
        print("    1D Coupling =               ", str(self.coupling_1D))
        print("    0D Coupling =               ", str(self.coupling_0D))
        print("    Interpreted functions =     ", str(self.interpreted_func))

# ***********************
# END OF CLASS DEFINITION
# ***********************

#***************************************************************************************
def create_syrcase(casedir):
#***************************************************************************************

    retval = 0

    # Check syrthes environment and set paths
    if s4home == None:
        retval = check_and_load_env()
        if retval != 0:
            err = '\n   Error during checking and loading environment.\n'
            sys.stderr.write(err)
            return 1

    if casedir == None or len(casedir) == 0:
        err = '\n  No path defined for creating SYRTHES case!\n     -> Error !\n'
        sys.stdout.write(err)
        return 1

    # Call shell script to create a new case
    try:
        s4case = os.path.join(s4home, 'bin', 'syrthes_create_case')
        subprocess.Popen([s4case, casedir]).communicate()
    except:
        err = '   -> Error !\n   Cannot create syrthes case directory.\n'
        sys.stderr.write(err)
        return 1

    return retval

#-------------------------------------------------------------------------------
# Processes the passed command line arguments and defined its relative syrthes
# case object
#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """

    # Parse command line

    parser = OptionParser(usage="usage: %prog 4.0 [options].\nType %prog -h")

    parser.add_option("-d", "--data",
                      action="store", type="string",
                      dest="data_file",
                      help="Name of the SYRTHES data file")
    parser.add_option("-n", "--nbprocs",
                      action="store", type="int",
                      dest="n_procs",
                      help="Number of processors for SYRTHES computation")
    parser.add_option("-t", "--toolpart",
                      action="store", type="string",
                      dest="part_tool_name",
                      help="Tool for partitionning [metis] or [scotch]")
    parser.add_option("-r",
                      action="store", type="int",
                      dest="n_procs_r",
                      help="Number of processors for SYRTHES radiation computation")
    parser.add_option("--name",
                      action="store", type="string",
                      dest="syr_name",
                      help="List of coupled app. codes")
    parser.add_option("-l", "--log",
                      action="store", type="string",
                      dest="logfile",
                      help="Name of the execution log file")
    parser.add_option("-p","--no-prepro",
                      action="store_false",
                      dest="prepro",
                      help="Don't run the SYRTHES pre-processor")
    parser.add_option("-v","--visu",
                      action="store", type="string",
                      dest="post_mode",
                      help="Define postprocessing format [ensight, med]")
    parser.add_option("-g","--debug",
                      action="store_true",
                      dest="debug",
                      help="Run SYRTHES in debug")
    parser.add_option("--exec-dir",
                      action="store", type="string",
                      dest="exec_dir",
                      help="Execution directory (default: current directory)")
    parser.add_option("--data-dir",
                      action="store", type="string",
                      dest="data_dir",
                      help="Data directory (default: current directory)")
    parser.add_option("--src-dir",
                      action="store", type="string",
                      dest="src_dir",
                      help="Source directory (default: current directory)")

    parser.set_defaults(n_procs=1)
    parser.set_defaults(n_procs_r=-1)
    parser.set_defaults(part_tool_name="metis")
    parser.set_defaults(syr_name="SYR")
    parser.set_defaults(prepro=True)
    parser.set_defaults(post_mode=None)
    parser.set_defaults(debug=False)
    parser.set_defaults(exec_dir=None)
    parser.set_defaults(data_dir=None)
    parser.set_defaults(src_dir=None)

    (options, args) = parser.parse_args(argv)

    # Error management
    if not options.data_file:
        parser.print_help()
        err = "\n Incorrect command line definition.\n Define a data file.\n"
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    if not os.path.isfile(options.data_file):
        parser.print_help()
        err = '\n  You must define the name of the data file\n' \
              '  Usage : syrthes.py  -n nb_proc  -d syrthes_data.syd\n\n'
        err += '  ' + options.data_file + ' is not a valid file.\n\n'
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    if int(options.n_procs) < 1:
        parser.print_help()
        err = '\n  You have to define the number of processors : -n [1..n]\n' \
              '  Usage : syrthes.py  -n nb_proc  -d syrthes.data\n\n'
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    if options.post_mode != None:
        if options.post_mode == "med" or options.post_mode == "MED":
            post_mode = 'med'
        elif options.post_mode == "ens" or options.post_mode == "ensight" or options.post_mode == "ENSIGHT":
            post_mode = 'ens'
        else:
            parser.print_help()
            err = '\n  Incorrect argument for post-processing mode\n'
            err += '  Usage : syrthes.py  -v ens\n'
            err += '  Usage : syrthes.py  -v ensight\n'
            err += '  Usage : syrthes.py  -v med\n'
            sys.stderr.write(err)
            sys.exit('Stop SYRTHES execution.')
    else:
        post_mode = None

    # Create a new syrthes case object and initialize by default

    syr_case = SyrthesCase(name = options.syr_name,
                           data_file = options.data_file,
                           data_dir = options.data_dir,
                           exec_dir = options.exec_dir,
                           src_dir = options.src_dir,
                           n_procs = options.n_procs,
                           n_procs_ray = options.n_procs_r,
                           part_tool_name = options.part_tool_name,
                           post_mode = post_mode,
                           debug = options.debug,
                           prepro = options.prepro)

    # Update syrthes case parameters

    if options.logfile:
        syr_case.set_logfile(options.logfile)

    if syr_case.n_procs_ray > syr_case.n_procs:
        err = '\n  Number of processors for radiation must be less than the' \
              ' total number of processors : -r [1..n]\n' \
              '  Usage : syrthes.py  -n nb_proc  -d syrthes.data -r nb_proc_rad\n'
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    return syr_case


#***************************************************************************************
def read_syrthes_param(data_file = "syrthes.data"):
#***************************************************************************************

    return SyrthesParam(name = data_file)


#***************************************************************************************
def build_syrthes(parall = True, cfd_coupling = False, fname = None,
                  debug = False, srcdir = None, destdir = None):
#***************************************************************************************

    # Display comment on the current step
    sys.stdout.write(' Building the executable file syrthes.. \n')

    call_dir = os.getcwd()
    if srcdir != None:
        if not os.path.isdir(srcdir):
            err = "   " + srcdir + " is not a valid directory !"
            sys.stderr.write(err)
            return 1
        os.chdir(srcdir)

    # Check if Makefile is in the current working directory
    makefile=os.path.join(os.getcwd(), "Makefile")
    if not os.path.isfile(makefile):
        err = '   Makefile not found in the current working directory: ' + os.getcwd()
        sys.stderr.write(err)
        return 1

    # Output file
    if fname != None:
        logfname = os.path.abspath(fname)
        log = open(logfname, 'w', encoding='utf-8')

    # Append command line
    cmd='make'
    if parall == False and cfd_coupling == False:
        cmd += ' MPI=no'
    if cfd_coupling == True:
        cmd += ' CFD=yes'
    if debug == True:
        cmd += ' DEBUG=yes'
    if fname != None:
        cmd += ' 2>&1 >| ' + logfname

    retval = 0
    try:
        retval = os.system(cmd)
    except:
        sys.stderr.write('\n  Error during the compilation stage\n')
        retval = 1
    else:
        if retval != 0:
            sys.stderr.write('\n  Error during the compilation stage\n')

    if fname != None:
        log.close()

    if retval != 1 and destdir != None:
        if not os.path.isdir(destdir):
            sys.stderr.write('\n    ' + destdir + ' is not a valid directory !')
            retval = 1
        else:
            if destdir != os.getcwd():
                exec_name = os.path.join(destdir, "syrthes")
                shutil.move(os.path.abspath("syrthes"), exec_name)

    os.chdir(call_dir)

    if retval == 0:
        sys.stdout.write('\n  *****  SYRTHES compilation and link completed *****\n')

    return retval

#***************************************************************************************
def check_and_create_dir(dirname):
#***************************************************************************************

    retval = 0

    if dirname == None or len(dirname) == 0:
        err = '\n  No name defined for the directory to create!\n     -> Error !\n'
        sys.stdout.write(err)
        return 1

    else:
        if os.path.exists(dirname) == False:
            try:
#                os.makedirs(dirname, mode=755)
                os.makedirs(dirname)
                return 0
            except:
                err = '\n   -> Error !\n   Unable to create ' + dirname + ' directory.\n'
                sys.stderr.write(err)
                return 1
        elif not os.path.isdir(dirname):
            err = '   -> Error !\n   A regular file named ' + dirname + \
                  ' is existing.\n' + '  Unable to create the directory.\n'
            sys.stderr.write(err)
            return 1
        elif os.access(dirname, os.W_OK) == False:
            err = '   -> Error !\n   Cannot write into existing directory: ' + dirname + '\n'
            sys.stderr.write(err)
            return 1

    return retval

#***************************************************************************************
def check_and_load_env():
#***************************************************************************************

    global s4home, arch, mpidir
    global syrthes_pp, syrthes_ppfonc, syrthes_post, syrthes_ensight, syrthes_grouphist, syrthes_med, syrthes_ppfonc_tmc

    # Check SYRTHES is in the $PATH
    if s4home == None:
        try:
            out = subprocess.check_output("which syrthes_create_case",shell=True)
            out=str(out.decode("utf-8"))
            s4bin = os.path.dirname(out)
            s4home = os.path.abspath(os.path.join(s4bin,'..'))
        except:
            err = '   -> Error !\n   check_and_load_env : Cannot find syrthes home directory.\n'
            sys.stderr.write(err)
            return 1

    # Load environnement
    sys.stdout.write(' SYRTHES home directory: ' + s4home + '\n')

    # Find which MPI to use
    try:
        fd = open(os.path.join(s4home, 'bin', 'syrthes.profile'), 'r')
        content = fd.readlines()
        matches = fnmatch.filter(content, 'SYRTHES_MPIPATH*')
        last_match = matches[len(matches)-1]
        mpidir = re.split('=', re.split('\n', last_match)[0])[1]
        fd.close()
    except:
        err = '   -> Error !\n   Cannot find a MPI librairy.\n'
        sys.stderr.write(err)
        return 1

    sys.stdout.write(' MPI home directory: ' + mpidir + '\n')

    # Define syrthes directories
    syrthes_pp = os.path.join(s4home, 'bin','syrthes-pp')
    syrthes_ppfonc = os.path.join(s4home, 'bin', 'syrthes-ppfunc')
    syrthes_ppfonc_tmc = os.path.join(s4home, 'bin', 'syrthes-ppfunc-tmc')
    syrthes_post = os.path.join(s4home, 'bin', 'syrthes-post')
    syrthes_ensight = os.path.join(s4home, 'bin', 'syrthes2ensight')
    syrthes_grouphist = os.path.join(s4home, 'bin', 'syrthes-group-histo')
    syrthes_med = os.path.join(s4home, 'bin', 'syrthes2med')

    return 0

#***************************************************************************************
#***************************************************************************************
def main(argv):
#***************************************************************************************
#***************************************************************************************

    # Check syrthes environment and set paths
    retval = check_and_load_env()
    if retval != 0:
        err = '\n   Error during checking and loading environment.\n'
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    # Build a new syrthes case object relative to the command line
    syr_case = process_cmd_line(argv)

    # Read data file and store parameters
    syr_case.read_data_file()

    # Initialize output file if needed
    syr_case.logfile_init()

    # Prepare SYRTHES execution
    # -------------------------

    print("  -----------------------------------")
    print("  Prepare SYRTHES execution directory")
    print("  -----------------------------------\n")

    retval = syr_case.prepare_run()
    if retval != 0:
        err = '\n   Error during the prepration step (errno: '
        err += str(retval) + '\n'
        sys.stderr.write(err)
        sys.exit('Stop SYRTHES execution.')

    # Sumary of the parameters
    syr_case.dump()

    # SYRTHES preprocessing
    # ---------------------

    # Pre-processing (only if SYRTHES computation is done in parallel)
    retval = syr_case.preprocessing()
    if retval != 0:
        sys.stderr.write('\n  Error during the preprocessing step\n')
        sys.exit('Stop SYRTHES execution.')

    # SYRTHES execution
    # -----------------

    print( "\n  -------------------------")
    print( "  Start SYRTHES computation")
    print( "  -------------------------\n")

    # Run SYRTHES
    retval = syr_case.run()
    if retval != 0:
        sys.stderr.write('\n  Error while running syrthes\n')
        sys.exit('Stop SYRTHES execution.')
    else:
        sys.stdout.write('\n  ***** Successful execution ***** \n\n')

    # Post-processing
    if syr_case.post_mode != None:
        retval = syr_case.postprocessing(mode = syr_case.post_mode)
        if retval != 0:
            sys.stderr.write('\n  Error while post-processing result files.\n')
            sys.exit('Stop SYRTHES execution.')
        else:
            sys.stdout.write('   -> Ok\n')

    sys.stdout.write("\n  *** End ***\n")


if __name__ == "__main__":
    main(sys.argv[1:])
