"""
This python script generates MESA run directories based on a template directory
for a specified range of initial masses.
"""
from distutils.dir_util import copy_tree
import shutil
from pathlib import Path
import numpy as np
import paths # showyourwork provides this

def makedir(path):
    if not path.exists():
        print('creating directory {}'.format(path))
        path.mkdir()
    else:
        print('{} already exists'.format(path))



masses_01  = np.arange(1.2,5.,0.1)
minDmix_01 = [10]*masses_01.size
masses_02  = np.arange(5,10,0.25)
minDmix_02 = [100]*masses_02.size
masses_03  = np.arange(10,41,1.)
minDmix_03 = [1000]*masses_03.size
masses  = np.concatenate([masses_01, masses_02, masses_03])
minDmix = np.concatenate([minDmix_01, minDmix_02, minDmix_03])

dbconvpen = [True, False]
Z_values = [0.014]

template_dir = paths.data / 'MESA_input/'
print('template_dir:', template_dir)
grid_dir = Path('./carbon_dep_models/')
makedir(grid_dir)

replacements = dict()
replacements['LOGS_DIR']  = "'LOGS'"
for do_pen in dbconvpen:
    for Zv in Z_values:
        Z_str = '{:.3f}'.format(Zv)
        replacements['METALLICITY'] = Z_str
        Z_dir = grid_dir.joinpath('Z{}/'.format(Z_str))
        makedir(Z_dir)
        for M, Dm in zip(masses, minDmix):

            if do_pen:
                this_template_dir = template_dir / 'template_dbcp/'
                run_tag = 'dbcp_'
            else:
                this_template_dir = template_dir / 'template_standard/'
                run_tag = 'standard_'
            run_tag += 'carbon_dep_M{:>05.1f}'.format(M)
            run_dir = Z_dir.joinpath('work_{}/'.format(run_tag))
            copy_tree(str(this_template_dir), str(run_dir))

            replacements['MASS'] = '{:.1f}'.format(M)
            replacements['DMIX_VALUE'] = '{}'.format(Dm)

            for filename in ['clean', 'mk', 're', 'rn', 'stash.py', \
                    'history_columns.list', 'profile_columns.list', 'inlist_xtra_coeff_mesh']:
                template_path = template_dir.joinpath(filename)
                run_path = run_dir.joinpath(filename)
                shutil.copy(str(template_path), str(run_path))

            for directory in ['make', 'src']:
                template_path = template_dir.joinpath(directory)
                run_path = run_dir.joinpath(directory)
                if not run_path.exists():
                    shutil.copytree(str(template_path), str(run_path))

            #copy inlist_project_temp w/ appropriate changes
            for inlist_name in ['inlist_specifics_zams', 'inlist_specifics_tams', 'inlist_specifics_hedep', 'inlist_specifics_cdep']:
                tmp_path    = run_dir.joinpath('tmp')
                inlist_path = run_dir.joinpath(inlist_name)
                #move inlist template to tmp
                shutil.move(str(inlist_path), str(tmp_path))
                template_inlist = open(str(tmp_path), 'r')
                run_inlist = open(str(inlist_path), 'w')
                for line in template_inlist.readlines():
                    for key, string in replacements.items():
                        if key in line:
                            line = line.replace(key, string)
                    run_inlist.write(line)
                template_inlist.close()
                run_inlist.close()

            #copy evan's pleiades job submission file (may need to change for other people)
            ### TODO: if you're someone else, you just need to make sure you have an appropriate script to run the grid.
            ### My pleiades job script just loads MESA 22.11.1 and then does a "./clean", a "./mk" and a "./rn inlist_standard"
            if template_dir.joinpath('evan_pleiades_job_submit').exists():
                template_job = open(str(template_dir.joinpath('evan_pleiades_job_submit')), 'r')
                run_job = open(str(run_dir.joinpath('evan_pleiades_job_submit')), 'w')
                for line in template_job.readlines():
                    if "#PBS -N mesa_standard" in line:
                        line = "#PBS -N mesa_{}_Z{}\n".format(run_tag, Z_str)
                    if 'inlist_standard' in line:
                        line = line.replace('inlist_standard', inlist_name)
                    run_job.write(line)
                template_job.close()
                run_job.close()
            else:
                print('cannot find "evan_pleiades_job_submit"; not copying it')
