"""
    Fragon utils.py 

    Copyright (C) 2017-2018 University of York
    Author: Huw Jenkins

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this program.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""
import argparse
import json
import os
import shutil
import sys
import errno
# import logging
from iotbx import pdb
from lxml import etree
from fragon.data import mtz_output
from fragon.version import __version__

def print_header(log):
  log.info('\n     -------------------------------------------------------')
  log.info('     |                                                     |')
  log.info('     |                Fragon version %10s            |'  % __version__)
  log.info('     |                                                     |')
  log.info('     |  Please send bug reports to huw.jenkins@york.ac.uk  |')
  log.info('     |                                                     |')
  log.info('     -------------------------------------------------------\n')

# read command line arguments
def parse_command_line(args):
  parser = argparse.ArgumentParser( prog='fragon', description='Places fragments with Phaser followed by density modification with ACORN')
  optional = parser.add_argument_group('additional options')
  optional.add_argument('--mtz', required=False, default=None, type=str, metavar='data.mtz',
                        help='input MTZ file')
  optional.add_argument('--log', required=False, default=None, type=str, metavar='Fragon_.log',
                        help='logfile')
  optional.add_argument('--i2', required=False,  action='store_true', default=False,
                        help='ccp4i2 mode')
  optional.add_argument('--seq', required=False, default=None, type=str, metavar='protein.seq', help='sequence')
  optional.add_argument('--helix', required=False, default=None, type=int, metavar='n',
                        help='helix length')
  optional.add_argument('--fragment', required=False, default=None, type=str, metavar='fragment.pdb',
                        help='input fragment')
  optional.add_argument('--ensemble', required=False, default=None, type=str, metavar='ensemble.pdb',
                        help='input ensemble')
  optional.add_argument('--rescore_strands', required=False, action='store_true', default=False,\
                        help='rigid body refine each strand in ensemble separately')
  optional.add_argument('--rescore_residues', required=False, action='store_true', default=False,\
                        help='rigid body refine half strands in ensemble separately')
  optional.add_argument('--rescore_all', required=False, action='store_true', default=False,\
                        help='rescore all solutions (not up to number of requested solutions)')
  optional.add_argument('--rescore_models', required=False, action='store_true', default=False,\
                        help='rescore models in ensemble separately')
  optional.add_argument('--fixed', required=False, default=None, type=str, metavar='partial_solution.pdb',
                        help='partial solution postion will be fixed')
  optional.add_argument('--solutions', required=False, type=int, metavar='10', default=10,
                        help='The maximum no. of solutions to try and test, default 10')
  optional.add_argument('--test_all', required=False, action='store_true', default=False,
                        help='test all solutions after ACORN CC indicates success')
  optional.add_argument('--ACORN_definitive_CC', required=False, type=float, metavar='0.300',default=0.30000,
                        help='Stop testing solutions when ACORN CC is above this value')
  optional.add_argument('--ACORN_CC_diff', required=False, type=float, metavar='0.15',default=0.15000,
                        help='Stop testing solutions when difference between best and worse ACORN CC is above this value')
  optional.add_argument('--ncs_copies', required=False, type=int, metavar='1', default=1,
                        help='no. of molecules in asymmetric unit, default = 1')
  optional.add_argument('--copies', required=False, type=int, metavar='1', default=1,
                        help='no. of copies of helix to place, default = 1')
  optional.add_argument('--test_all_plausible_sg', required=False, action='store_true', default=False,
                        help='Test all space groups in same laue group as input data')
  optional.add_argument('--test_sg_list', required=False, type=str, metavar='P 21 21 21, P 21 21 2', default=None,
                        help='list of space groups_to_test')
  optional.add_argument('--name', required=False, default=None, type=str, metavar='My_Protein',
                        help='Name to add to all output files, default = same as mtz file name')
  optional.add_argument('--I', required=False, default=None, type=str, help='I in mtzin')
  optional.add_argument('--SIGI', required=False, default=None, type=str, help='SIGI in mtzin')
  optional.add_argument('--FP', required=False, default=None, type=str, help='FP in mtzin')
  optional.add_argument('--SIGFP', required=False, default=None, type=str, help='SIGFP in mtzin')
  optional.add_argument('--solvent', required=False, type=float, metavar='0.nn', help='solvent content (fraction) to override that calculated from sequence and no. of NCS copies')
  optional.add_argument('--search_highres', metavar='1.5', required=False, default=None,
                        type=float, help='High resolution for MR search, default = high resolution of data')
  optional.add_argument('--search_lowres', metavar='50', required=False, default=None,
                        type=float, help='Low resolution for MR search, default = low resolution of data')
  optional.add_argument('--rotation_peaks', metavar='1000', required=False, default=None,
                        type=float, help='Select n peaks from FRF')
  optional.add_argument('--rotation_cluster_off', required=False, action='store_true', default=False,
                        help='Turn off peak clustering in FRF')
  optional.add_argument('--rotation_sampling', required=False, metavar='1.5', default=None,
                        type=float, help='Sampling for FRF, default = decided by Phaser')
  optional.add_argument('--translation_sampling', required=False, metavar='0.5', default=None,
                        type=float, help='Sampling for FTF, default = decided by Phaser')
  optional.add_argument('--search_down_percent', required=False, type=int, metavar='15',default=None,
                        help='Percentage of extra solutions for deep search')
  optional.add_argument('--purge',  type=int, metavar='100',default=None,
                        help='for multicopy searches purge all but top N solutions after refinement of each component')
  optional.add_argument('--definitive_TFZ', required=False, metavar='8.0', default=8.0,
                        type=float, help='TFZ indicating solution, default = 8')
  optional.add_argument('--RMS', required=False, metavar='0.2', default=0.2,
                        type=float, help='RMS error of helix, default = 0.2')
  optional.add_argument('--results_json', required=False, type=str, metavar='/path/to/Fragon_N/name_fragment_results.json', default=None,
                        help='restart from this file')
  optional.add_argument('--input_only', required=False, action='store_true', default=False,
                        help='Do not test both hands of enantiomorphic space groups (to speed up debugging)')
  optional.add_argument('--no_tncs', required=False, action='store_true', default=False,
                        help='Turn off tNCS correction')
  optional.add_argument('--debug', required=False, action='store_true', default=False,
                        help='turn on debugging messages')
  optional.add_argument('--nproc', required=False, type=int, metavar='1', default=1,
                        help='no. threads')
  if len(args) == 0:
    parser.print_usage()
    sys.exit()
  return parser.parse_args(args)

def setup_restart(results_json, log):
  results = read_results_json(results_json)
  run_dir, results_file = os.path.split(results_json)
  log.debug('DEBUG results read from JSON %s' % results)
  restart = {'start_point':results['stage'], 'name':results['name'], 'root':results['root'],'mtzin':results['mtzin']}
  log.info('Restarting from %s' % results_file)
  log.info('')
  log.info('Run_dir: %s  start point: %s' % (run_dir, restart['start_point']))
  log.info('')
  log.debug('DEBUG reading json: %s' % results['solutions'])
  logfile = results['root'] +'_restart.log'
  run_dir = os.path.abspath(run_dir)
  restart.update({'run_dir':run_dir})
  if not os.path.isfile(os.path.join(run_dir, results['mtzin'])):
    sys.exit('File not found: %s' % results['mtzin'])
  if not os.path.isdir(run_dir):
    raise RuntimeError('run dir %s does not exist' % run_dir)
  tempdir = os.path.join(run_dir,'TEMP')
  restart.update({'tempdir':tempdir})
  if os.path.isdir(tempdir):
    for tempfile in os.listdir(tempdir):
      file_path = os.path.join(tempdir, tempfile)
      if os.path.isfile(file_path):
        os.unlink(file_path)
  return restart, results['solutions'], results['search'], results['scoring'], logfile

def setup_run(mtzin, seqin, pdbin, pdbin_fixed):
  run_dir_base='Fragon_'
  run_number = 1
  while run_number <= 100:
    run_dir_name = run_dir_base + str(run_number)
    try:
      os.mkdir(run_dir_name)
      run_dir = os.path.abspath(run_dir_name)
      break
    except OSError as exception:
      if exception.errno == errno.EEXIST:
        run_number += 1
      elif exception.errno == errno.EACCES:
        raise RuntimeError('No write permission to this directory')
  else:
    raise RuntimeError('Maximum run number exceded')
  tempdir_name = os.path.join(run_dir,'TEMP')
  try:
    os.mkdir(tempdir_name)
    tempdir = os.path.abspath(tempdir_name)
  except OSError as exception:
    if exception.errno == errno.EACCES:
      raise RuntimeError('No write permission to this directory')

  # copy files
  shutil.copy(mtzin, run_dir)
  if seqin is not None:
    shutil.copy(seqin, run_dir)
  if pdbin is not None:
    if pdbin[:-4].split('-')[0] == 'Helix' and pdbin[:-4].split('-')[1] in [str(n) for n in range(1,71)]:
      shutil.move(pdbin, run_dir)
    else:
      shutil.copy(pdbin, run_dir)
  if pdbin_fixed is not None:
    shutil.copy(pdbin_fixed, run_dir)
  return run_dir, tempdir

def define_search(args):
  search = {}
  if args.mtz is None and args.results_json is None:
    sys.exit('--mtz is required if not restarting')
  elif args.results_json is None:
    mtzin = args.mtz
  if not os.path.isfile(mtzin):
    sys.exit('File not found: %s' % mtzin)
  if args.name is not None:
    name = args.name
    if len(name.split()) > 1:
      sys.exit('Sorry Name cannot have spaces')
  else:
    name = os.path.split(os.path.abspath(args.mtz))[1][:-4]
  search['have_ensemble'] = False
  search['num_models'] = 1
  search['copies'] = args.copies
  if search['copies'] > 1:
    search['fragment'] = '1 copy of '
  else:
    search['fragment'] = '%d copies of ' % search['copies']
  if args.helix is not None:
    helix_length = args.helix
    make_helix(helix_length)
    search['pdbin'] = 'Helix-' + str(helix_length) + '.pdb'
    root = name + '_helix' + str(helix_length)
    logfile = name + '_Fragon_' + str(search['copies']) + 'x_helix' +  str(helix_length) + '.log'
    search['fragment'] += '%d residue ideal helix' % helix_length
  elif args.fragment is not None:
    search['pdbin'] = args.fragment
    if not os.path.isfile(search['pdbin']):
      sys.exit('File not found: %s' % search['pdbin'])
    fragment_name = os.path.split(os.path.abspath(search['pdbin']))[1][:-4]
    root = name + '_' + fragment_name
    logfile = name + '_Fragon_' + str(search['copies']) + 'x_' + fragment_name + '.log'
    search['fragment'] += 'fragment %s' % fragment_name
  elif args.ensemble is not None:
    search['pdbin'] = args.ensemble
    if args.i2:
      search['pdbin'] = os.path.join(os.environ['CCP4'], 'share', 'fragon', 'include', 'fragments', search['pdbin'])
    if not os.path.isfile(search['pdbin']):
      sys.exit('File not found: %s' % search['pdbin'])
    fragment_name = os.path.split(os.path.abspath(search['pdbin']))[1][:-4]
    root = name + '_' + fragment_name
    logfile = name + '_Fragon_' + str(search['copies']) + 'x_' + fragment_name + '.log'
    search['have_ensemble'] = True
    search['num_models'] = num_models(search['pdbin'])
    search['fragment'] += 'ensemble %s of %d models' % (fragment_name, search['num_models'])
  else:
    sys.exit('You need to set either --helix_length or supply a fragment file')
  search['pdbin_fixed'] = args.fixed
  search['rescore_strands'] = args.rescore_strands
  search['rescore_residues'] = args.rescore_residues
  search['rescore_all'] = args.rescore_all
  search['rescore_models'] = args.rescore_models
  if args.rescore_models:
    search['rescore_all'] = True
  if args.seq is not None:
    search['seqin'] = args.seq
    search['solvent'] = None
    if not os.path.isfile(search['seqin']):
      sys.exit('File not found: %s' % search['seqin'])
  else:
    if args.solvent is None:
      sys.exit("You must set solvent content with --solvent if you don't supply a sequence")
    else:
      search['solvent'] = args.solvent
      search['seqin'] = None
  search['ncs_copies'] = args.ncs_copies
  search['num_solutions'] = args.solutions
  # non default Phaser options
  search['rot_peaks'] = args.rotation_peaks
  search['purge'] = args.purge
  # turn on by default now (as in CCP4i2 wrapper)
  if search['purge'] == None and search['copies'] > 1:
    search['purge'] = 100
  search['rot_cluster_off'] = args.rotation_cluster_off
  search['rot_samp'] = args.rotation_sampling
  search['tra_samp'] = args.translation_sampling
  search['search_down'] = args.search_down_percent
  search['sgalternative'] = args.test_all_plausible_sg
  search['input_hand'] = args.input_only
  search['tncs'] = args.no_tncs
  search['sg_list'] = args.test_sg_list
  search['tfz_solved'] = args.definitive_TFZ
  search['highres'] = args.search_highres
  search['lowres'] = args.search_lowres
  search['rms'] = args.RMS
  search['fragment'] += ' with %0.1f Angstrom RMS error' % search['rms']
  return name, root, mtzin, search, logfile

def get_labels(args):
  if args.I is not None:
    i = args.I
    if args.SIGI is not None:
      sigi = args.SIGI
      fp = None
      sigfp = None
    else:
      sys.exit('You must set --SIGI if you set --I')
  elif args.FP is not None:
    fp = args.FP
    if args.SIGFP is not None:
      sigfp = args.SIGFP
      i = None
      sigi = None
    else:
      sys.exit('You must set --SIGFP if you set --FP')
  # check user hasn't set I and FP
  if i is not None or sigi is not None:
    if fp is not None or sigfp is not None:
      sys.exit('Only set --I, --SIGI (preferably) or --FP, --SIGFP not both')
  else:
    return i, sigi, fp, sigfp

def define_scoring(args, log):
  test_all = args.test_all
  if test_all:
    #hard coded
    acornCC_solved = 0.20000
  else:
    acornCC_solved = args.ACORN_definitive_CC
  acornCC_diff = args.ACORN_CC_diff
  return {'acornCC_solved':acornCC_solved, 'acornCC_diff':acornCC_diff, 'test_all':test_all}

# function to format times nicely
def print_time(process, seconds, log):
  days = int(seconds/60/60/24)
  hours = int(seconds/60/60-24*days)
  mins= int(seconds/60-24*60*days-60*hours)
  secs = seconds-24*60*60*days-60*60*hours-60*mins
  time_string = '%d days %d hours %d mins %0.2f seconds (%0.2f s)' % (days, hours, mins, secs, seconds)
  log.info('%s: %s' % (process, time_string))
  return time_string

# make helix from 70-mer ideal helix
def make_helix(helix_length):
  if helix_length < 1 or helix_length > 70:
    sys.exit('Helix length must be between 1 and 70 residues')
  theor_helix = os.path.join(os.environ['CCP4'], 'share', 'fragon', 'include', 'fragments', 'ideal_helix.pdb')
  output = 'Helix-' + str(helix_length) + '.pdb'
  num_lines = 0
  if os.path.isfile(output):
    os.remove(output)
  with open(theor_helix, 'r') as pdbin:
    for line in pdbin:
      if num_lines < 5*helix_length:
        with open(output, 'a') as pdbout:
          pdbout.write(line)
        num_lines += 1
  with open(output, 'a') as pdbout:
    pdbout.write('END')

def num_models(ensemble):
  pdbin =  pdb.hierarchy.input(ensemble)
  return len(pdbin.hierarchy.models())

def get_solutions(phaser_solutions, root, num_solutions):
  solutions = []
  solution_number = 0
  if phaser_solutions != 'No solutions':
    for phaser_solution in phaser_solutions:
      solution_number = phaser_solution.NUM
      if solution_number <= num_solutions:
        solution = {}
        solution_id = root + '.' + str(solution_number)
        solution_sg = phaser_solution.getSpaceGroupName()
        solution_llg = phaser_solution.LLG
        solution_tfz = phaser_solution.TFZ
        solution.update({'number': solution_number, 'id':solution_id, 'sg':solution_sg, 'llg':solution_llg, 'tfz':solution_tfz})
        solutions.append(solution)
      else:
        break
  return solutions

def write_results_json(version=None, results_json=None, name=None, root=None, mtzin=None, stage=None, solutions=None, search=None, scoring=None):
  if version is not None:
    output = dict(Fragon=version, name=name, root=root, mtzin=mtzin, 
                  stage=stage, solutions=solutions, search=search, scoring=scoring)
  else:
    output = solutions
  with open(results_json, 'w') as results:
    print >> results, json.dumps(output, sort_keys=True, indent=2, separators=(',', ': '))

def write_output(items, json_file=None, xml_file=None, xmlroot=None, output=None):
  # in non-i2 mode items are added to the output dictionary which is dumped to json
  if json_file is not None:
    if 'result' in items.keys():
      result = items['result']
      for solution in output['solutions']:
        if solution['id'] == result['id']:
          solution.update({'acornCC':result['acornCC']})
    else:
      output.update(items)
    temp_filename = json_file + '.tmp'
    with open(temp_filename, 'w') as jsonfile:
      print >> jsonfile, json.dumps(output, sort_keys=True, indent=2, separators=(',', ': '))
    os.rename(temp_filename, json_file)
    return output
  elif xmlroot is None:
    xmlroot = etree.Element('Fragon')
    return xmlroot
  elif xml_file is not None:
    # in i2 mode new items are added to the etree as this preserves the order in the xml
    for key in items.keys():
      if key == 'Fragon':
        version_node = etree.SubElement(xmlroot, 'Version')
        version_node.text = output['Fragon']
      elif key == 'callback':
        callback = items['callback']
        if callback[0] == 'progress':
          try:
            progress_node = xmlroot.xpath('//Fragon/phaser_progress')[0]
          except IndexError:
            progress_node = etree.SubElement(xmlroot, 'phaser_progress')
          progress_node.text = callback[1]
        elif callback[0] == 'Best LLG':
          best_llg_node = etree.SubElement(xmlroot, 'best_llg')
          best_llg_node.text = callback[1]
        elif callback[0] == 'Best TFZ':
          best_tfz_node = etree.SubElement(xmlroot, 'best_tfz')
          best_tfz_node.text = callback[1]
      elif key == 'solutions':
        solutions = items['solutions']
        try:
          solutions_node = xmlroot.xpath('//Fragon/solutions')[0]
        except IndexError:
          solutions_node = etree.SubElement(xmlroot, 'solutions')
        if len(solutions) > 0:
          solutions_node.text = json.dumps(solutions)
      else:
        node = etree.SubElement(xmlroot, key)
        node.text = items[key].__str__()
    temp_filename = 'program.xml.tmp'
    with open(temp_filename, 'w') as xmlfile: xmlfile.write(etree.tostring(xmlroot, pretty_print=True))
    os.rename(temp_filename, xml_file)

def read_results_json(results_json):
  with open(results_json, 'r') as results:
    loaded_results = json.load(results)
  # need to convert unicode solution_id and name to string otherwise:
  #   Boost.Python.ArgumentError: Python argument types in
  #       object.write(object, unicode)
  #   did not match C++ signature:
  #       write(iotbx::mtz::object {lvalue}, char const* file_name)
  if loaded_results['name'] is not None:
    loaded_results['name'] = str(loaded_results['name'])
    loaded_results['mtzin'] = str(loaded_results['mtzin'])
    for solution in loaded_results['solutions']:
      solution['id'] = str(solution['id'])
  return loaded_results

def count_tested(solutions, key):
  solutions_tested = 0
  number_solutions = 0
  for solution in solutions:
    number_solutions += 1
    if solution[key] is not None:
      solutions_tested += 1
  assert number_solutions == len(solutions)
  return solutions_tested

def write_output_files(name, best_solution_id, mtzin, log):
  log.info('ACORN CC greater than 0.2 suggests map is useful for autobuilding\n')
  log.info('Copying %s to %s_phaser_solution.pdb' % (best_solution_id+'.pdb', name))
  shutil.copy(best_solution_id+'.pdb', name+'_phaser_solution.pdb')
  log.info('Copying %s to %s_phaser_solution.mtz' % (best_solution_id+'.mtz', name))
  shutil.copy(best_solution_id+'.mtz', name+'_phaser_solution.mtz')
  acorn_mtz = best_solution_id+'.acorn.mtz'
  mtzout = name + '_acorn_phases.mtz'
  log.info('Copying %s to %s' % (acorn_mtz, mtzout))
  mtz_output(mtzin, acorn_mtz, mtzout)

def print_refs(log):
  log.info('\nIf you solve a structure with Fragon please cite:')
  log.info('\nPhaser crystallographic software')
  log.info('McCoy AJ, Grosse-Kunstleve RW, Adams PD, Winn MD, Storoni LC & Read RJ.')
  log.info('J. Appl. Cryst. (2007). 40, 658-674 (http://doi.org/10.1107/S0021889807021206)')
  log.info('\nA modified ACORN to solve protein structures at resolutions of 1.7 A or better')
  log.info('Yao J-X, Woolfson MM, Wilson KS & Dodson EJ.')
  log.info('Acta Cryst. (2005). D61, 1465-1475 (http://dx.doi.org/10.1107/S090744490502576X)\n')
