"""
    Fragon score.py

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
from __future__ import absolute_import, division, print_function
import glob
import os
import logging
import shlex
import subprocess
from fragon.data import massage_data, tidy_data, extend_ecalc
from fragon.utils import write_output, write_results_json, read_results_json
from libtbx import easy_mp

log = logging.getLogger(__name__)

class setup_solution_tester(object):
  def __init__(self, scores, tempdir, test_all, acornCC_solved, acornCC_diff, i, sigi, fp, sigfp, lowres, highres, solvent, formfactors):
    self.scores = scores
    self.tempdir = tempdir
    self.test_all = test_all
    self.acornCC_solved = acornCC_solved
    self.acornCC_diff = acornCC_diff
    self.i = i
    self.sigi = sigi
    self.fp = fp
    self.sigfp = sigfp
    self.lowres = lowres
    self.highres = highres
    self.solvent = solvent
    self.formfactors = formfactors

  def __call__(self, solution):
    if len(self.scores) == 0: # first callback
      for acorn_json in glob.glob('*.acorn.json'):
        try:
          self.scores.append(read_results_json(acorn_json)['acornCC'])
        except IOError:
          self.log = logging.getLogger(__name__)
          self.log.debug('DEBUG Error reading file: %s' % acorn_json)
    elif not_definitive(self.scores, self.acornCC_solved, self.acornCC_diff): # avoid needless file reads
      try:
        self.scores.extend([ read_results_json(acorn_json)['acornCC'] for acorn_json in glob.glob('*.acorn.json') if read_results_json(acorn_json)['acornCC'] not in self.scores])
      except IOError:
        self.log = logging.getLogger(__name__)
        self.log.debug('DEBUG Error reading file: %s' % acorn_json)
    return test_solution(tempdir=self.tempdir,
                         scores=self.scores,
                         solution=solution,
                         test_all=self.test_all,
                         acornCC_solved=self.acornCC_solved,
                         acornCC_diff=self.acornCC_diff,
                         i=self.i,
                         sigi=self.sigi,
                         fp=self.fp,
                         sigfp=self.sigfp,
                         lowres=self.lowres,
                         highres=self.highres,
                         solvent=self.solvent,
                         formfactors=self.formfactors)

def run_acorn(mtzin, xyzin, mtzout, lowres, highres, solvent, formfactors, acorn_script, acorn_logfile):
  with open(acorn_script,'w') as com:
    print('labin E=E_ISO FP=F_ISO SIGFP=SIGF_ISO', file=com)
    print('RESO %0.2f %0.2f' % (lowres,highres), file=com)
    if highres > 1.0:
      print('EXTEND 50.0 1.0', file=com)
    else:
      print('EXTEND 50.0 %0.2f' % highres, file=com)
    print('ESTRO 0.8', file=com)
    print('POSI 1', file=com)
    print('SOLV %s' % solvent, file=com)
    print('NTRY 1', file=com)
    print('PSFINISH 0.25', file=com)
    print('CCFIN 0.9', file=com)
    print('END', file=com)

  if formfactors == 'electron':
    atomsf = os.path.join(os.environ['CLIBD'], 'atomsf_electron.lib')
    acorn_command = shlex.split('acorn atomsf %s hklin %s xyzin %s hklout %s' % (atomsf, mtzin, xyzin, mtzout))
  else:
    acorn_command = shlex.split('acorn hklin %s xyzin %s hklout %s' % (mtzin, xyzin, mtzout))
  with open(acorn_script, 'r') as com, open(acorn_logfile, 'w') as acorn_log:
    acorn = subprocess.Popen(acorn_command, stdin=com, stdout=subprocess.PIPE, universal_newlines=True)
    per_cycle = []
    while True:
      line = acorn.stdout.readline()
      if line != '':
        acorn_log.write(line)
        if 'R-factor & Corr for medium E' in line:
          temp = line.split()
          cc = float(temp[10])
          cycle = int(temp[14])
          per_cycle.append({'cycle':cycle, 'CCs':cc})
        if ' Output phase set ' in line:
          cc_final = cc
      else:
        break
  return cc_final, per_cycle

def not_definitive(scores, acornCC_solved, acornCC_diff):
  if len(scores) == 0:
    return True
  elif sorted(scores, reverse=True)[0] >= acornCC_solved:
    return False
  elif -1 * (sorted(scores)[0] - sorted(scores, reverse=True)[0]) > acornCC_diff:
    return False
  else:
    return True

def test_solution(tempdir, scores, solution, test_all, acornCC_solved, acornCC_diff,
                  i, sigi, fp, sigfp, lowres, highres, solvent, formfactors):
  solution_id = solution['id']
  if i is not None and sigi is not None:
    have_intensities = True
  else:
    have_intensities = False
  mtzin = solution_id + '.mtz'
  log.debug('DEBUG %s %s' % (solution_id, scores))
  if os.path.isfile(solution_id + '.acorn.json'):
    # try and read results from json if restarting
    try:
      cc = read_results_json(solution_id + '.acorn.json')['acornCC']
    except KeyError:
      cc = None
  elif test_all or not_definitive(scores, acornCC_solved, acornCC_diff):
    if have_intensities:
       massage_data(solution_id, i, sigi)
    else:
      tidy_data(solution_id, fp, sigfp)
    extend_ecalc(solution_id=solution_id, i=i, sigi=sigi, fp=fp, sigfp=sigfp, highres=highres, tempdir=tempdir)
    mtzin = solution_id + '.aniso.ecalc.mtz'
    xyzin = solution_id + '.pdb'
    mtzout = solution_id + '.acorn.mtz'
    acorn_script = solution_id + '_acorn.com'
    acorn_logfile = solution_id + '_acorn.log'
    cc, cc_per_cycle = run_acorn(mtzin=mtzin, xyzin=xyzin, mtzout=mtzout,
                                 lowres=lowres, highres=highres, solvent=solvent, formfactors=formfactors,
                                 acorn_script=acorn_script, acorn_logfile=acorn_logfile)
    write_results_json(version=None, results_json=solution_id + '.acorn.json', solutions={'name':None, 'id':solution['id'], 'acornCC':cc, 'acornCC_per_cycle':cc_per_cycle})
  else:
    cc = None
  solution['acornCC'] = cc
  return solution

def test_solutions(json_file, xml_file, xmlroot, docid, output, formfactors, solutions, num_solutions,
                   test_all, acornCC_solved, acornCC_diff, i, sigi, fp, sigfp,
                   lowres, highres, solvent, tempdir, nproc):
  # this complicated arrangement ensures every time a solution is tested the function receives the latest scores
  scores = []
  if xml_file is not None or docid is not None: # i2/rvapi
    # set up initial running solutions to dump to xml
    for solution in solutions:
      if solution['number'] <= nproc:
        solution['acornCC'] = 'Running'
      else:
        solution['acornCC'] = '-'
    write_output({'solutions':solutions},
                json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, docid=docid, output=output)

  def callback(result):
    if result['acornCC'] is not None:
      log.info('Solution %d of %d ACORN CC = %0.5f\n' % (result['number'], num_solutions, result['acornCC']))
      scores.append(result['acornCC'])
    else:
      log.info('Solution %d of %d not tested as definite solution found\n' % (result['number'], num_solutions))
    if xml_file is not None or docid is not None: #i2/rvapi
      assert solutions[result['number'] - 1]['id'] == result['id']
      solutions[result['number'] - 1]['acornCC'] = result['acornCC']
      running = 0
      for solution in solutions:
        if solution['acornCC'] == 'Running':
          running += 1
        elif solution['acornCC'] is None:
          # defitive solution found
          break
        elif solution['acornCC'] == '-' and running < nproc:
          solution['acornCC'] = 'Running'
          running += 1
      write_output({'solutions':solutions},
                  json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, docid=docid, output=output)
    else:
      write_output({'result':result},
                  json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, docid=docid, output=output)
  solution_tester = setup_solution_tester(scores=scores,
                                          tempdir=tempdir,
                                          test_all=test_all,
                                          acornCC_solved=acornCC_solved,
                                          acornCC_diff=acornCC_diff,
                                          i=i, sigi=sigi, fp=fp, sigfp=sigfp,
                                          lowres=lowres, highres=highres,
                                          solvent=solvent,
                                          formfactors=formfactors)
  results = easy_mp.parallel_map(func=solution_tester, iterable=solutions, callback=callback, preserve_order=False, processes=nproc)
  return results
