import glob, os, subprocess
from fragon.data import massage_data, tidy_data, extend_ecalc
from fragon.utils import write_output, write_results_json, read_results_json
from libtbx import easy_mp

# writing acorn script
def write_acorn_script(scriptname, mtzin, xyzin, mtzout, lowres, highres, solvent):
  # we write a script to run acorn so need to make sure sh is in path
  for path in os.environ.get('PATH', '').split(':'):
    if os.path.exists(os.path.join(path, 'sh')) and not os.path.isdir(os.path.join(path, 'sh')):
      sh = os.path.join(path, 'sh')
  with open(scriptname, 'w') as acorn_script:
    print >> acorn_script, '#!' + sh
    print >> acorn_script, 'acorn hklin %s xyzin %s hklout %s <<EOF' % (mtzin, xyzin, mtzout)
    print >> acorn_script, 'labin E=E_ISO FP=F_ISO SIGFP=SIGF_ISO'
    print >> acorn_script, 'RESO %0.2f %0.2f' % (lowres,highres)
    if highres > 1.0:
      print >> acorn_script, 'EXTEND 50.0 1.0'
    else:
      print >> acorn_script, 'EXTEND 50.0 %0.2f' % highres
    print >> acorn_script, 'ESTRO 0.8'
    print >> acorn_script, 'POSI 1'
    print >> acorn_script, 'SOLV %s' % solvent
    print >> acorn_script, 'NTRY 1'
    print >> acorn_script, 'PSFINISH 0.25'
    print >> acorn_script, 'CCFIN 0.9'
    print >> acorn_script, 'END'
    print >> acorn_script, 'EOF'
  os.chmod(scriptname,0700)

def run_acorn(acorn_script, acorn_logfile):
  acorn_script = os.path.abspath(acorn_script)
  acorn_log = open(acorn_logfile, 'w')
  acorn = subprocess.Popen(acorn_script, stdout=subprocess.PIPE)
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

def test_solution(log, tempdir, scores, solution, test_all, acornCC_solved, acornCC_diff,
                  i, sigi, fp, sigfp, lowres, highres, solvent):
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
    extend_ecalc(solution_id=solution_id, i=i, sigi=sigi, fp=fp, sigfp=sigfp, highres=highres, tempdir=tempdir, log=log)
    mtzin = solution_id + '.aniso.ecalc.mtz'
    xyzin = solution_id + '.pdb'
    mtzout = solution_id + '.acorn.mtz'
    acorn_script = solution_id + '_acorn.sh'
    acorn_logfile = solution_id + '_acorn.log'
    write_acorn_script(scriptname=acorn_script, mtzin=mtzin, xyzin=xyzin, mtzout=mtzout,
                       lowres=lowres, highres=highres, solvent=solvent)
    cc, cc_per_cycle = run_acorn(acorn_script, acorn_logfile)
    write_results_json(version=None, results_json=solution_id + '.acorn.json', solutions={'name':None, 'id':solution['id'], 'acornCC':cc, 'acornCC_per_cycle':cc_per_cycle})
  else:
    cc = None
  solution['acornCC'] = cc
  return solution

def test_solutions(log, json_file, xml_file, xmlroot, output, solutions, num_solutions,
                   test_all, acornCC_solved, acornCC_diff, i, sigi, fp, sigfp,
                   lowres, highres, solvent, tempdir, nproc):
  # this complicated arrangement ensures every time a solution is tested the function receives the latest scores
  scores = []
  if xml_file is not None: # i2
    # set up initial running solutions to dump to xml
    for solution in solutions:
      if solution['number'] <= nproc:
        solution['acornCC'] = 'Running'
      else:
        solution['acornCC'] = '-'
    write_output({'solutions':solutions},
                json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

  def solution_tester(solution):
    # the next process is launched before the callback updates scores so we have to get the CC from a file
    if len(scores) == 0: # first callback
      for acorn_json in glob.glob('*.acorn.json'):
        try:
          scores.append(read_results_json(acorn_json)['acornCC'])
        except IOError: 
          log.debug('DEBUG Error reading file: %s' % acorn_json)
    elif not_definitive(scores, acornCC_solved, acornCC_diff): # avoid needless file reads
      try:
        scores.extend([ read_results_json(acorn_json)['acornCC'] for acorn_json in glob.glob('*.acorn.json') if read_results_json(acorn_json)['acornCC'] not in scores])
      except IOError:
        log.debug('DEBUG Error reading file: %s' % acorn_json)
    return test_solution(log=log, tempdir=tempdir, scores=scores, solution=solution,
                         test_all=test_all, acornCC_solved=acornCC_solved, acornCC_diff=acornCC_diff,
                         i=i, sigi=sigi, fp=fp, sigfp=sigfp,
                         lowres=lowres, highres=highres, solvent=solvent)

  def callback(result):
    if result['acornCC'] is not None:
      log.info('Solution %d of %d ACORN CC = %0.5f\n' % (result['number'], num_solutions, result['acornCC']))
      scores.append(result['acornCC'])
    else:
      log.info('Solution %d of %d not tested as definite solution found\n' % (result['number'], num_solutions))
    if xml_file is not None: #i2
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
                  json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    else:
      write_output({'result':result},
                  json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

  results = easy_mp.parallel_map(func=solution_tester, iterable=solutions, callback=callback, preserve_order=False, processes=nproc)
  return results
