from __future__ import division
import os
import sys
import errno
import shutil
import time
import logging
from fragon import data
from fragon import utils
from fragon import place
from fragon import score
from fragon import version
__version__ = version.__version__

def run():

  # start logging to stdout
  log = logging.getLogger()
  log.setLevel(logging.INFO)
  fmt = logging.Formatter('%(message)s')
  ch = logging.StreamHandler(sys.stdout)
  ch.setFormatter(fmt)
  log.addHandler(ch)

  args = utils.parse_command_line(sys.argv[1:])

  # are we run from CCP4i2?
  i2 = args.i2
  if i2:
    xml_file = 'program.xml'
    logfile = 'fragon.log'
    json_file = None
    output = None
    xmlroot = utils.write_output({'Fragon':__version__},xml_file=xml_file)

  debug = args.debug
  if debug:
    log.setLevel(logging.DEBUG)
    log.debug('DEBUG On')

   # restarting
  if args.results_json is not None:
    restart, solutions, search, scoring, logfile = utils.setup_restart(args.results_json, log)
    start_point = restart['start_point']
    run_dir = restart['run_dir']
    tempdir = restart['tempdir']
    name = restart['name']
    root = restart['root']
    mtzin = restart['mtzin']
  else:
    start_point = None
    run_dir = None
    # What to search for
    name, root, mtzin, search, logfile = utils.define_search(args)

  # user specified logfile
  if args.log is not None:
    logfile = os.path.abspath(logfile)
  # now logfile is defined
  logfile = os.path.abspath(logfile)
  if os.path.isfile(logfile):
    os.unlink(logfile)
  log.info('Writing logfile to %s' % logfile)
  fh = logging.FileHandler(logfile)
  fh.setFormatter(fmt)
  log.addHandler(fh)

  if not i2:
    json_file = logfile[:-4] + '.json'
    output = utils.write_output({'Fragon':__version__}, json_file=json_file, output={})
    xml_file = None
    xmlroot = None

  # print banner
  utils.print_header(log=log)

  # scoring
  if start_point is None:
    scoring = utils.define_scoring(args, log=log)

  #parallel
  if args.nproc is not None:
    nproc = args.nproc
  else:
    nproc = 1

  # initialise all timers to start time
  start_time = time.time()
  phaser_time = acorn_time = start_time
  log.info('Start time: %s\n' % str(time.asctime(time.localtime(time.time()))))
  utils.write_output({'start_time':str(time.asctime(time.localtime(time.time()))), 'name':name, 'root':root, 'nproc':nproc},
                    json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  log.info('Process ID: %s\n' % os.getpid())

  # set up versioned working directory
  if run_dir is None and not i2:
    run_dir, tempdir = utils.setup_run(mtzin=mtzin, seqin=search['seqin'], pdbin=search['pdbin'], pdbin_fixed=search['pdbin_fixed'])
  elif i2:
    run_dir = os.path.abspath(os.getcwd())
    os.mkdir('TEMP')
    tempdir = os.path.abspath('TEMP')
    if not args.helix:
      shutil.copy(os.path.join(os.environ['CCP4'], 'share', 'fragon', 'include', 'fragments', search['pdbin']), os.getcwd())

  log.debug('DEBUG start_point : %s  \n' % start_point)
  log.info('Working directory is: %s\n' % run_dir)
  os.chdir(run_dir)

  if args.name is None:
    log.info('Filename root set to %s override this with --name option\n' % name)
  else:
    log.info('Filename root set to %s\n' % name)

  # User supplied column labels
  if not any([args.I, args.SIGI, args.FP, args.SIGFP]): # all unset
    data_info = data.mtz_info(mtzin, log=log)
  else:
    i, sigi, fp, sigfp = utils.get_labels(args)
    data_info = data.mtz_info(mtzin, i=i, sigi=sigi, fp=fp, sigfp=sigfp, log=log)

  # check if we merged I(+)/I(-)
  if data_info['anomalous_merged']:
    mtzin = mtzin[:-4] + '_Imean.mtz'

  # prepare data for use by Phaser
  if data_info['i'] is not None and data_info['sigi'] is not None:
    if args.FP is None:
      i, sigi, fp, sigfp = data_info['i'], data_info['sigi'], None, None
      utils.write_output({'mtzin':mtzin, 'i':i, 'sigi':sigi},
                         json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    else:
      log.warning('*** Warning you have chosen to use structure factors ***\n')
      log.warning('It would be better to use intensities\n')
      i, sigi, fp, sigfp = None, None, data_info['fp'], data_info['sigfp']
  else:
    i, sigi, fp, sigfp = None, None, data_info['fp'], data_info['sigfp']
    utils.write_output({'mtzin':mtzin, 'f':fp, 'sigf':sigfp},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  log.debug('DEBUG i:%s sigi:%s fp:%s sigfp %s\n' % (i,sigi,fp,sigfp))
  data_logfile = root + '_prepare_data.log'
  log.info('Preparing data, logfile is: %s\n' % os.path.abspath(data_logfile))
  data_array = place.prepare_data(mtzin=mtzin, i=i, sigi=sigi, fp=fp, sigfp=sigfp, logfile=data_logfile, log=log)
  
  if start_point is None:
    # solvent content
    if search['solvent'] is not None:
      # over-ride Phaser calculation
      if search['solvent'] > 1:
        scoring['solvent'] = search['solvent']/100
        phaser_solvent = search['solvent']
      else:
        phaser_solvent = 100.00*search['solvent']
        scoring['solvent'] = search['solvent']
    else:
      cca_logfile = root + '_solvent.log'
      calc_z, calc_vm, z, vm = place.calculate_solvent(root=root, data=data_array,
                                                       seqin=search['seqin'],
                                                       ncs_copies=search['ncs_copies'],
                                                       highres=data_info['highres'],
                                                       logfile=cca_logfile,
                                                       log=log)
      scoring['solvent'] = round(1.0 - (1.232/vm), 2)
      if calc_z > 1:
        calc_solvent = round(1.0 - (1.232/calc_vm), 2)
        utils.write_output({'Warning':'Solvent', 'calc_z':calc_z, 'calc_solvent':calc_solvent},
                           json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

      # let Phaser calculate
      phaser_solvent = None
    utils.write_output({'ncs_copies':search['ncs_copies'], 'solvent':scoring['solvent']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    log.debug('DEBUG search solvent: %s solvent %0.2f\n' % (search['solvent'], scoring['solvent']))

  if start_point not in ['Scored solutions']:
    log.debug('DEBUG start_point : %s  \n' % start_point)
    if not scoring['test_all']:
      acornCC_solved, acornCC_diff = scoring['acornCC_solved'], scoring['acornCC_diff']
      utils.write_output({'acornCC_solved':acornCC_solved,'acornCC_diff':acornCC_diff},
                         json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      if acornCC_solved < 0.2:
        log.warning('************************************************ WARNING!!! ************************************************\n')
        log.warning('Fragon will stop testing solutions when ACORN CC above %0.5f - this is too low to guarantee success\n' % acornCC_solved)
        log.warning('************************************************************************************************************\n')
      else:
        if search['have_ensemble']:
          log.info('Testing a maximum of %d solution(s) for ensemble of %d models until ACORN CC is above %0.5f\n' %
                   (search['num_solutions'], search['num_models'],  scoring['acornCC_solved']))
        else:
          log.info('Testing a maximum of %d solution(s) until ACORN CC is above %0.5f\n' %
                   (search['num_solutions'], scoring['acornCC_solved']))
    else:
      utils.write_output({'acornCC_solved':'Test All'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      if search['have_ensemble']:
        log.info('Testing a maximum of %d solution(s) for ensemble of %d models\n' %
                 (search['num_solutions'], search['num_models']))
      else:
        log.info('Testing a maximum of %d solution(s)\n' % search['num_solutions'])
    utils.write_output({'num_solutions':search['num_solutions'], 'num_models':search['num_models']},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  if start_point is None:
    log.debug('DEBUG start_point : %s  \n' % start_point)
    if search['sgalternative']:
      log.info('Testing all possible space groups compatible with input space group %s\n' % data_info['sg'])
      utils.write_output({'sgalternative':'All'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    elif search['input_hand']:
      log.info('Testing only input space group %s\n' % data_info['sg'])
      utils.write_output({'sgalternative':'Input Only'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    elif search['sg_list'] is not None:
      log.info('Testing the following space group(s): %s\n' % search['sg_list'])
      utils.write_output({'sgalternative':search['sg_list']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    else :
      if data_info['enantiomorph'] is None:
        log.info('Testing input space group %s\n' % data_info['sg'])
      else:
        log.info('Testing input space group %s and enantiomorph (%s)\n' %(data_info['sg'], data_info['enantiomorph']))
      utils.write_output({'sgalternative':'Enantiomorph'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

    log.info('Data resolution is              %0.2f - %0.2f Angstrom' % (data_info['lowres'],  data_info['highres']))
    utils.write_output({'lowres':data_info['lowres'], 'highres':data_info['highres']},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if data_info['highres'] >= 1.71 and i2:
      utils.write_output({'Sorry':'Resolution'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      sys.exit('Sorry Fragon requires high resolution data.')
    log.info('Search resolution is            %0.2f - %0.2f Angstrom' % (search['lowres'] if search['lowres'] is not None else data_info['lowres'],
                                                                         search['highres'] if search['highres'] is not None else data_info['highres']))
    utils.write_output({'search_lowres':search['lowres'], 'search_highres':search['highres']},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if search['rot_samp'] is None:
      log.info('Rotation Search sampling        decided by Phaser')
    else:
      log.info('Rotation Search sampling is     %0.2f degrees' % search['rot_samp'])
      utils.write_output({'rot_samp':search['rot_samp']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if search['tra_samp'] is None:
      log.info('Translation Search sampling     decided by Phaser\n')
    else:
      log.info('Translation Search sampling is  %0.2f Angstrom\n' % search['tra_samp'])
      utils.write_output({'trans_samp':search['tra_samp']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if search['tfz_solved'] < 8.0:
      log.warning('TFZ for definitive solution set to %0.1f WARNING this is less than Phaser default!\n' % search['tfz_solved'])
    utils.write_output({'tfz_solved':search['tfz_solved']},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if search['search_down'] is not None:
      if search['search_down'] > 75:
        sys.exit('Search down percentage cannot be greater than 75 %')
      elif search['search_down'] > 15:
        log.info('Search down raised to add extra %d %% WARNING this is more than Phaser default (15 %%) and so will increase search time!\n' % search['search_down'])
      else:
        log.info('Search down lowered to only add extra %d %% WARNING this is less than Phaser default (15 %%) and so may decrease success!\n' % search['search_down'])
      utils.write_output({'search_down':search['search_down']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    if search['purge'] is not None:
      log.info('Purging solution list to keep only %d solutions after refinement of each copy\n' % search['purge'])
    utils.write_output({'purge':search['purge']}, json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  if start_point is None:
    log.debug('DEBUG start_point : %s  \n' % start_point)
    log.info('Attempting to place %s\n' % search['fragment'])
    if search['pdbin_fixed'] is not None:
      log.info('%s will be fixed in place\n' % search['pdbin_fixed'])
    if args.helix is not None:
      utils.write_output({'helix_length':args.helix},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    elif args.fragment is not None:
      utils.write_output({'fragment':search['pdbin']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    elif args.ensemble is not None:
      utils.write_output({'ensemble':search['pdbin']},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    utils.write_output({'copies':search['copies'], 'rms':search['rms']},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
                      
    # If highres and lowres are resolution limits of data set to 0 and 1000 to avoid rounding error truncating data in Phaser
    log.debug('DEBUG lowres %s search_lowres %s highres %s search_highres %s\n' % (data_info['lowres'], search['lowres'], data_info['highres'], search['highres']))
    log.debug('DEBUG xmlfile %s xmlroot %s' % (xml_file, xmlroot))
    
    phaser_solutions = place.place_fragment(root=root, xml_file=xml_file, xmlroot=xmlroot, log=log,                                            
                                            data=data_array,
                                            pdbin=search['pdbin'],
                                            copies=search['copies'],
                                            rms=search['rms'],
                                            pdbin_fixed=search['pdbin_fixed'],
                                            seqin=search['seqin'],
                                            ncs_copies=search['ncs_copies'],
                                            phaser_solvent=phaser_solvent,
                                            sgalternative=search['sgalternative'],
                                            sg_list=search['sg_list'],
                                            input_hand=search['input_hand'],
                                            tncs=search['tncs'],
                                            search_lowres=1000 if search['lowres'] is None else search['lowres'],
                                            search_highres=0 if search['highres'] is None else search['highres'],
                                            rot_peaks=search['rot_peaks'],
                                            rot_cluster_off=search['rot_cluster_off'],
                                            rot_samp=search['rot_samp'],
                                            tra_samp=search['tra_samp'],
                                            search_down=search['search_down'],
                                            tfz_solved=search['tfz_solved'],
                                            solutions=search['num_solutions'],
                                            purge=search['purge'],
                                            rescore_strands=search['rescore_strands'],
                                            rescore_residues=search['rescore_residues'],
                                            rescore_models=search['rescore_models'],
                                            rescore_all=search['rescore_all'],
                                            nproc=nproc)

    num_solutions = len(phaser_solutions) if len(phaser_solutions) < search['num_solutions'] else search['num_solutions']
    utils.write_output({'num_phaser_solutions':len(phaser_solutions)},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
                      
    # now set up solution list
    solutions = utils.get_solutions(phaser_solutions, root, num_solutions)
    # if Phaser failed to generate solutions this is an empty list
    if len(solutions) > 0:
      # write results to json for easy parsing
      results_json = root + '_phaser_solutions.json'
      utils.write_results_json(version=__version__, results_json=results_json, name=name, root=root, mtzin=mtzin, 
      stage='Phaser solutions', solutions=solutions, search=search, scoring=scoring)
      log.debug('DEBUG writing json after Phaser: %s\n' % solutions)
      if json_file is not None:
        # not for i2 as xml is updated when first solution is tested
        utils.write_output({'solutions':solutions},
                          json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
    else:
      utils.write_output({'Sorry':'No Solutions'},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      sys.exit('Phaser failed to find any solutions.')

    log.info('\nTime now: %s' % str( time.asctime(time.localtime(time.time()))))
    phaser_time = time.time()
    time_string = utils.print_time('Time for Phaser', phaser_time - start_time, log)
    utils.write_output({'phaser_time':time_string},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

  if start_point in [None, 'Phaser solutions']:
    log.debug('DEBUG start_point : %s\n' % start_point)
    if start_point == 'Phaser solutions':
      utils.write_output({'solutions':solutions},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      num_solutions = len(solutions) if len(solutions) < search['num_solutions'] else search['num_solutions']
    if not scoring['test_all']:
      log.info('Running ACORN density modification for up to %d solutions or until solution with CC above %0.5f found' % (num_solutions, scoring['acornCC_solved']))
    else:
      log.info('Running ACORN density modification for %d solutions' % num_solutions)

    # at atomic resolution don't stop until CC >= CCsolved
    acornCC_diff = acornCC_solved if data_info['highres'] <= 1.20 else scoring['acornCC_diff']
    log.debug('DEBUG acornCC_diff: %0.4f search acornCC_diff: %0.4f data highres: %0.2f\n' % (acornCC_diff, scoring['acornCC_diff'], data_info['highres']))

    solutions = score.test_solutions(log=log, json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output,
                                     solutions=solutions, 
                                     num_solutions=num_solutions,
                                     test_all=scoring['test_all'], 
                                     acornCC_solved=scoring['acornCC_solved'], 
                                     acornCC_diff=acornCC_diff,
                                     i=i, sigi=sigi, fp=fp, sigfp=sigfp,
                                     lowres=data_info['lowres'], 
                                     highres=data_info['highres'], 
                                     solvent=scoring['solvent'],
                                     tempdir=tempdir,  
                                     nproc=nproc)

    # write results to json for easy parsing
    results_json = root + '_scored_solutions.json'
    utils.write_results_json(version=__version__, results_json=results_json, name=name, root=root, mtzin=mtzin, 
                            stage='Scored solutions', solutions=solutions, search=search, scoring=scoring)
    log.debug('DEBUG writing json after ACORN: %s\n' % solutions)

    log.info('\nTime now: %s\n' % str(time.asctime(time.localtime(time.time()))))
    acorn_time = time.time()
    time_string = utils.print_time('Time for ACORN density modification', acorn_time - phaser_time, log)
    utils.write_output({'acorn_time':time_string},
                      json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)

  if start_point in [None, 'Phaser solutions', 'Scored solutions']:
    log.debug('DEBUG start_point : %s\n' % start_point)
    if start_point == 'Scored solutions':
      utils.write_output({'solutions':solutions},
                        json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
      num_solutions = len(solutions) if len(solutions) < search['num_solutions'] else search['num_solutions']

    # count how many we tested and get best
    solutions_tested = utils.count_tested(solutions, 'acornCC')
    cc_best = sorted(solutions, key=lambda r: r['acornCC'], reverse=True)[0]['acornCC']
    best_solution_id = sorted(solutions, key=lambda r: r['acornCC'], reverse=True)[0]['id']
    
  log.info('\n\nAll done')
  log.info('Time now: %s\n' % str(time.asctime(time.localtime(time.time()))))
  end_time = time.time()
  time_string = utils.print_time('Total time taken', end_time - start_time, log)
  utils.write_output({'total_time':time_string, 'finish_time':str(time.asctime(time.localtime(time.time())))},
                    json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  log.info('Best solution %s.pdb , CC = %0.5f\n' % (best_solution_id, cc_best))
  utils.write_output({'best_solution_id':best_solution_id, 'cc_best':cc_best},
                    json_file=json_file, xml_file=xml_file, xmlroot=xmlroot, output=output)
  if i2:
    if i is not None and sigi is not None:
      mtzin = os.path.join(tempdir, best_solution_id + '.fiso.mtz')
    else:
      mtzin = os.path.join(tempdir, best_solution_id + '.aniso.mtz')
    acorn_mtz = best_solution_id + '.acorn.mtz'
    phaser_pdb = best_solution_id + '.pdb'
    shutil.copy(phaser_pdb, name + '_xyzout_fragon.pdb')
    data.minimtz_output(name, mtzin, acorn_mtz)
  else:
    if cc_best > 0.20000:
      utils.write_output_files(name=name, best_solution_id=best_solution_id, mtzin=mtzin, log=log)
  if not debug:
    for tempfile in os.listdir(tempdir):
      file_path = os.path.join(tempdir, tempfile)
      if os.path.isfile(file_path):
        os.unlink(file_path)
  utils.print_refs(log)
if __name__ == '__main__':
  run()
