"""
    Fragon data.py 

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
from __future__ import print_function
import sys
import shutil
import logging
from iotbx import mtz
from cctbx.miller import array_info
# for anisotropy correction
import phaser

log = logging.getLogger(__name__)

def mtz_info(mtzin, i=None, sigi=None, fp=None, sigfp=None):
  mtz_file = mtz.object(mtzin)
  mtz_low_res = mtz_file.max_min_resolution()[0]
  mtz_high_res = mtz_file.max_min_resolution()[1]
  cell = []
  for parameter in range(6):
    cell.append(mtz_file.crystals()[0].unit_cell_parameters()[parameter])
  sg = mtz_file.space_group().type().lookup_symbol()
  if mtz_file.space_group().type().is_enantiomorphic():
    enantiomorph = mtz_file.space_group_info().change_hand().type().lookup_symbol()
  else:
    enantiomorph = None
  # convert cctbx R 3 :H etc
  if sg == 'R 3 2 :H':
    sg = 'H 3 2'
  elif sg == 'R 3 :H':
    sg = 'H 3'
  results = {'cell':cell, 'sg':sg, 'enantiomorph':enantiomorph}
  freeR_flag = None
  i_found = False
  sigi_found = False
  fp_found = False
  sigfp_found = False
  anomalous_merged = False
  for column in mtz_file.column_labels():
    if column in ['FreeR_flag', 'FREE', 'Free']:
      freeR_flag = column
    elif column in ['I', 'IMEAN'] or column[0:6] == 'IMEAN_' or column[0:2] == 'I_' and column[-3:] not in ['(+)', '(-)']:
      if i_found:
        sys.exit('Multiple I coulumns in input file, please assign input labels with --I and --SIGI')
      else:
        i = column
        i_found = True
    elif column in ['SIGI', 'SIGIMEAN'] or column[0:9] == 'SIGIMEAN_' or column[0:5] == 'SIGI_' and column[-3:] not in ['(+)', '(-)']:
      if sigi_found:
        sys.exit('Multiple SIGI columns in input file, please assign input labels with --I and --SIGI')
      else:
        sigi = column
        sigi_found = True
    elif column in ['F', 'FP'] or column[0:2] == 'F_' and column[-3:] not in ['(+)', '(-)']:
      if fp_found:
        sys.exit('Multiple FP coulumns in input file, please assign input labels with --FP and --SIGFP')
      else:
        fp = column
        fp_found = True
    elif column in ['SIGF', 'SIGFP'] or column[0:5] == 'SIGF_' and column[-3:] not in ['(+)', '(-)']:
      if sigfp_found:
        sys.exit('Multiple SIGFP columns in input file, please assign input labels with --FP and --SIGFP')
      else:
        sigfp = column
        sigfp_found = True
  # if we can't find IMEAN try and find I(+) and I(-), warn (unless i2) and merge.
  if i is None or sigi is None:
    miller_arrays = mtz_file.as_miller_arrays()
    iplus_iminus_found = False
    for array in miller_arrays:
      if str(array.info().labels[0]) == freeR_flag:
        freer = array
      elif array.is_xray_intensity_array():
        if iplus_iminus_found:
          sys.exit('Multiple I+/I- columns in input file, unable to proceed')
        else:
          iobs = array
          iplus_iminus_found = True
    if iplus_iminus_found:
      anomalous_merged = True
      log.info('Using I and SIGI from merged I(+) and I(-) in input file\n')
      iobs = iobs.average_bijvoet_mates()
      if freeR_flag is not None:
        iobs_mtz = freer.as_mtz_dataset(column_root_label='FreeR_flag')
        iobs_mtz.add_miller_array(iobs, column_root_label='I')
      else:
        iobs_mtz = iobs.as_mtz_dataset(column_root_label='I')
      i = 'I'
      sigi = 'SIGI'
      mtzout = mtzin[:-4] + '_Imean.mtz'
      iobs_mtz.mtz_object().write(mtzout)
  results.update({'anomalous_merged':anomalous_merged})
  if not any([i, sigi, fp, sigfp]): # all unset
    sys.exit('Could not assign input labels')
  elif i is not None and sigi is not None:
    log.info('Using I=%s SIGI=%s from input file %s\n' % (i, sigi, mtzin))
    label = "['" + i + "', '" + sigi + "']"
  elif fp is not None and sigfp is not None:
    log.warning('*** Warning could not find intensities in input file ***\n')
    log.info('Using FP=%s SIGFP=%s from input file %s\n' % (fp, sigfp, mtzin))
    label = "['" + fp + "', '" + sigfp + "']"
  results.update({'i':i, 'sigi':sigi, 'fp':fp, 'sigfp':sigfp, 'freeR_flag':freeR_flag})
  for array in mtz_file.as_miller_arrays():
    if str(array.info().labels) == label:
      data_low_res = array.d_max_min()[0]
      data_high_res = array.d_max_min()[1]
  if anomalous_merged:
    data_low_res = iobs.d_max_min()[0]
    data_high_res = iobs.d_max_min()[1]
  mtz_low_res = round(mtz_low_res, 2)
  mtz_high_res = round(mtz_high_res, 2)
  data_low_res = round(data_low_res, 2)
  data_high_res = round(data_high_res, 2)
  if data_low_res != mtz_low_res:
    results.update({'lowres':data_low_res, 'highres':data_high_res})
  elif data_high_res != mtz_high_res:
    results.update({'lowres':data_low_res, 'highres':data_high_res})
  else:
    results.update({'lowres':mtz_low_res, 'highres':mtz_high_res})

  return results

# use cctbx french_wilson() to massage the intensities

def massage_data(solution_id, i, sigi):
  if i[-4:] == '_ISO':
    mtzin = solution_id + '.aniso.mtz'
  else:
    mtzin = solution_id + '.mtz'
  label = "['" + i + "', '" + sigi + "']"
  mtz_file = mtz.object(mtzin)
  miller_arrays = mtz_file.as_miller_arrays()
  for array in miller_arrays:
    if str(array.info().labels) == label:
      iobs = array
  fobs = iobs.french_wilson()
  if i[-4:] == '_ISO':
    fobs_mtz = fobs.as_mtz_dataset(column_root_label='F_ISO')
    mtzout = solution_id + '.fiso.mtz'
  else:
    fobs_mtz = fobs.as_mtz_dataset(column_root_label='F')
    mtzout = solution_id + '.french-wilson.mtz'
  fobs_mtz.mtz_object().write(mtzout)

# this was originally to workaround (now fixed) Phaser anisotropy correction bug but also renames FP/SIGFP

def tidy_data(solution_id, fp, sigfp):
  mtzin = solution_id + '.mtz'
  label = "['" + fp + "', '" + sigfp + "']"
  mtz_file = mtz.object(mtzin)
  miller_arrays = mtz_file.as_miller_arrays()
  for array in miller_arrays:
    if str(array.info().labels) == label:
      fobs = array
  fobs_mtz = fobs.as_mtz_dataset(column_root_label='F')
  mtzout = solution_id + '.F.mtz'
  fobs_mtz.mtz_object().write(mtzout)

def aniso_correct(solution_id, mtzin, i, sigi, fp, sigfp, logfile):
  from place import CallbackObject
  input = phaser.InputMR_DAT()
  input.setHKLI(mtzin)
  if i is not None and sigi is not None:
    input.setLABI_I_SIGI(i, sigi)
  else:
    input.setLABI_F_SIGF(fp, sigfp)
  input.setMUTE(True)
  data = phaser.runMR_DAT(input)
  with open(logfile, 'w') as aniso_log:
    print(data.logfile(), file=aniso_log)
  mtzout = solution_id + '.aniso'
  input = phaser.InputANO()
  input.setSPAC_HALL(data.getSpaceGroupHall())
  input.setCELL6(data.getUnitCell())
  input.setREFL_DATA(data.getDATA())
  input.setHKLI(mtzin)
  input.setROOT(mtzout)
  input.setMUTE(True)
  aniso = phaser.runANO(input)
  with open(logfile, 'a') as aniso_log:
    print(aniso.logfile(), file=aniso_log)
  if data.Success():
    success = True
  else:
    log.critical('Job exit status FAILURE')
    log.critical('%s ERROR : %s' % (data.ErrorName(), data.ErrorMessage()))

# This does anisotropy correction with Phaser, extends to 1.0 A (or completes to highres of input if higher) and calcs Es

def extend_ecalc(solution_id, i, sigi, fp, sigfp, highres, tempdir):
  if i is not None and sigi is not None:
    have_intensities = True
  else:
    have_intensities = False
  tempfiles = []
  # aniso correct
  if have_intensities:
    mtzin = solution_id + '.mtz'
  else:
    mtzin = solution_id + '.F.mtz'
  # now F and SIGF are defined
  fp = 'F'
  sigfp = 'SIGF'
  aniso_logfile = solution_id + '_anisotropy_correction.log'
  tempfiles.append(aniso_logfile)
  aniso_correct(solution_id, mtzin, i, sigi, fp, sigfp, aniso_logfile)
  if have_intensities:
    i_iso = i + '_ISO'
    sigi_iso = sigi + '_ISO'
    massage_data(solution_id, i_iso, sigi_iso)
    tempfiles.append(solution_id + '.aniso.mtz')
  else:
    tempfiles.append(mtzin)
  # extend
  if highres > 1.0:
    highres = 1.00000
  if have_intensities:
    ipfile = solution_id + '.french-wilson.mtz'
    ipfile_fiso = solution_id + '.fiso.mtz'
  else:
    ipfile = solution_id + '.aniso.mtz'
  mtz_file = mtz.object(ipfile)
  miller_arrays = mtz_file.as_miller_arrays()
  for array in miller_arrays:
    if str(array.info().labels) == "['F', 'SIGF']":
      fobs = array 
    elif str(array.info().labels) == "['F_ISO', 'SIGF_ISO']":
      fiso = array
  if have_intensities:
    mtz_file = mtz.object(ipfile_fiso)
    miller_arrays = mtz_file.as_miller_arrays()
    for array in miller_arrays:
      if str(array.info().labels) == "['F_ISO', 'SIGF_ISO']":
        fiso = array
    tempfiles.append(ipfile_fiso)
  else:
    tempfiles.append(ipfile)
  
  # Extend. Extended reflections have F=0.0 SIGF=0.0 which means the iotbx mtz reader will treat them as missing.
  fobs = fobs.complete_array(d_min_tolerance=1e-06, d_min=highres, d_max=None, new_data_value=0.0, new_sigmas_value=0.0)    
  # Calculate Es from F_ISO
  fiso.setup_binner(reflections_per_bin=200)
  eiso = fiso.quasi_normalize_structure_factors().quasi_normalized_as_normalized()
  # write output
  mtzout = fobs.as_mtz_dataset(column_root_label='F')
  mtzout.add_miller_array(fiso, column_root_label='F_ISO')
  mtzout.add_miller_array(eiso, column_root_label='E_ISO')
  opfile = solution_id+'.aniso.ecalc.mtz'
  mtzout.mtz_object().write(opfile)
# cleanup
  for file in tempfiles:
    shutil.move(file,tempdir)


def mtz_output(mtzin, acorn_mtz, mtzout=None, minimtz_phifom=None, minimtz_fphi=None):
  freer = None
  freer_labels = ['FreeR_flag', 'FREE', 'Free']
  mtz_file = mtz.object(mtzin)
  miller_arrays = mtz_file.as_miller_arrays()
  for array in miller_arrays:
    if str(array.info().labels[0]) in freer_labels:
      freer = array
  acorn_mtz_file = mtz.object(acorn_mtz)
  acorn_miller_arrays = acorn_mtz_file.as_miller_arrays()
  for array in acorn_miller_arrays:
    if str(array.info().labels) == "['F', 'SIGF']":
      fobs = array
    elif str(array.info().labels) == "['F_ISO', 'SIGF_ISO']":
      fiso = array
    elif str(array.info().labels) ==  "['EOEXT', 'PHIOUT']":
       phi = array.customized_copy(data=array.phases(deg=True).data())
    # catch inverted files
    elif str(array.info().labels) ==  "['PHIOUT']":
       phi = array
    elif str(array.info().labels) == "['WTOUT']":
      fom = array
  if minimtz_phifom is None and minimtz_fphi is None:
    if freer is not None:
      # catch situations where input and solution spacegroups differ
      if freer.space_group() != fobs.space_group():
        freer = freer.customized_copy(space_group_info=fobs.space_group_info())
      combined_mtz = freer.as_mtz_dataset(column_root_label='FreeR_flag')
      combined_mtz.add_miller_array(fobs, column_root_label='F')
    else:
      combined_mtz = fobs.as_mtz_dataset(column_root_label='F')
    combined_mtz.add_miller_array(phi, column_root_label='PHI', column_types='P')
    combined_mtz.add_miller_array(fom, column_root_label='FOM', column_types='W')
  else:
    fobs = fobs.map_to_asu()
    fobs, phi_fobs = fobs.common_sets(other=phi.map_to_asu())
    fobs, fom_fobs = fobs.common_sets(other=fom.map_to_asu())
    phifom = phi_fobs.as_mtz_dataset(column_root_label='PHI', column_types='P')
    phifom.add_miller_array(fom_fobs, column_root_label='FOM', column_types='W')
    phifom.mtz_object().write(minimtz_phifom)
  info = fiso.info()
  fiso = fiso.map_to_asu().set_info(info)
  fiso, fom = fiso.common_sets(other=fom.map_to_asu())
  f = fiso * fom
  f, phi = f.common_sets(other=phi.map_to_asu())
  map_coeffs = f.phase_transfer(phi, deg=True).customized_copy(sigmas=None).set_info(info)
  if minimtz_phifom is None and minimtz_fphi is None:
    info = array_info(labels=['FWT', 'PHWT'])
    map_coeffs.set_info(info)
    combined_mtz.add_miller_array(map_coeffs, column_root_label='FWT')
    combined_mtz = combined_mtz.mtz_object()
    for column in combined_mtz.columns():
      if column.label()=='PHIFWT':
        column.set_label('PHWT')
    combined_mtz.write(mtzout)
  else:
    fphi = map_coeffs.as_mtz_dataset(column_root_label='F')
    fphi_mtz = fphi.mtz_object()
    for column in fphi_mtz.columns():
      if column.label()=='PHIF':
        column.set_label('PHI')
    fphi_mtz.write(minimtz_fphi)
