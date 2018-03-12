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
import sys
import shutil
import logging
from iotbx import mtz
from cctbx.miller import array_info
# for anisotropy correction
import phaser
import clipper

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
  output = phaser.Output()
  aniso_log = open(logfile, 'w')
  output_object = CallbackObject(None, None)
  output.setPhenixPackageCallback(output_object)
  output.set_file_object(aniso_log)
  data = phaser.runMR_DAT(input, output)
  mtzout = solution_id + '.aniso'
  input = phaser.InputANO()
  input.setSPAC_HALL(data.getSpaceGroupHall())
  input.setCELL6(data.getUnitCell())
  input.setREFL_DATA(data.getDATA())
  input.setHKLI(mtzin)
  input.setROOT(mtzout)
  aniso = phaser.runANO(input, output)
  aniso_log.close()
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
    highres = 1.0
  if have_intensities:
    ipfile = solution_id + '.french-wilson.mtz'
    ipfile_fiso = solution_id + '.fiso.mtz'
  else:
    ipfile = solution_id + '.aniso.mtz'
  # set up Clipper data objects
  mtzin = clipper.CCP4MTZfile()
  mtzin.open_read(ipfile)
  sg = mtzin.spacegroup()
  cell = mtzin.cell()
  reso = clipper.Resolution(highres)
  hkls = clipper.HKL_info(sg, cell, reso, True)
  fsig = clipper.HKL_data_F_sigF_float(hkls)
  fsigiso = clipper.HKL_data_F_sigF_float(hkls)
  # read data
  mtzin.import_hkl_data(fsig, '*/*/[F, SIGF]')
  if not have_intensities:
    mtzin.import_hkl_data(fsigiso, '*/*/[F_ISO, SIGF_ISO]')
    tempfiles.append(ipfile)
  mtzin.close_read()
  if have_intensities:
    mtzin_fiso = clipper.CCP4MTZfile()
    mtzin_fiso.open_read(ipfile_fiso)
    mtzin_fiso.import_hkl_data(fsigiso, '*/*/[F_ISO, SIGF_ISO]')
    mtzin_fiso.close_read()
    tempfiles.append(ipfile_fiso)
  # Calculate Es from F_ISO
  esig = clipper.HKL_data_E_sigE_float(hkls)
  esig.compute_from_fsigf(fsigiso)
  # now calculate scaling
  nparm = 12 # as in cecalc
  initial_params = clipper.DoubleVector(nparm, 1.0)
  basis_f = clipper.BasisFn_spline(esig, nparm, 2.0)
  target_f = clipper.TargetFn_scaleEsq_E_sigE(esig)
  escale = clipper.ResolutionFn(hkls, basis_f, target_f, initial_params)
  # apply scaling
  esig.scaleBySqrtResolution(escale)
  # write output
  mtzout = clipper.CCP4MTZfile()
  opfile = solution_id+'.aniso.ecalc.mtz'
  mtzout.open_write(opfile)
  mtzout.export_hkl_info(fsig.hkl_info())
  mtzout.export_hkl_data(fsig, '*/*/[F, SIGF]')
  mtzout.export_hkl_data(fsigiso, '*/*/[F_ISO, SIGF_ISO]')
  mtzout.export_hkl_data(esig, '*/*/[E_ISO, SIGE_ISO]')
  mtzout.close_write()
# cleanup
  for file in tempfiles:
    shutil.move(file,tempdir)


def mtz_output(mtzin, acorn_mtz, mtzout):
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
  info = fiso.info()
  fiso = fiso.map_to_asu().set_info(info)
  fiso, fom = fiso.common_sets(other=fom.map_to_asu())
  f = fiso * fom
  f, phi = f.common_sets(other=phi.map_to_asu())
  map_coeffs = f.phase_transfer(phi, deg=True).customized_copy(sigmas=None).set_info(info)
  info = array_info(labels=['FWT', 'PHWT'])
  map_coeffs.set_info(info)
  combined_mtz.add_miller_array(map_coeffs, column_root_label='FWT')
  combined_mtz = combined_mtz.mtz_object()
  for column in combined_mtz.columns():
    if column.label()=='PHIFWT':
      column.set_label('PHWT')

  combined_mtz.write(mtzout)

def minimtz_output(name, mtzin, acorn_mtz):
  mtzref = clipper.CCP4MTZfile()
  mtzphifom = clipper.CCP4MTZfile()
  mtzfphi = clipper.CCP4MTZfile()
  hkls_ref = clipper.HKL_info()
  # open the file from Phaser that (may have been F-W and) has been corrected for anisotropy
  mtzref.open_read( mtzin )
  sg = mtzref.spacegroup()
  cell = mtzref.cell()
  mtzref.import_hkl_info(hkls_ref, False)
  fsig_ref = clipper.HKL_data_F_sigF_float(hkls_ref)
  mtzref.import_hkl_data(fsig_ref, '*/*/[F_ISO,SIGF_ISO]')
  mtzref.close_read()

  mtzacorn = clipper.CCP4MTZfile()
  #  read f_iso, and phi/fom from acorn
  mtzacorn.open_read( acorn_mtz )
  hkls = clipper.HKL_info()
  mtzacorn.import_hkl_info(hkls)
  reso = mtzacorn.resolution()
  fsigiso = clipper.HKL_data_F_sigF_float(hkls)
  phifom = clipper.HKL_data_Phi_fom_float(hkls)
  mtzacorn.import_hkl_data(fsigiso, '*/*/[F_ISO,SIGF_ISO]')
  mtzacorn.import_hkl_data(phifom, '*/*/[PHIOUT,WTOUT]')
  mtzacorn.close_read()

  #  create a list of hkls present in both unextended and extended mtz
  matched = clipper.HKLVector()
  hkls_matched = clipper.HKL_info(sg, cell, reso, False)
  clipper.PopulateMatchesF_sigF_float(fsig_ref, fsigiso, matched)
  hkls_matched.add_hkl_list(matched)
  fsigiso_out = clipper.HKL_data_F_sigF_float(hkls_matched)
  phifom_out = clipper.HKL_data_Phi_fom_float(hkls_matched)
  fphi_out = clipper.HKL_data_F_phi_float(hkls_matched)

  # and copy data for non-missing reflections
  clipper.CopyIfF_sigFRefNotMissingF_sigF_float(fsigiso, fsigiso_out, fsig_ref)
  clipper.CopyIfF_sigFRefNotMissingPhi_fom_float(phifom, phifom_out, fsig_ref)

  # calculate map coefficients
  fphi_out.compute_from_fsigf_phifom(fsigiso_out, phifom_out)

  # write out mini-MTZs
  acorn_phases = name + '_phsout_fragon.mtz'
  map_coeffs = name + '_fphiout_fragon.mtz'

  mtzphifom.open_write(acorn_phases)
  mtzphifom.export_hkl_info(phifom_out.hkl_info())
  mtzphifom.export_hkl_data(phifom_out, '*/*/[PHI,FOM]')
  mtzphifom.close_write()
  mtzfphi.open_write(map_coeffs)
  mtzfphi.export_hkl_info(fphi_out.hkl_info())
  mtzfphi.export_hkl_data(fphi_out, '*/*/[F,PHI]')
  mtzfphi.close_write()
