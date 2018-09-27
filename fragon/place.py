"""
    Fragon place.py

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
import os
import logging
import phaser
from iotbx import pdb
from fragon.utils import write_output

log = logging.getLogger(__name__)

class CallbackObject(object):
  def __init__(self, xml_file=None, xmlroot=None, docid=None, output=None):
    self.xml_file = xml_file
    self.xmlroot = xmlroot
    self.docid = docid
    self.output = output

  def startProgressBar (self, label, size):
    if self.xml_file is None and self.xmlroot is None and self.docid is None:
      pass # for command-line mode
    else:
      from fragon.utils import write_output
      if label != 'Generating Statistics':
        callback = ['progress', label]
        write_output({'callback':callback}, json_file=None, xml_file=self.xml_file, xmlroot=self.xmlroot, docid=self.docid, output=self.output)
  def incrementProgressBar (self):
    pass
  def endProgressBar (self):
    pass
  def warn (self, message):
    pass
  def loggraph (self, title, data):
    pass
  def call_back (self, message, data):
    if self.xml_file is None and self.xmlroot is None and self.docid is None:
      pass # for command-line mode
    else:
      from fragon.utils import write_output
      if message == 'summary':
        if data.startswith('** New Best LLG'):
          llg = data.split('=')[1].split()[0]
          callback = ['Best LLG', llg]
          write_output({'callback':callback}, json_file=None, xml_file=self.xml_file, xmlroot=self.xmlroot, docid=self.docid, output=self.output)
      elif message == 'current best solution':
        for item in data.split():
          if item[0:5] == 'TFZ==':
            tfz = item.split('==')[1]
            callback = ['Best TFZ', tfz]
            write_output({'callback':callback}, json_file=None, xml_file=self.xml_file, xmlroot=self.xmlroot, docid=self.docid, output=self.output)


def prepare_data(mtzin, i, sigi, fp, sigfp, logfile):
  input = phaser.InputMR_DAT()
  input.setHKLI(mtzin)
  if i is not None and sigi is not None:
    input.setLABI_I_SIGI(i, sigi)
  else:
    input.setLABI_F_SIGF(fp, sigfp)
  input.setMUTE(True)
  data = phaser.runMR_DAT(input)
  with open(logfile, 'w') as data_log:
    print(data.logfile(), file=data_log)
  if data.Success():
    return data
  else:
    log.critical('Job exit status FAILURE')
    log.critical('%s ERROR : %s' % (data.ErrorName(), data.ErrorMessage()))

def calculate_solvent(root, data, seqin, ncs_copies, highres, logfile):
  input = phaser.InputCCA()
  input.setSPAC_HALL(data.getSpaceGroupHall())
  input.setCELL6(data.getUnitCell())
  input.addCOMP_PROT_SEQ_NUM(seqin, ncs_copies)
  input.setRESO(highres, 10000)
  input.setMUTE(True)
  cca = phaser.runCCA(input)
  with open(logfile, 'w') as cca_log:
    print(data.logfile(), file=cca_log)
  # check if solvent content is higher than average
  if cca.getBestZ() > 1:
    log.warning ('*** Warning solvent content estimated by Phaser is %0.2f' % (1.0 - 1.232/cca.getBestVM()))
    if ncs_copies > 1:
      log.warning('    solvent content calculated from %d copies of sequence is %0.2f\n' % (ncs_copies, (1.0 - 1.232/cca.getVM()[0])))
    else:
      log.warning('    solvent content calculated from 1 copy of sequence is %0.2f\n' % (1.0 - 1.232/cca.getVM()[0]))
  elif ncs_copies > 1:
    log.info('Solvent content calculated from %d copies of sequence is %0.2f\n' % (ncs_copies, (1.0 - 1.232/cca.getVM()[0])))
  else:
    log.info('Solvent content calculated from 1 copy of sequence is %0.2f\n' % (1.0 - 1.232/cca.getVM()[0]))

  return cca.getBestZ(), cca.getBestVM(), cca.getZ()[0], cca.getVM()[0]

def run_phaser(root, xml_file, xmlroot, docid, output, logfile, data, pdbin, copies, rms, pdbin_fixed,
               seqin, ncs_copies, phaser_solvent, sgalternative, sg_list, input_hand, tncs,
               search_lowres, search_highres, rot_peaks, rot_cluster_off,
               rot_samp, tra_samp, search_down, tfz_solved, solutions, purge,
               rescore_strands, rescore_residues, rescore_models, nproc):
  input = phaser.InputMR_AUTO()
  input.setRESO(search_highres, search_lowres)
  # test OpenMP
  if nproc > 1:
    input.setJOBS(nproc)
  else:
    input.setJOBS(1)
  if rot_samp is not None:
    input.setSAMP_ROTA(rot_samp)
  if tra_samp is not None:
    input.setSAMP_TRAN(tra_samp)
  input.setSPAC_HALL(data.getSpaceGroupHall())
  input.setCELL6(data.getUnitCell())
  input.setREFL_DATA(data.getDATA())
  input.addENSE_PDB_RMS(root, pdbin, rms)
  input.setROOT(root)
  # fixed solution
  if pdbin_fixed is not None:
    fixed_root = pdbin_fixed[:-4]
    input.addENSE_PDB_RMS(fixed_root, pdbin_fixed, rms)
    input.addSOLU_6DIM_ENSE(fixed_root, [0.0, 0.0, 0.0], False, [0.0, 0.0, 0.0], 0.0, False, False, False, 1.0, 1.0)
  if phaser_solvent is not None:
    input.setCOMP_SOLV()
    input.setCOMP_PERC(phaser_solvent)
  else:
    input.addCOMP_PROT_SEQ_NUM(seqin, ncs_copies)
  input.addSEAR_ENSE_NUM(root, copies)
  if rot_peaks is not None:
    input.setPEAK_ROTA_SELE('NUMBER')
    input.setPEAK_ROTA_CUTO(rot_peaks)
  if rot_cluster_off:
    input.setPEAK_ROTA_CLUS(False)
  if sgalternative:
    input.setSGAL_SELE('ALL')
  elif input_hand:
    input.setSGAL_SELE('NONE')
  elif sg_list is not None:
    input.setSGAL_SELE('LIST')
    for sg in sg_list.split(','):
      input.addSGAL_TEST(sg)
  else:
    input.setSGAL_SELE('HAND')
  input.setTOPF(solutions)
  if rescore_strands or rescore_residues or rescore_models:
    # we don't need to write out the files here
    input.setXYZO(False)
    input.setHKLO(False)
  input.setXYZO_ENSE(False)

  # turn off tNCS correction if required
  if tncs:
    # to turn restore default MR_AUTO behaviour needs this
    input.setTNCS_PATT_PERC(99.9)
  # in desparate cases search further down the RF and TF list
  if search_down is not None:
    input.setSEAR_DOWN_PERC(search_down)
    input.setPURG_TRAN_PERC(75.00 - search_down)
  # turn purging on if requested for multicopy searches
  if purge is not None:
    input.setPURG_RNP_ENAB(True)
    input.setPURG_RNP_NUMB(purge)
   # turn it off otherwise (so more solutions are output)
  else:
    input.setPURG_RNP_ENAB(False)
  # these are not default but fragments should never clash.
  input.setPACK_KEEP_HIGH_TFZ(False)
  input.setSEAR_PRUN(False)
  input.setTRAN_PACK_CUTO(0.0)
  input.setPACK_CUTO(0.0)
  # end non defaults
  input.setZSCO_USE(True)
  input.setZSCO_SOLV(tfz_solved)
  if xml_file is not None and xmlroot is not None:
    input.setOUTP_LEVE('LOGFILE') #i2
  elif docid is not None:
    input.setOUTP_LEVE('LOGFILE') #rvapi
  else:
    input.setOUTP_LEVE('SUMMARY') #command line
  phaser_output = phaser.Output()
  phaser_log = open(logfile, 'w')
  output_object = CallbackObject(xml_file, xmlroot, docid, output)
  phaser_output.setPhenixPackageCallback(output_object)
  phaser_output.set_file_object(phaser_log)
  mr = phaser.runMR_AUTO(input, phaser_output)
  phaser_log.close()
  if mr.Success():
    if mr.foundSolutions() :
      log.info('\nPhaser has found %d solution(s)' % mr.numSolutions())
      if not rescore_strands or rescore_residues:
        log.info('Top solution PDB file = %s' % mr.getTopPdbFile())
      log.info('Top solution LLG = %f' % mr.getTopLLG())
      return mr.getDotSol()
    else:
      log.info('\nPhaser failed to find any solutions')
      return 'No solutions'
  else:
      log.info('\nPhaser failed to find any solutions')
      return 'No solutions'

def split_chains(pdbin, copies):
  chain_ids =['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
  ensemble = pdb.hierarchy.input(pdbin)
  num_models = len(ensemble.hierarchy.models())
  name = pdbin[:-4]
  selection = pdb.hierarchy.ext.root()
  models = []
  for copy in range(copies):
    selection = pdb.hierarchy.ext.root()
    for model_no in range(num_models):
      selected_model = ensemble.hierarchy.models()[model_no].detached_copy()
      selection.append_model(selected_model)
      chain_length = selection.models()[model_no].chains()[0].residue_groups_size()
      chain1 = selection.models()[model_no].chains()[0]
      chain1.id = chain_ids[copy*2]
      for residue in selection.models()[model_no].chains()[0].residue_groups():
        if int(residue.resseq) > chain_length/2:
          selection.models()[model_no].chains()[0].remove_residue_group(residue)
      chain2 = ensemble.hierarchy.models()[model_no].chains()[0].detached_copy()
      chain2.id = chain_ids[copy*2 + 1]
      for residue in chain2.residue_groups():
        if int(residue.resseq) <= chain_length/2:
          chain2.remove_residue_group(residue)
      selection.models()[model_no].append_chain(chain2)
    model = name + '.chains%s%s' % (chain_ids[copy*2], chain_ids[copy*2 + 1])
    models.append(model)
    with open(model+'.pdb', 'w') as pdbfile:
      print(selection.as_pdb_string(),'END', file=pdbfile)
  return models

def split_strands(pdbin, copies):
  chain_ids =['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
  ensemble = pdb.hierarchy.input(pdbin)
  num_models = len(ensemble.hierarchy.models())
  chain_length = ensemble.hierarchy.models()[0].chains()[0].residue_groups_size()
  name = pdbin[:-4]
  selection = pdb.hierarchy.ext.root()
  models = []
  for copy in range(copies):
    selection = pdb.hierarchy.ext.root()
    for model_no in range(num_models):
      for segment in range(4):
        atoms_to_keep = []
        if segment in [0,2]:
          for resn in range(1, chain_length + 1):
            if resn <= int(chain_length/4) + int(segment*(chain_length/4.)) and resn > int(segment*(chain_length/4.)):
              atoms_to_keep.extend(['%4d N  ' % resn, '%4d CB ' % resn, '%4d CA ' % resn, '%4d C  ' % resn, '%4d O  ' % resn])
            elif resn == int(chain_length/4) + int(segment*(chain_length/4.)) + 1:
              atoms_to_keep.append('%4d N  ' % resn)
          if segment == 0: # first one
            selected_model = ensemble.hierarchy.models()[model_no].detached_copy()
            res_groups = selected_model.chains()[0].residue_groups()
          else:
            next_chain = ensemble.hierarchy.models()[model_no].chains()[0].detached_copy()
            next_chain.id = chain_ids[copy*4 + segment]
            res_groups = next_chain.residue_groups()
          for residue_group in res_groups:
            for atom in residue_group.atom_groups()[0].atoms():
              if atom.parent().parent().resseq + atom.name not in atoms_to_keep:
                residue_group.atom_groups()[0].remove_atom(atom)
          if segment == 0: # first one
            selection.append_model(selected_model)
          else:
            selection.models()[model_no].append_chain(next_chain)
        elif segment in [1,3]:
          for resn in range(1, chain_length + 1):
            if resn >= (int(chain_length/4) + int(segment*(chain_length/4.))) and resn <= int(chain_length/4. + segment*(chain_length/4.)) :
             atoms_to_keep.extend(['%4d N  ' % resn, '%4d CB ' % resn, '%4d CA ' % resn, '%4d C  ' % resn, '%4d O  ' % resn])
            elif resn == int(segment*(chain_length/4.)) + 1:
              atoms_to_keep.extend(['%4d C  ' % resn, '%4d O  ' % resn])
          next_chain = ensemble.hierarchy.models()[model_no].chains()[0].detached_copy()
          next_chain.id = chain_ids[copy*4 + segment]
          for residue_group in next_chain.residue_groups():
            for atom in residue_group.atom_groups()[0].atoms():
              if atom.parent().parent().resseq + atom.name not in atoms_to_keep:
                residue_group.atom_groups()[0].remove_atom(atom)
          selection.models()[model_no].append_chain(next_chain)

    model = name + '.chains%s%s%s%s' % (chain_ids[copy*4], chain_ids[copy*4 + 1],chain_ids[copy*4 + 2],chain_ids[copy*4 + 3])
    models.append(model)
    with open(model+'.pdb', 'w') as pdbfile:
      print(selection.as_pdb_string(),'END', file=pdbfile)
  return models

def split_models(pdbin):
  ensemble = pdb.hierarchy.input(pdbin)
  num_models = len(ensemble.hierarchy.models())
  name = pdbin[:-4]
  models = []
  for model_no in range(num_models):
    selection = pdb.hierarchy.ext.root()
    selected_model = ensemble.hierarchy.models()[model_no].detached_copy()
    selection.append_model(selected_model)
    model = name + '.model%d' % (model_no + 1)
    models.append(model)
    with open(model + '.pdb','w') as pdbfile:
      print(selection.as_pdb_string(),'END', file=pdbfile)
  return models

def run_rnp(root, xml_file, xmlroot, docid, output, logfile, data, tncs, pdbin, models,
            solutions, num_solutions, rescore_all, rescore_models, nproc):
  input = phaser.InputMR_RNP()
  # test OpenMP
  if nproc > 1:
    input.setJOBS(nproc)
  else:
    input.setJOBS(1)
  input.setSPAC_HALL(data.getSpaceGroupHall())
  input.setCELL6(data.getUnitCell())
  input.setREFL_DATA(data.getDATA())
  for model in models:
    input.addENSE_PDB_RMS(model, model + '.pdb', solutions[0].VRMS.values()[0][0])
  if not rescore_all:
    solutions = solutions[:num_solutions]
  for solution in solutions:
    solution_id = root + '.' + str(solution.NUM)
    if rescore_models:
      for model in models:
        input.addSOLU_SET(solution_id + '_' + model.split('.')[1])
        for mr in solution.KNOWN:
          rot = mr.getEuler()
          tra = mr.TRA
          input.addSOLU_6DIM_ENSE(model,[rot[0],rot[1],rot[2]],True,[tra[0],tra[1],tra[2]],mr.BFAC,False,False,False,1.0,1.0)
    else:
      input.addSOLU_SET(solution_id)
      copy = 0
      for mr in solution.KNOWN:
        rot = mr.getEuler()
        tra = mr.TRA
        log.debug('DEBUG copy: %d, len(models): %d, rot: %s, tra %s' %(copy, len(models), rot, tra))
        input.addSOLU_6DIM_ENSE(models[copy],[rot[0],rot[1],rot[2]],True,[tra[0],tra[1],tra[2]],mr.BFAC,False,False,False,1.0,1.0)
        copy += 1
  if tncs:
  # to restore default MR_AUTO behaviour needs this
    input.setTNCS_PATT_PERC(99.9)
  input.setPURG_RNP_ENAB(False)
  input.setOUTP_LEVE('SUMMARY')
  if rescore_models:
    input.addMACM(False,False,False,True,False,False,False,50,'BFGS')
  else:
    input.addMACM(True,True,True,True,False,False,False,50,'BFGS')
    input.setMACM_CHAI(True)
  input.setROOT(root)
  input.setXYZO_ENSE(False)
  input.setTOPF(num_solutions)
  if xml_file is not None and xmlroot is not None:
    input.setOUTP_LEVE('LOGFILE') #i2
  elif docid is not None:
    input.setOUTP_LEVE('LOGFILE') #rvapi
  else:
    input.setOUTP_LEVE('SUMMARY') #command line
  phaser_output = phaser.Output()
  rescore_log = open(logfile, 'w')
  output_object = CallbackObject(xml_file, xmlroot, docid, output)
  phaser_output.setPhenixPackageCallback(output_object)
  phaser_output.set_file_object(rescore_log)
  rnp = phaser.runMR_RNP(input, phaser_output)
  rescore_log.close()
  if rnp.Success():
    log.info('After rescoring models in ensemble there are  %d solution(s)' % rnp.numSolutions())
    log.info('Top solution LLG = %f\n' % rnp.getTopLLG())
    return rnp.getDotSol()
  else:
    return 'No solutions'


def place_fragment(root, xml_file, xmlroot, docid, output, data, pdbin, copies, rms, pdbin_fixed,
                   seqin, ncs_copies, phaser_solvent, sgalternative, sg_list, input_hand, tncs,
                   search_lowres, search_highres, rot_peaks, rot_cluster_off, rot_samp, tra_samp,
                   search_down, tfz_solved, solutions, purge, rescore_strands,
                   rescore_residues, rescore_all, rescore_models, nproc):

  phaser_logfile = root + '_Phaser.log'
  if nproc > 1:
    log.info('    Running Phaser with %d threads, logfile is: %s' % (nproc, os.path.abspath(phaser_logfile)))
  else:
    log.info('    Running Phaser, logfile is: %s' % (os.path.abspath(phaser_logfile)))
  if tncs:
    log.info('\n    tNCS correction not applied')

  phaser_solutions = run_phaser(root=root, xml_file=xml_file, xmlroot=xmlroot, docid=docid, output=output, logfile=phaser_logfile,
                                data=data, pdbin=pdbin, copies=copies, rms=rms, pdbin_fixed=pdbin_fixed,
                                seqin=seqin, ncs_copies=ncs_copies, phaser_solvent=phaser_solvent,
                                sgalternative=sgalternative, sg_list=sg_list, input_hand=input_hand, tncs=tncs,
                                search_lowres=search_lowres, search_highres=search_highres,
                                rot_peaks=rot_peaks, rot_cluster_off=rot_cluster_off, rot_samp=rot_samp,
                                tra_samp=tra_samp, search_down=search_down, tfz_solved=tfz_solved,
                                solutions=solutions, purge=purge, rescore_strands=rescore_strands,
                                rescore_residues=rescore_residues, rescore_models=rescore_models, nproc=nproc)

  if rescore_strands or rescore_residues or rescore_models:
    phaser_logfile = root + '_rescore.log'
    if nproc > 1:
      log.info('    Rescoring models in ensemble with %d threads, logfile is: %s' % (nproc,os.path.abspath(phaser_logfile)))
    else:
      log.info('    Rescoring models in ensemble, logfile is: %s' % (os.path.abspath(phaser_logfile)))
    if tncs:
      log.info('\n    tNCS correction not applied')
    copies_placed = len(phaser_solutions[0].KNOWN)
    log.debug('DEBUG copies %d len(solutions)[0].KNOWN %d' % (copies, len(phaser_solutions[0].KNOWN)))
    if len(phaser_solutions[0].KNOWN) > copies:
      log.info('Phaser placed more copies of the search fragment than requested increasing copies from %d to %d' % (copies, copies_placed))
      copies = copies_placed
    if rescore_strands:
      models = split_chains(pdbin, copies)
    elif rescore_residues:
      models = split_strands(pdbin, copies)
    elif rescore_models:
      models = split_models(pdbin)
    rescored_solutions = run_rnp(root=root, xml_file=xml_file, xmlroot=xmlroot, docid=docid, output=output, logfile=phaser_logfile,
                data=data, tncs=tncs, pdbin=pdbin, models=models, solutions=phaser_solutions,
                num_solutions=solutions, rescore_all=rescore_all, rescore_models=rescore_models, nproc=nproc)

    return rescored_solutions
  else:
    return phaser_solutions
