cimport TABIPBstruct
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free

import logging

_log = logging.getLogger()

cdef class TABIPB_Solver:

  cdef double density_ #density for mesh
  cdef double probe_radius_ #probe radius
  cdef double epsp_
  cdef double epsw_
  cdef double bulk_strength_
  cdef int order_
  cdef int maxparnode_
  cdef double theta_
  cdef double temp_
  cdef int mesh_flag_

  cdef int _pid

  def __cinit__(self):
    self.density_ = 0
    self.probe_radius_ = 0
    self.epsp_ = 0
    self.epsw_ = 0
    self.bulk_strength_ = 0
    self.order_ = 0
    self.maxparnode_ = 0
    self.theta_ = 0
    self.temp_ = 0
    self.mesh_flag_ = 0

  def __init__(self, density, probe_radius, epsp, epsw, bulk_strength,
      order, maxparnode, theta, temp, mesh_flag):
    self.density_ = density
    self.probe_radius_ = probe_radius
    self.epsp_ = epsp
    self.epsw_ = epsw
    self.bulk_strength_ = bulk_strength
    self.order_ = order
    self.maxparnode_ = maxparnode
    self.theta_ = theta
    self.temp_ = temp
    self.mesh_flag_ = mesh_flag

    import os
    self._pid = os.getpid()

  cdef _run_tabipb(self, molecules):
    cdef int natm, i
    cdef double[:,:] xyzr
    cdef TABIPBstruct.TABIPBparm tabipbparm
    cdef TABIPBstruct.TABIPBvars tabipbvars

    natm = len(molecules['atoms'])

    tabipbparm.density = str(self.density_);
    tabipbparm.probe_radius = str(self.probe_radius_);
    tabipbparm.temp = self.temp_;
    tabipbparm.epsp = self.epse_;
    tabipbparm.epsw = self.epsw_;
    tabipbparm.bulk_strength = self.bulk_strength_
    tabipbparm.order = self.order_;
    tabipbparm.maxparnode = self.maxparnode_;
    tabipbparm.theta = self.theta_;
    tabipbparm.mesh_flag = self.mesh_flag_;
    tabipbparm.number_of_lines = natm;

    tabipbvars.chrpos = <double*> malloc(3 * natm * sizeof(double));
    tabipbvars.atmchr = <double*> malloc(natm * sizeof(double));
    tabipbvars.atmrad = <double*> malloc(natm * sizeof(double));

    for i, atom in enumerate(molecules['atoms']):
      tabipbvars.chrpos[i*3] = atom['pos'][0]
      tabipbvars.chrpos[i*3+1] = atom['pos'][1]
      tabipbvars.chrpos[i*3+2] = atom['pos'][2]
      tabipbvars.atmchr[i] = atom['charge']
      tabipbvars.atmrad[i] = atom['radius']

    TABIPBstruct.sphinx2tabipb(&tabipbparm, &tabipbvars)

    free(tabipbvars.chrpos)
    free(tabipbvars.atmchr)
    free(tabipbvars.atmrad)

    return {'Solvation energy': tabipbvars.soleng,
          'Max Surface potential': tabipbvars.max_vert_ptl,
          'Min Surface potential': tabipbvars.min_vert_ptl}

  def run_solve(self, molecules):
    '''Method called from plugin.py, seems to be wrapping
       to internal stuff'''

    results = []

    _log.info("({}): TABIPB flow solver started.".format(self._pid))

    results.append(self._run_tabipb(molecules))

    _log.info("({}): TABIPB flow solver done.".format(self._pid))
    return results
