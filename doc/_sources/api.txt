Full API Documentation
===========================

oPDF module
-------------------------
.. automodule:: oPDF

Variables and Constants
*************************
.. py:data:: Globals
   
   Collection of global variables of the module, of class :py:class:`globals_t`. It controls numerical precision, internal units, and cosmology. 
   
.. py:data:: VirTypes
   
   Collection of virial definitions. It contains 
    - ``VirTH``: the spherical collapse prediction (i.e, Bryan & Norman 98 fitting).
    - ``VirB200``: the 200 times mean density deifinition.
    - ``VirC200``: the 200 times critical density definition.
      
.. py:data:: HaloTypes
   
   Collection of halo parametrizations. It contains
    -  ``NFWMC``: NFW halo parametrized by :math:`(M,c)`
    -  ``NFWRhosRs``: NFW, :math:`(\rho_s,r_s)`
    -  ``NFWPotsRs``: NFW, (:math:`\psi_s, r_s`\ ), with :math:`\psi_s=4\pi G\rho_s r_s^2`\ .
    -  ``CorePotsRs``: Cored Generalized NFW Potential (inner density slope=0), parametrized by (:math:`\psi_s,r_s`\ )
    -  ``CoreRhosRs``: Cored GNFW, :math:`(\rho_s,r_s)`
    -  ``TMPMC``: Template profile, :math:`(M,c)` parametrization
    -  ``TMPPotScaleRScale``: Template, :math:`\psi_s/\psi_{s0}, r_s/r_{s0}`
      
.. py:data:: Estimators
   
   Collection of dynamical estimators. It contains
    - ``RBinLike``: binned radial likelihood. 
	Use ``RBinLike.nbin`` (``integer``) and ``RBinLike.logscale`` (``True`` or ``False``) to control the number and scale of bins.
    - ``AD``: Anderson-Darling distance.
    - ``MeanPhaseRaw``: Normalized mean phase deviation :math:`\bar{\Theta}=(\bar{\theta}-0.5)/\sigma_{\theta}`\ , to be compared to a standard normal variable.
    - ``MeanPhase``: :math:`\bar{\Theta}^2`\ , to be compared to a chi-square variable.
   
Classes
*************************
Global Parameters
~~~~~~~~~~~~~~~~~~
.. autoclass:: globals_t
   :members:

Halo
~~~~~~~~~~~~~~~~~~
.. autoclass:: Halo
   :members:

Tracer
~~~~~~~~~~~~~~~~~~
.. autoclass:: Tracer
   :members:
   
Utility functions
------------------------
.. automodule:: myutils
   :members: