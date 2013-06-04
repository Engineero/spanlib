#################################################################################
#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2013  Stephane Raynaud
# Contact: stephane dot raynaud at gmail dot com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#################################################################################

import numpy as npy
from util import Logger, broadcast, SpanlibIter, dict_filter
from analyzer import Analyzer, default_missing_value
import pylab as P

class Filler(Logger):
    """Class to fill missing value with a MSSA filtered version of the data
    
    The initialization automatically call the :meth:`fill` method.
    """
    
    def __init__(self, data, run=True, logger=None, loglevel=None, **kwargs):
                
        # Logger
        Logger.__init__(self, logger=logger, loglevel=loglevel, **dict_filter(kwargs, 'log_'))
        self._kwargs = kwargs
        self._kwargs['keep_invalids'] = True
        self._kwargs.setdefault('nvalid', 1)
        self._kwargs.setdefault('quiet', True)
        self._data = data
#        self.filled = None
#        self.filtered = None
        self.nstep = 0
        self.cv = 100
        self._analyzes = []   
        self._kwfill = {}
        self._ana = 'pca'
        
        # Setup analyzer
        span = self.analyze(**kwargs)
        
        # Checks
        for invalid in span.invalids:
            if invalid is not None and invalid.any():
                break
        else:
            self.warning("No gap to fill")
            
        # Keep original data safe
        self._set_field_(span.stacked_data, 'orig')
        
        # Start filling?
        if run: self.fill(**kwargs)
        
    def analyze(self, data=None, npcamax=20, nmssamax=15, **kwargs):
        # Setup analyzer
#            kwargs['keep_invalids'] = False
#            kwargs['nvalid'] = None
        if data is None: data = self._data
        kwargs['keep_invalids'] = True
        kwargs.setdefault('nvalid', 1)
        kwargs.setdefault('quiet', True)
        kwargs['npca'] = npcamax
        kwargs['nmssa'] = nmssamax
        span = self.span = Analyzer(data, sequential=False, logger=self.logger, **kwargs)
        span.update_params()
        return span
        

    def fill(self, nitermax=20, errchmax=-0.01, fillmode='normal', testmode='crossvalid', 
        mssa=True, full=True, cvregen=False, skipat=2, nreana=3, **kwargs):
        """Run the filler with a convergence loop
        
        Results are accessible in the following attributes:
        
        .. attribute:: filtered
        
            Filtered data, result from the convergence loop.
            
        .. attribute:: filled
        
            Data filled with :attr:`filtered`
            
        :Parameters:
        
            - **fillmode**: "zeros" or "normal"
            - **nitermax**: Maximal number of iterations
            - **cvmax**: Convergence criterion (%)
            - **npca**: Number of PCA modes (see :class:`Analyzer`)
            - **nmssa**: Number of MSSA modes (see :class:`Analyzer`)
            - Other parameters are passed to :class:`Analyzer`
            
        :Returns:
        
            - :attr:`filled`
        
        """
    
        
        # Parameters
        self._kwfill.update(nitermax=nitermax, errchmax=errchmax, fillmode=fillmode, 
            testmode=testmode, mssa=mssa, full=full, cvregen=cvregen, 
            skipat=skipat, nreana=nreana, **kwargs)
        kwgencv = dict_filter(kwargs, 'cvfield_')
        span = self.span
        if fillmode==0:
            fillmode = "none"
        fillmode = str(fillmode).lower()
        kwargs['zerofill'] = 0 if fillmode.startswith('n') else 2
        kwargs['prepca'] = True # always PCA
        self.debug('Filling mode: %s (%i)'%(fillmode, kwargs['zerofill']))
        span.update_params(**kwargs)
        
        # which analyzes types?
        analyzes = []
#        if not mssa or span.prepca: 
        analyzes.append('pca')
        if mssa: analyzes.append('mssa')
        self._analyzes = analyzes
        
        import pylab as P
        # Which modes?
        testmode = str(testmode).lower()
        if testmode.startswith('c'):
            kwgencv['regen'] = cvregen
            anamodes = [2, 0]
        else:
            anamodes = [1]
            
        # Loop on analysis modes
        self._nmodes = {}
        self._errors = {}
        for anamode in anamodes:
            self.debug('Analysis mode: '+['NORMAL', 'SELF-VALIDATION', 'CROSS-VALIDATION'][anamode])
            
            # Loop on analysis types (PCA, MSSA)
            self._errors[anamode] = {}
            for self._ana in analyzes:
                self.debug(' Running '+self._ana.upper())
                
                # Update the number of pre-PCA modes for MSSA
                if self._ana=='mssa':
                    span.update_params(prepca = imode+1)
                    
                # Link to appropriate data
                
                # - reference data
                self._link_field_('orig' if self._ana=='pca' else 'pcs', 'ref') 
                rmask = self._refm.mask
                saxis = int(self._ana=='mssa')
                if rmask is npy.ma.nomask or not (rmask.all(axis=saxis)|rmask.any(axis=saxis)).any():
                    self.warning(' No gap to fill -> just analyzing with all modes')
                    self._get_func_()()
                    self._nmodes.setdefault(self._ana, [[self.span.nmssa]])
                    break
                
                # - data to fill
                if anamode==2: # cross-validation
                    self._gen_cvfield_(**kwgencv)
                    self._link_field_('cvfield', 'current')
                else: # normal, self-validation
                    self._link_field_('ref', 'current')
                    
                if self._ana=='mssa':
                    pass
                    
                # Initialize raw data
                self._set_raw_(self._currentm.data)
                
                # Reanalyses loop
                self._errors[anamode][self._ana] = []
                for irf in xrange(nreana):
                    
                    self.debug('  Analysis (%i/%i)'%(irf+1, nreana))                        
                
                    # Run analysis to get EOFs
                    self._get_func_()(force=True)
                    
                    # CV loop on EC estimation (not for PCA with T-EOF)?
                    ecloop = self._ana!= 'pca' or not self.span.useteof
                    niterec = nitermax if ecloop else 1
                    
                    # Number of modes to retain
                    self._nmodes.setdefault(self._ana, [])
                    amodes = self._nmodes[self._ana]
                    if len(amodes)<irf+1:
                        amodes.append(range(getattr(self.span, 'n'+self._ana))) # test all
                    
                    # Loop on the number of modes 
                    last_mode_err = None
                    self._last_pcs_mode = {}
                    self._errors[anamode][self._ana].append([])
                    for imode in amodes[irf]:
                        verb = '   Reconstructing' if anamode==0 else '   Trying'
                        self.debug(verb+' with %i mode%s'%(imode+1, 's'*(imode>0)))
                        
                        # Inits
                        self._recm = default_missing_value
                        if hasattr(self, '_necmiss'): del self._necmiss
                        if hasattr(self, '_necmiss_old'): del self._necmiss_old
                        last_iter_err = None
                        skiplast = False
                        if anamode==1 and imode>0:
                            self._last_pcs_mode[self._ana] = \
                                getattr(self.span, '_%s_raw_pc'%self._ana)
                        self._last_pcs_iter = {}
                        
                        # Convergence loop for expansion coefficients
                        self._errors[anamode][self._ana][-1].append([])
                        for istep in xrange(niterec):
                            if ecloop: self.debug('    EC convergence step: %i'%istep)
                            
                            # Reconstruct
                            if anamode==1 and istep>0:
                                self._last_pcs_iter[self._ana] = \
                                    getattr(self.span, '_%s_raw_pc'%self._ana)
                            self._rec_(imode)
        
                            # Current error
                            err = self._get_error_(anamode)
                            self._errors[anamode][self._ana][-1][-1].append(err)
                            
                            # Check MSSA full filling
                            if self._ana=='mssa' and full:
                                nem, nemch = self._get_necmiss_()
                                if nem and (nemch or nemch is None):
                                    self.debug('    Still %i missing values in %s PCs'%
                                        (nem, self._ana.upper()))
                                    last_iter_err = err
                                    continue
        
                            # Check convergence error for EC
                            if ecloop:
                                self.debug('    Current error: %.1f%%'%err)
                                if istep>0 and last_iter_err is not None:
                                    errch = err - last_iter_err
                                    self.debug('    Error change: %g%%'%errch)
                                    if errch>=errchmax:
                                        if errch>0:
                                            self.debug('  Error change > 0: unstable mode -> step skipped')
                                            err = last_iter_err
                                            if anamode==1:
                                                setattr(self.span, '_%s_raw_pc'%self._ana, 
                                                    self._last_pcs_iter[self._ana])
                                            self.debug('    Recovered PCs from last step')
                                            self._errors[anamode][self._ana][-1][-1] *= -1
        #                                    skiplast = True
                                        else:
                                            self.debug('    Error change > threshold -> stopping convergence'%errch)
                                        break
                            else:
                                self.debug('    Error: %.1f%%'%err)
                                break
                            last_iter_err = err
                        
                        else:
                            self.debug('  Reached max number of iterations for EC convergence loop')
                                        
                        # Check mode truncature error
                        if anamode!=0 and imode>0:
                            if skiplast: 
                                errch = 1
                            else:
                                errch = err-last_mode_err
                                self.debug('   Error change between %i and %i modes: %g'%
                                    (imode, imode+1, errch))
                            if errch>errchmax:
                                imode -= 1
                                self.debug('   Best number of %s modes: %i'%(self._ana.upper(), imode+1))
                                if anamode==1:
                                    setattr(self.span, '_%s_raw_pc'%self._ana, 
                                        self._last_pcs_mode[self._ana])
                                    self.debug('   Recovered PCs from last mode')
                                    if self._errors[anamode][self._ana][-1][-1] >0:
                                        self._errors[anamode][self._ana][-1][-1] *= -1
                                break    
                        last_mode_err = err
                    else:
                        if imode>0:
                            self.debug('   Reached max number of %s modes (%i)'%(self._ana.upper(), imode+1))
                 
                    # Refill
                    if nreana>0:
                        self._set_raw_(self._currentm.filled(self._recm.data))
       
                    # Store optimal number of modes for normal analysis after cross-validation
                    if anamode==2:
                        self._nmodes[self._ana][irf] = [imode]

                    # -> NMODES

                # -> REANA

                # Store PCA pcs for MSSA
                if self._ana=='pca' and 'mssa' in analyzes:
                    self.span.prepca = imode+1
                    self._set_field_(self.span._pca_raw_pc[:, :imode+1].T, 'pcs', mean=True, std=False)
                
            # -> PCA/MSSA
            
        # -> NORM/SELF/CV
               
        
    def get_filtered(self, mssa=None, **kwargs):
        if not self._analyzes:
            kw = self._kwfill.copy()
            kw.update(kwargs, mssa=mssa is not False)
            self.fill(**kw)
        if mssa is False:
            ana = 'pca'
        elif mssa is None:
            ana = self._ana
        elif 'mssa' not in self._analyzes:
            kw = self._kwfill.copy()
            kw.update(kwargs, mssa=True)
            self.fill(mssa=True, **kwargs)
            ana = 'mssa'
        rec_meth = getattr(self.span, ana+'_rec')
        return rec_meth(modes=-self._nmodes[ana][-1][0])
        
    filtered = property(fget=get_filtered, doc='Filtered data')
        
    def get_filled(self, mssa=None, **kwargs):
        filtered = self.get_filtered(mssa=mssa, **kwargs)
        return self.span.fill_invalids(self._data, self.filtered)
        
    filled = property(fget=get_filled, doc='Data filled with filtered data where missing')


    def _set_field_(self, field, name, mean=True, std=False):
        """Put a field in self._<name>*"""
        fieldm = npy.ma.masked_values(field, default_missing_value, copy=False)
        setattr(self, '_%sm'%name, fieldm)
        if std:
            setattr(self, '_%s_std'%name, fieldm.std())
        if mean:
            setattr(self, '_%s_mean'%name, fieldm.mean(axis=1).reshape(-1, 1))

    def _link_field_(self, name, to='current'):
        """
        1.) Put a stored field and its attributes in self._current
        2.) Put data it in the current field to be analyzed (stacked data or PCs)
        """
        # Put it in self._current*
        if self._ana=='mssa':
            pass
        for att in 'm', '_mask', '_std', '_mean':
            if hasattr(self, '_%s'%name+att):
                setattr(self, '_'+to+att, getattr(self, '_'+name+att))

    def _set_raw_(self, data):
        if self._ana=='mssa':
            data = data.T
            raw_name = '_pca_raw_pc'
        else:
            raw_name = 'stacked_data'
        if npy.ma.isMA(data): data = data.data
        setattr(self.span, raw_name, data)
        

#        # Put it in analyzed field
#        if self._ana=="mssa":
#            self.span._pca_raw_pc = self._current
#        else:
#            self.span.stacked_data = self._current
    
        
    def _get_error_(self, anamode):
        """Get current reconstruction error"""
        diffm = self._refm-self._recm
        fieldm = self._refm
        if anamode==2:
            diffm = npy.ma.masked_where(self._cvfield_kept, diffm, copy=False)
            fieldm = npy.ma.masked_where(self._cvfield_kept, fieldm, copy=False)
        return 100*diffm.compressed().std()/fieldm.compressed().std()
        
     
    def _rec_(self, imode):
        """Get PCA or MSSA recontruction
        
        For PCA, output is same as raw input.
        For MSSA, output is same as PCA PCs input
        """
        
        # Expansion coefficients
        self._last_pcs = {} #getattr(self.span, '_%s_raw_pc'%self._ana)[0]
        if self._ana=='mssa' or not self.span.useteof: 
            
            # Fill masked point with reconstructed data
            rec = self._recm.data if hasattr(self._recm, 'data') else self._recm
            data =  self._currentm.filled(rec)
            
            if self._ana=='mssa': # MSSA rec
            
                kw = dict(xdata=data, xraw=True) if not self.span.prepca else {}
                self.span.mssa_ec(raw=True, replace=True, demean=False, **kw)
                
            else: # PCA rec
            
                self.span.pca_ec(raw=True, xraw=True, xdata=data, replace=True, demean=False)
                
            del data
    
        # Reconstruction (masked)
        recfunc = self._get_func_('rec')
        self._recm = recfunc(modes=-imode, raw=1, rescale=False, unmap=False)
        if hasattr(self, '_current_mean'):
            self._recm += self._current_mean
    
    def _get_func_(self, suf=None):
        """Get a PCA or MSSA related generic function (method)"""
        ana = self._ana
        if suf is not None:
            ana += '_'+suf
        return getattr(self.span, ana)
        
    def _get_necmiss_(self):
        """Get the number of missing values of current expansion coefficents"""
        if not hasattr(self,  '_ana'): return None, None
        if hasattr(self, '_necmiss'): self._necmiss_old = self._necmiss
        pc0 = npy.ma.masked_values(getattr(self.span, '_%s_raw_pc'%self._ana)[:, 0], 
            default_missing_value)
        self._necmiss =  pc0.size-npy.ma.count(pc0)
        self._necmiss =  pc0.size-npy.ma.count(pc0)
        del pc0
        if not hasattr(self, '_necmiss_old'):  
            necmissch = None
        else:
            necmissch = self._necmiss-self._necmiss_old
        return self._necmiss, necmissch
        
    def _gen_cvfield_(self, mask=None, level=5, sscale=1, tscale=1, minrelscale=0.2, regen=False):
        """Generate false missing values
        
        Field used for analysis is stored using `set_field('cvfield')`.
        Retained points are saved as logical array in :attr:`_cvfield_kept`.
        
        :Params:
        
            - **mask**, optional: Mask set to True where values are not retained for
              analysis. These rejected values are used to cross-validation.
              If not set, it is randomly generated.
        """
        if not regen and hasattr(self, '_cvfield') and self._cvfield_ana==self._ana: return

        # Generate the mask
        if mask is None:
            
            # Smooth scale
            ns, nt = self._refm.shape
            transp = self._ana=='mssa'
            if transp: ns, nt = nt, ns
            if sscale<0:
                sscale = -int(ns*sscale/100.)
            if sscale > minrelscale*ns:
                sscale = 1
            else:
                sscale = sscale/2*2+1
            if tscale<0:
                tscale = -int(nt*tscale/100.)
            if tscale > minrelscale*nt:
                tscale = 1
            else:
                tscale = tscale/2*2+1
            
            # Random values
            tmp = npy.random.random(self._refm.shape).astype('f')
            
            # Apply smoothing
            if tscale>1 or sscale>1:
                ww = npy.ones((ns, nt), 'f')
                from numpy import convolve
                kernel = npy.ones(sscale, 'f')
                if sscale>1:
                    for i in xrange(nt):
                        thist = (slice(None), i)[::1-2*reverse]
                        tmp[thist] = convolve(tmp[thist], kernel, 'same') / convolve(ww[thist], kernel, 'same')
                kernel = npy.ones(tscale, 'f')
                if tscale>1:
                    for i in xrange(ns):
                        thiss = (i, slice(None))[::1-2*reverse]
                        tmp[thiss] = convolve(tmp[thiss], kernel, 'same') / convolve(ww[thiss], kernel, 'same')
                del ww
            
            # Get max value for noisy field
            sorted = npy.sort(tmp.ravel())
            nlevel = int(npy.round(tmp.size*level/100.))
            maxval = sorted[max(0, nlevel-1)]
            del sorted
            
            # Mask
            mask = tmp<maxval
            del tmp
        
        # Apply mask
        cvfield = npy.where(mask, default_missing_value, self._refm)
        self._set_field_(cvfield, 'cvfield', mean=True, std=False)
        self._cvfield_kept = ~mask
        self._cvfield_ana = self._ana
        
        
