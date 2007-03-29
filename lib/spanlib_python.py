#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006  Charles Doutiraux, Stephane Raynaud
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

import spanlib_fort,MV,Numeric as N,genutil.statistics


def stackData(*data):
	""" Takes several data files, of same time and stacks them up together

	Description:::
	  This fonction concatenates several dataset that have the
	  same time axis. It is useful for analysing for example
	  several variables at the same time.
	  It takes into account weights, masks and axes.
	:::

	Usage:::
	dout, weights, mask, axes = stackData(data1[, data2...])

	  *data   :: One or more data objects to stack.
				 They must all have the same time length.
	:::

	Output:::
	  dout	:: Stacked data
	  weights :: Associated stacked weights
	  masks   :: Associated stacked masks
	  axes	:: Associated stacked axes
	:::
	"""
	len_time=None
	axes=[]
	dout=None # data output
	for d in data:
		d = MV.array(d)
		t=d.getTime()
		if t is None:
			t = d.getAxis(0)
			t.designateTime()
			d.setAxis(0,t)
			#raise 'Error, all data must have a time dimension'
		if len_time is None:
			len_time=len(t)
		elif len_time!=len(t):
			raise 'Error all datasets must have the same time length!!!!'

		if d.getAxis(0)!=t:
			d=d(order='t...')

		axes.append(d.getAxisList())
		tdata,w,m=pack(d)
		if dout is None: # Create
			dout=tdata
			weights=w
			masks=[m]
		else: # Append
			dout=N.concatenate((dout,tdata))
			weights=N.concatenate((weights,w))
			masks.append(m)

	return N.transpose(dout),weights,masks,axes

def unStackData(din,weights,masks,axes):
	""" Unstack data in the form returned from stackData

	Description:::
	  This function is the reverse operation of stakData.
	  It splits stacked datasets into a list.
	:::

	Usage:::
	dout = unStackData(din,weights,mask,axes)

	  din	 :: Stacked data (see stackData function)
	  weights :: Associated stacked weights
	  masks   :: Associated stacked masks
	  axes	:: Associated stacked axes
	:::

	Output:::
	  dout	:: List of unstacked data
	:::
	"""
	nvar=len(axes)

	if nvar!=len(masks):
		raise 'Error masks and input data length not compatible'

	totsize=0
	for m in masks:
		totsize+=int(MV.sum(MV.ravel(m)))
	if totsize!=din.shape[1]:
		raise 'Error data and masks are not compatible in length!!!! (%s) and (%s)' \
		      % (totsize,din.shape[1])

	istart=0
	dout=[]
	for i in range(nvar):
		m=masks[i]
		mlen=int(MV.sum(MV.ravel(m)))
		iend=istart+mlen
		data=N.transpose(din[:,istart:iend])
		w=weights[istart:iend]
		up=spanlib_fort.chan_unpack(m,data,1.e20)
		unpacked = MV.transpose(MV.array(up))
		unpacked = MV.masked_where(N.equal(N.resize(N.transpose(m),unpacked.shape),0),unpacked)
		unpacked.setAxisList(axes[i])
		istart+=mlen
		dout.append(unpacked)
	return dout


def pack(data,weights=None):
	""" Pack a dataset and its weights according to its mask

	Description:::
	  This function packs a dataset in 2D space by removing
	  all masked points and returning a space-time array.
	  It performs this operation also on the weights.
	  It is used for removing unnecessary points and
	  simplifying the input format for analysis functions.
	:::

	Usage:::
	packed_data, packed_weights, mask = pack(data,weights)

	  data	  :: Flatten in space an [x,y,t] array by removing
				 its masked point
	  weights :: Weights to be flatten also
	:::

	Output:::
	  packed_data	 :: Space-time packed array
	  packed_weights :: Packed weights that were guessed or used
	  mask		     :: Mask that were guessed or used
	:::
	"""

	# Total number of channels
	data = MV.array(data)
	nt = data.shape[0]
	nstot = N.multiply.reduce(data.shape[1:])

	# Is it already packed? Check the mask...
	mask = data.mask()
	if not N.sum(mask.flat):
		if weights is None:
			weights = N.ones(nstot,typecode='f')
		else:
			weights = N.array(weights)
		packed_data = N.transpose(MV.reshape(data,(nt,nstot)))
		packed_weights = N.ravel(weights)
		mask = None
		return packed_data,packed_weights,mask

	sh=list(data.shape)

	# Weights ?
	if weights is None:
		try:
			import cdutil
			tmp=data
			while tmp.rank()>2:
				tmp=tmp[0]
			weights=cdutil.area_weights(tmp).raw_data()
			del(tmp)
		except Exception,err:
			weights=MV.ones(nstot,typecode='f')
	else:
		weights = MV.array(weights)

	# Mask
	# - First from data
	mask = data.mask()
	if mask is None:
		mask = MV.zeros(sh,typecode='f')
	else:
		mask = mask.astype('i')
	mask = MV.sum(mask,axis=0)
	# - Now add the ones from the weights
	mask = mask.filled()+MV.equal(weights,0.).filled(1)
	# - >=1 means masked, Fortran "mask": 1 means data ==> 1-mask
	mask = 1.-MV.greater_equal(mask,1).filled()
	mask = N.ravel(mask)

	# Number of points in spatial dimension
	ns = int(MV.sum(N.ravel(mask)))

	# Ok now calls fortran, but need to transpose first
	# and reduce dimensions
	data_to_pack = N.transpose(N.reshape(data.filled(1.e20),
													 (nt,nstot)))

	# Dummy 1D for time for tmask
	# - Pack data
	packed_data = spanlib_fort.chan_pack(data_to_pack,mask,ns)
	del data_to_pack
	# - Pack weights
	#weights=MV.reshape(weights.filled(0),(nstot,1))
	#tweights=N.transpose(weights.filled(0)) # FIXME
	#tweights=N.ones(tweights.shape,'f')
	tweights=N.ones((nstot,1),'f')
	packed_weights = spanlib_fort.chan_pack(tweights,mask,ns)[:,0].astype('f')
	return packed_data,packed_weights,mask


def getPairs(pcs,minCorr=0.95,maxDistance=3):
	""" Get pairs in MSSA results

	Description:::
	  This function detects pairs of mode in principal
	  components (typical from MSSA). Modes of a pair
	  have their PCs (and EOFs) in phase quadrature in time.
	  The algorithm uses correlations between a pc and the 
	  derivative of the other pc.
	:::

	Usage:::
	pairs = getPairs(pcs)

	  pcs         :: Principal components [modes,time]
	  minCorr     :: Minimal correlation [default:0.95]
	  maxDistance :: Maximal difference of indices between a pair of modes [default: 3]
	:::

	Output:::
	  pairs :: List of two-element lists if indices
	:::
	"""

	pairs = []
	flat = []
	minCorr = N.clip(0.,1.)
	maxDistance = N.clip(1,len(pcs)-1)
	for i in xrange(len(pcs)):
		if i in flat:
			continue
		ii = [i]
		for j in xrange(1,maxDistance+1):
			ii.append(i+j)
			corr = 0.
			for k in 0,1:
				pc = pcs[ii[k],1:-1]
				dpc = pcs[ii[1-k],2:] - pcs[ii[1-k],:-2]
				corr += float(genutil.statistics.correlation(pc,dpc))
			if corr >= 2*minCorr:
				break
		else:
			continue
		pairs.append([i,j])
		flat.extend([i,j])

	return pairs


def computePhases(data,nphases=8,offset=.5,firstphase=0):
	""" Phase composites for oscillatory fields

	Description:::
	  This computes temporal phase composites of a spatio-temporal
	  dataset. The dataset is expected to be oscillatory in time.
	  It corresponds to a reoganisation of the time axis to
	  to represents the dataset over its cycle in a arbitrary
	  number of phases. It is useful, for example, to have a
	  synthetic view of an reconstructed MSSA oscillation.
	:::

	Usage:::
	phases = computePhases(data,nphases,offset,firstphase)

	  data	   :: Space-time data oscillatory in time data.shape is rank 2 and dim 0 is space
	  nphases	:: Number of phases (divisions of the cycle)
	  offset	 :: Normalised offset to keep higher values only [default:
	  firstphase :: Position of the first phase in the 360 degree cycle
	:::

	Output:::
	  phases :: Space-phase array
	:::
	"""



class SpAn(object):
	def __init__(self,data,weights=None,npca=None,window=None,nmssa=None,nsvd=None,relative=False):
		""" Prepare the Spectral Analysis Object

		Description:::
		  This function creates an object for future analyses.
		  It optionally initializes some parameters.
		:::

		Usage:::
		analysis_object = SpAn(data,weights=None,npca=None,window=None,nmssa=None)

		  data	:: List of data on which to run the PC Analysis
					 Last dimensions must represent the spatial dimensions.
					 Analysis will be run on the first dimension.
		  weights :: If you which to apply weights on some points,
					 set weights to "0" where you wish to mask.
					 The input data mask will be applied,
					 using the union of all none spacial dimension mask.
					 If the data are on a regular grid, area weights
					 will be generated, if the cdutil (CDAT) module is available.
					 [default: 1. everywhere]
		  npca	:: Number of principal components to return [default: 10]
		  nmssa   :: Number of MSSA modes retained [default: 4]
		  nsvd	:: Number of SVD modes retained [default: 10]
		  window  :: MSSA window parameter [default: time_length/3.]
		:::

		Output:::
		  analysis_object :: SpAn object created for further analysis
		:::
		"""

		## Sets all values to None
		self.clean()

		## Before all makes sure data is list of data
		if not isinstance(data,(list,tuple)):
			data=[data,]
		if weights is None:
			weights=[None,] * len(data)
		elif not isinstance(weights,(list,tuple)):
			weights = [weights,]

		## First pack our data, prepare the weights and mask for PCA
		self.pdata=[]
		self.weights=[]
		self.mask=[]
		self.attributes=[]
		self.shapes=[]
		for i in range(len(data)):
			d=MV.array(data[i])
			w=weights[i]
			tmp = pack(d,w)
			tmpdata,tmpweights,tmpmask = tmp
			self.pdata.append(tmpdata)
			self.weights.append(tmpweights)
			self.mask.append(tmpmask)
			self.attributes.append(d.attributes)
			self.shapes.append(d.shape)

		## Store axes for later
		self.axes=[]
		self.varname=[]
		self.grids=[]
		for d in data:
			self.axes.append(d.getAxisList())
			self.varname.append(d.id)
			self.grids.append(d.getGrid())

		## Figures out length of time dimension
		self.nt = data[0].shape[0]

		for d in data:
			if d.shape[0] != self.nt:
				raise Exception, 'Error your dataset are not all consistent in time length'

		if npca is None:
			self.npca=10
		else:
			self.npca=npca

		self.ns=[]
		for p in self.pdata:
			self.ns.append(p.shape[0])

		if nmssa is None:
			self.nmssa = 4
		else:
			self.nmssa = nmssa

		if nsvd is None:
			self.nsvd = 10
		else:
			self.nsvd = nsvd

		if window is None:
			self.window = int(self.nt/3.)
		else:
			self.window = window
		print "Object built!"

	def pca(self,npca=None,get_ev_sum=False,relative=False):
		""" Principal Components Analysis (PCA)

		Descriptions:::
		  This function performs a PCA on the analysis objects
		  and returns EOF, PC and eigen values.
		  EOF are automatically unpacked.
		:::

		Usage:::
		eof, pc, ev = pca(npca=None,weights=None,relative=False)

		OR

		eof, pc, ev, ev_sum = pca(npca=None,weights=None,get_ev_sum=True,relative=False)

		  npca	:: Number of principal components to return [default: 10]
		  get_ev_sum  :: Also return sum of all eigen values [default: False]
		  relative :: Egein values are normalized to their sum (%) [default: False]
		:::

		Output:::
		  eof	:: List EOF array (one per data input when created SpAn object)
		  pc	 :: List of Principal Components array
		  ev	 :: List of Eigein Values array
		  ev_sum :: Sum of all eigen values (even thoses not returned).
					Returned ONLY if get_ev_sum is True.
					It can also be retreived with <SpAn_object>.ev_sum.
		:::
		"""

		if npca is not None:
			self.npca = npca


		## Calls Fortran pca
		self.eof=[]
		self.pc=[]
		self.ev=[]
		eof=[]
		pc=[]
		ev=[]
		self.ev_sum=[]

		for i in range(len(self.pdata)):

			# Compute PCA
			pdat=self.pdata[i]
			w=self.weights[i]
			nteof,ntpc,ntev,ntev_sum = spanlib_fort.pca(pdat,self.npca,w,1)

			# Get axes and dimensions
			Axes=list(self.axes[i][1:])
			sh = [self.npca,]
			sh.extend(self.shapes[i][1:])

			# Recover physical dimensions for EOFs
			if self.mask[i] is not None:
				# Unpack eofs before...
				teof = MV.transpose(MV.array(spanlib_fort.chan_unpack(self.mask[i],nteof,1.e20)))
				teof = MV.reshape(teof,sh)
				tmask = N.resize(N.transpose(self.mask[i]),teof.shape)
				teof = MV.masked_where(N.equal(tmask,0),teof)
				del tmask
			else:
				# Just recover dimensions
				teof = MV.reshape(MV.transpose(nteof),sh)
			ax=teof.getAxis(0)
			ax.id='pca_mode'
			ax.standard_name='PCA Modes in decreasing order'
			Axes.insert(0,ax)
			teof.setAxisList(Axes)
			teof.setGrid(self.grids[i])
			teof.id=self.varname[i]+'_pca_eof'
			teof.name = teof.id
			teof.standard_name='PCA Empirical Orthogonal Functions'
			teof.long_name = teof.standard_name

			# PCs
			tpc=MV.transpose(MV.array(ntpc,axes=[self.axes[i][0],ax]))
			tpc.id=self.varname[i]+'_pca_pc'
			tpc.standard_name='PCA Principal Components'

			# EVs
			tev=MV.array(ntev,id=self.varname[i]+'_pca_ev',axes=[ax])
			tev.standard_name='PCA Eigen Values'
			if relative:
				tev[:] = tev[:] * 100. / ntev_sum
				tev.units = '%'

			# Standard attributes
			for var in tpc,tev,ax:
				var.name = var.id
				var.long_name = var.standard_name


			self.pc.append(ntpc)
			self.eof.append(nteof)
			eof.append(teof)
			pc.append(tpc)
			ev.append(tev)
			self.ev_sum.append(ntev_sum)

		if len(eof)==1:
			ret =  [eof[0],pc[0],ev[0]]
			self.ev_sum = self.ev_sum[0]
		else:
			ret =  [eof,pc,ev]

		if get_ev_sum:
			ret.append(self.ev_sum)

		return ret


	def mssa(self,nmssa=None,pca=None,window=None,get_ev_sum=False,relative=False):
		""" MultiChannel Singular Spectrum Analysis (MSSA)

		Descriptions:::
		  This function performs a MSSA on the analysis objects
		  and returns EOF, PC and eigen values.
		  Unless pca parameter is set to false, a pre
		  PCA is performed to reduced the number of d-o-f
		  if already done and if the number of channels is
		  greater than 30.
		:::

		Usage:::
		eof, pc, ev = mssa(nmssa,pca,relative=False)

		OR

		eof, pc, ev, ev_sum = mssa(nmssa,pca,get_ev_sum=True,relative=False)

		  nmssa  :: Number of MSSA modes retained
		  window :: MSSA window parameter
		  pca	:: If True, performs a preliminary PCA
		  get_ev_sum  :: Also return sum of all eigen values (default: False)
		  relative :: Egein values are normalized to their sum (%) [default: False]

		Output:::
		  eof :: EOF array
		  pc  :: Principal Components array
		  ev  :: Eigen Values  array
		  ev_sum :: Sum of all eigen values (even thoses not returned).
					Returned ONLY if get_ev_sum is True.
					It can also be retreived with <SpAn_object>.stev_sum.
		:::
		"""

		## Check for default values for mssa and pca if not passed by user
		if pca is None:
			if self.pc ==[] and max(0,self.ns) > 30: # Pre-PCA needed
				print '[mssa] The number of valid points is greater than',30,' so we perform a pre-PCA'
				pca = True
			elif self.pc is not None:
				pca = True
			else:
				pca = False

		if pca is True: # From PCA to MSSA
			nspace = [self.npca,]*len(self.pdata)
			if self.pc ==[]: # Still no PCA done
				self.pca()
		else:
			nspace = self.ns

		if nmssa is not None:
			self.nmssa = nmssa

		if window is not None:
			self.window = window

		self.steof = []
		self.stpc  = []

		eof=[]
		ev=[]
		pc=[]
		self.stev_sum=[]
		self.pairs = []

		for i in range(len(self.pdata)):
			if pca is True: # Pre-PCA case
				ntsteof, ntstpc, ntstev, ntev_sum = \
				  spanlib_fort.mssa(N.transpose(self.pc[i]), self.window, self.nmssa)
			else: # Direct MSSA case
				ntsteof, ntstpc, ntstev, ntev_sum = \
				  spanlib_fort.mssa(self.pdata[i], self.window, self.nmssa)


			teof = MV.transpose(MV.reshape(ntsteof,(self.window,nspace[i],self.nmssa)))
			teof.id=self.varname[i]+'_mssa_steof'
			teof.standard_name='Empirical Orthogonal Functions'

			ax0=teof.getAxis(0)
			ax0.id='mssa_chan'
			ax0.standard_name='MSSA Channels'

			ax1=teof.getAxis(1)
			ax1.id='mssa_mode'
			ax1.standard_name='MSSA Modes in decreasing order'

			ax2=teof.getAxis(2)
			ax2.id='mssa_window'
			ax2.standard_name='MSSA Window Axis'

			tpc=MV.transpose(MV.array(ntstpc))
			tpc.id=self.varname[i]+'_stpc'
			tpc.standard_name='Principal Components'
			tpc.setAxis(0,ax0)

			ax3 = tpc.getAxis(1)
			ax3.id='time'
			ax3.standard_name = 'Time'

			tev = MV.array(ntstev,id=self.varname[i]+'_mssa_ev',axes=[ax0])
			tev.standard_name = 'MSSA Eigen Values'
			if relative:
				tev[:] = tev[:] * 100. / ntev_sum
				tev.units = '%'

			for var in teof,tpc,tev,ax0,ax1,ax2,ax3:
				var.name = var.id
				var.long_name = var.standard_name

			self.stpc.append(ntstpc)
			self.steof.append(ntsteof)
			self.pairs.append(getPairs(tpc))
			eof.append(teof)
			pc.append(tpc)
			ev.append(tev)
			self.stev_sum.append(ntev_sum)

		if len(eof)==1:
			self.pairs = self.pairs[0]
			ret = [eof[0],pc[0],ev[0]]
			self.stev_sum = self.stev_sum[0]
		else:
			ret = [eof,pc,ev]

		if get_ev_sum:
			ret.append(self.ev_sum)

		return ret


	def svd(self,nsvd=None,pca=None):
		""" Singular Value Decomposition (SVD)

		Descriptions:::
		  This function performs a SVD
		  and returns EOF, PC and eigen values.
		  Unless pca parameter is set to false, a pre
		  PCA is performed to reduced the number of d-o-f
		  if already done and if the number of channels is
		  greater than 30.
		:::

		Usage:::
		eof, pc, ev = svd(nsvd,pca)

		  nsvd  :: Number of SVD modes retained
		  window :: MSSA window parameter
		  pca	:: If True, performs a preliminary PCA

		Output:::
		  eof :: EOF array
		  pc  :: Principal Components array
		  ev  :: Eigen Values  array

		:::
		"""

		## Check we have at least 2 variables!!
		## At the moment we will not use any more variable
		if len(self.pdata)<2:
			raise Exception,'Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa'

		## Check for default values for mssa and pca if not passed by user
		if pca is None:
			if self.pc ==[] and max(0,self.ns) > 30: # Pre-PCA needed
				print '[svd] The number of valid points is greater than',30,' so we perform a pre-PCA'
				pca = True
			elif self.pc is not None:
				pca = True
			else:
				pca = False

		if pca is True: # From PCA to MSSA
			nspace = [self.npca,]*len(self.pdata)
			if self.pc ==[]: # Still no PCA done
				self.pca()
		else:
			nspace = self.ns

		if nsvd is not None:
			self.nsvd = nsvd


		if pca is True: # Pre-PCA case
			lneof, rneof, lnpc, rnpc, nev = \
			  spanlib_fort.svd(N.transpose(self.pc[0]), 
			    N.transpose(self.pc[1]), self.nsvd)
		else: # Direct SVD case
			lneof, rneof, lnpc, rnpc, nev = \
			  spanlib_fort.svd(self.pdata[0], self.pdata[1], self.nsvd)

		self.svdeof = [lneof,rneof]
		self.svdpc = [lnpc,rnpc]

		eof=[]
		pc=[]

		for i in range(2):
			teof = MV.transpose(self.svdeof[i])
			teof.id=self.varname[i]+'_svd_eof'
			teof.standard_name='SVD Empirical Orthogonal Functions'

			ax0=teof.getAxis(0)
			ax0.id='svd_chan'
			ax0.standard_name='Channels'

			ax1=teof.getAxis(1)
			ax1.id=self.varname[i]+'_svd_mode'
			ax1.standard_name='SVD Modes in decreasing order'


			tpc=MV.transpose(MV.array(self.svdpc[i]))
			tpc.id=self.varname[i]+'_svd_pc'
			tpc.standard_name='SVD Principal Components'
			tpc.setAxis(0,ax0)

			ax3 = tpc.getAxis(1)
			ax3.id='time'

			tev=MV.array(ntstev,id=elf.varname[i]+'_svd_ev',axes=[ax0])
			tev.standard_name='SVD Eigen Values'

			eof.append(teof)
			pc.append(tpc)

		return eof[0],pc[0],eof[1],pc[1],ev

	def _reconstruct(self,function,imode,reof,rpc,rns,rnt,rwin=False):

		# Rearrange modes ([1,3,4,5,9] -> [1,1],[3,5],[9,9])
		imode.sort()
		imodes = []
		while im1 < len(imode):
			imode1 = imode2 = last_imode = imode[im1]
			im2 = im1
			imt = im2+1
			while imt < len(imode):
				if (last_imode-imode[imt]) < 2:
					last_imode = imode[imt]
					continue
				else:
					break
			if last_imode != imode2:
				imode2 = last_imode
			imodes.append([imode1,imode2])
		
		ffrec=[]
		for i in range(len(self.pdata)):
			args = [reof[i],rpc[i],rns[i],rnt]
			if rwin: args.append(rwin)
			for j,ims in enumerate(imodes):
				args.extend(ims)
				this_ffrec = function(*args)
				if not j:
					ffrec = this_ffrec
				else:
					ffrec += this_ffrec
				del this_ffrec
		return ffrec


	def reconstruct(self,imode=None,mssa=None,pca=None,phases=False,nphases=8,offset=.5,firstphase=0,svd=None):
		""" Reconstruct results from mssa or pca

		Description:::
		  This function performs recontructions to retreive the
		  the contribution of a selection of modes to the original field.
		  By default, it recontructs from available PCA and MSSA
		  results. Recontruction of MSSA modes also calls recontruction
		  from of pre-PCA to get back to the original space.
		  This function can optionally performs phase composites
		  (useful for pairs of MSSA modes = oscillations) on MSSA
		  recontructions.
		:::

		Usage:::
		ffrec = reconstruct(imode,mssa,pca)

		  imode :: Selection of modes [default: None]. If:
		    - None: all modes
		    - positive integer: only this mode
		    - negative integer: all modes until abs(imode)
		    - list of modes: use it directly
		  end   :: Last mode
		  mssa  :: Reconstruct MSSA if True
		  pca   :: Reconstruct PCA if True
		  phases :: Operate phase reconstruction True/False (default is False)
		:::

		Output:::
		  ffrec :: Reconstructed field
		:::
		"""

		if type(imode) is type(1):
			if imode < 0:
				imode = range(-imode)
			else:
				imode = [imode-1,]
		imode = (N.array(imode)+1).tolist()

		ntimes=self.nt
		comments = 'Reconstructed from'
		axes=list(self.axes)

		if mssa is True and self.steof == []: # Want MSSA and didn't run it!
			raise Exception, 'Error you did not run MSSA yet'

		if svd is True and self.svdeof == []:
			raise Exception, 'Error you did not run SVD yet'

		## Check for svd
		if svd is None:
			if self.svdeof == []:
				svd = False
			elif self.steof==[]:
				svd = True

		## Check for default values for mssa and pca if not passed by user
		if mssa is None:
			if self.steof ==[]:
				mssa = False
			elif svd is False:
				mssa = True
			else:
				mssa = False

		if pca is None:
			if self.pc == []:
				pca = False
			else:
				pca = True


		# Phase reconstruction for pca or mssa
		if phases and not pca and not mssa:
			raise 'Error you did not do any PCA or MSSA!\n To do a phases analysis only use the function %s in this module.\n%s' % ('computePhases',computePhases.__doc__)

		# MSSA reconstruction
		if mssa:
			comments+=' MSSA '
			if pca:
				nspace=[self.npca,]*len(self.pdata[0])
			else:
				nspace=[]
				for i in range(len(self.pdata)):
					nspace.append(self.pdata[i].shape[0])
			if imode is None:
				imode=range(1,self.nmssa+1)
			ffrec = self._reconstruct(spanlib_fort.mssa_rec,imode,
			  self.steof,self.stpc,nspace,self.nt, self.window)

		# SVD Reconstruction
		if svd:
			comments+=' SVD '
			if pca:
				nspace=[self.npca,self.npca]
			else:
				nspace=[]
				for i in range(2):
					nspace.append(self.pdata[i].shape[0])
			if imode is None:
				imode = range(1,self.nsvd+1)
			ffrec = self._reconstruct(spanlib_fort.pca_rec,imode,
			  self.svdeof,self.svdpc,nspace,self.nt)

		# Phase composites reconstuction
		if phases:
			comments+=' Phases'
			if mssa:
				for i in range(len(self.pdata)):
					ffrec[i] = computePhases(ffrec[i],nphases,offset,firstphase)
			else:
				ffrec=[]
				for i in range(len(self.pdata)):
					ffrec.append(computePhases(N.transpose(self.pc[i]),\
											   nphases,offset,firstphase))

			## Replace time axis with phases axis
			ntimes=nphases
			for j in range(len(self.pdata)):
				for i in range(len(self.axes[0])):
					if axes[j][i].isTime():
						axes[j][i]=ffrec[j].getAxis(1)


		# PCA reconstruction (alone or for SVD or MSSA)
		if svd:
			nloop=2
		else:
			nloop = len(self.pdata)

		if pca:
			comments+=' PCA'
			if mssa is True or phases is True or svd is True:
				pcreconstruct=[]
				for i in xrange(nloop):
					pcreconstruct.append(N.transpose(ffrec[i]))
				del(ffrec)
			else:
				pcreconstruct = self.pc
			if mssa or imode is None:
				imode = range(1,self.npca+1)
			ffrec = self._reconstruct(spanlib_fort.pca_rec,imode,self.eof,pcreconstruct)


		# Add meta information to output field
		for i in xrange(nloop):

			# Back to channel dimensions
			sh = [ffrec[i].shape[1],]
			sh.extend(self.shapes[i][1:])
			if self.mask[i] is not None:
				print ffrec[i].shape
				ffrec[i] = MV.transpose(spanlib_fort.chan_unpack(self.mask[i],ffrec[i],1.e20))
				print 'shapes:',self.shapes[i],ffrec[i].shape
				ffrec[i] = MV.reshape(ffrec[i],sh)
				tmask = N.resize(N.transpose(self.mask[i]),ffrec[i].shape)
				ffrec[i] = MV.masked_where(N.equal(tmask,0),ffrec[i])
				del tmask
			else:
				ffrec[i] = MV.transpose(ffrec[i])
				ffrec[i] = MV.reshape(ffrec[i],sh)
			ffrec[i].setAxisList(axes[i])
			ffrec[i].id=self.varname[i]+'_rec'
			ffrec[i].name = ffrec[i].id
			for att in 'units','long_name':
				if att in self.attributes[i].keys():
					if att is 'long_name':
						ffrec[i].attributes[att] = \
						  'Reconstruction of '+self.attributes[i][att]
					else:
						ffrec[i].attributes[att] = self.attributes[i][att]
			ffrec[i].comment=comments
			ffrec[i].modes = imode
			if not svd:
				ffrec[i].setGrid(self.grids[i])

		if len(ffrec)==1:
			return ffrec[0]
		else:
			return ffrec

	def clean(self):
		self.pc=[]
		self.eof=[]
		self.stpc=[]
		self.steof=[]
		self.svdpc=[]
		self.svdeof=[]
		self.l2r=[]
		self.pairs = []


class SVDModel(SpAn):

	def __init__(self,data,**kwargs):

		SpAn.__init__(self,data,**kwargs)

		# Perform an SVD between the first two datasets
		self.svd(nsvd=None,pca=None)

		# Compute the scale factors between the two datasets
		self.scale_factors = N.sqrt((N.average(self.svdpc[0]**2)/nsr - \
									 (N.average(self.svdpc[0])/nsr)**2) / \
									(N.average(self.svdpc[1]**2)/nsr - \
									 (N.average(self.svdpc[1])/nsl)**2))

	def __call__(self,data,nsvdrun=None):
		""" Run the SVD model """

		if nsvdrun is not None:
			self.nsvdrun = nsvdrun

		#TODO: finish the svd model man !
		print 'missing code'





