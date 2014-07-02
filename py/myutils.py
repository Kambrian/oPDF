""" utility functions """
import sys
import numpy as np
from scipy.stats import gaussian_kde,norm
import ctypesGsl as cgsl
import matplotlib
#matplotlib.user('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
matplotlib.rcParams.update({'font.size': 18, 'axis.labelsize': 20, 'legend.fontsize': 18, 'ps.fonttype' : 42 , 'pdf.fonttype' : 42 ,' image.origin': 'lower', 'image.interpolation': 'None'})

from matplotlib.ticker import MaxNLocator # added 

def create31fig(sharex=True, sharey=False, figsize=(8,8)):
    f,ax = plt.subplots(3, sharex=sharex, sharey=sharey, figsize=figsize)
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    nbins = len(ax[0].get_yticklabels())
    plt.setp([a.yaxis for a in ax[:-1]], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
    return ax

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
      Credit: http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
  
def plot_circle(cen=[0,0], r=1, **kwargs):
	phi=np.arange(0,2*np.pi+0.11,0.1)
	x=cen[0]+r*np.cos(phi)
	y=cen[1]+r*np.sin(phi)
	h=plt.plot(x,y,**kwargs)
	return h

def lnADPDF(lnAD):
  gpar=[[0.569,   -0.570,    0.511], [ 0.431,    0.227,    0.569]] #w, mu, sigma, bi-normal fit
  p=gpar[0][0]*norm.pdf(lnAD,gpar[0][1],gpar[0][2])+gpar[1][0]*norm.pdf(lnAD, gpar[1][1], gpar[1][2])
  return p
  
def ADSurvFunc(AD):
  gpar=[[0.569,   -0.570,    0.511], [ 0.431,    0.227,    0.569]] #w, mu, sigma, bi-normal fit
  lnAD=np.log(AD)
  p=gpar[0][0]*norm.sf(lnAD,gpar[0][1],gpar[0][2])+gpar[1][0]*norm.sf(lnAD, gpar[1][1], gpar[1][2])
  return p
  #return norm.sf(np.log(AD),loc=-0.22,scale=0.66);

def P2Sig(pval):
  """convert pval to sigma"""
  return norm.ppf(1.0-pval/2)

def AD2Sig(AD):
  """convert AndersonDarling TS to sigma"""
  AD=np.array(AD)
  sig=P2Sig(ADSurvFunc(AD))
  #sig[AD>5]=(np.log(AD[AD>5])+0.22)/0.66
  return sig


def fmin_gsl(func, x0, args=[], xtol=1e-3, ftolabs=0.01, xstep=1.0, maxiter=1000, full_output=0):
    '''
    minimize function with gsl_simplex method
    func(x [,args]): function to be minimized
    x0: [a,b,c...], initial parameter
    args: list of additional parameter if any
    xstep: initial simplex size
    mimics scipy.optimize.fmin() interface
    '''
    args=list(args)
    if args is []:
	  myfunc=lambda x,arg: func(x)
    else:
	  myfunc=lambda x,arg: func(x, *arg)
    F = cgsl.gsl_multimin_function(myfunc, len(x0), args)
    x = cgsl.vector(x0)
    T = cgsl.multimin_fminimizer_nmsimplex
    s = cgsl.multimin_fminimizer(T, F)
    s.init(x, cgsl.vector([xstep] * F.n))
    it = 0
    f1=F(x)
    while True:
        it += 1
        s.iterate()
        f0=f1
        f1=s.minimum()
        status = s.test_size(xtol)
        xx = s.x()
        if status and abs(f1-f0)<ftolabs:
            print "Optimization terminated successfully."
            print "\t Current function value: ", f1
            print "\t Iterations: ", it
            print "\t x abs err: ", s.size()
            print "\t", xx
            status=0
            break
        if it >= maxiter:
            status=1
            break
    x=np.array([xx[i] for i in xrange(F.n)])
    if full_output:
      return x,f1,it,status
    else:
      return x
  
class ProgressMonitor:
	"""monitor progress of your loops"""
	def __init__(self,total_steps, total_show=100, init_show=0):
		"""init_show: initial value progress percentage, set to 0 if no reason
		total_steps: maximum iteration steps to monitor
		total_show: number of revealing times
		"""
		self.current_show=int(init_show)
		self.total_steps=total_steps
		self.total_show=total_show
		print " %02d%%"%(self.current_show*100/self.total_show),
		sys.stdout.flush()
		
	def monitor_progress(self,current_step):
		"""put this inside loop to monitor progress, to print the percent of
		job finished."""
		#print when current_step first exceeds current show percent
		if current_step>=self.total_steps*self.current_show/self.total_show:
			print "\b\b\b\b\b %02d%%"%(self.current_show*100/self.total_show),
			sys.stdout.flush()
			self.current_show+=1
			
def percent2level(p,z):
    ''' convert percentiles to levels '''
    try:
      p=list(p)
    except: #single number
      p=[p]
    x=np.sort(np.array(z.ravel()))
    x=x[::-1]
    frac=x.cumsum()/x.sum()
    l=[x[abs(frac-pi).argmin()] for pi in p]
    return l
  
def percentile_contour(data, nbin=100, percents=0.683, colors=None, logscale=False, **kwargs):
    """
    plot contour at specific percentile levels
    
    percents can be a list, specify the contour percentile levels
    data should be shape [2,n] array
    logscale: default False; whether to plot in linear or logspace
    colors should be a tuple, e.g, (r,)
    **kwargs specify linestyles
    return a handle artist of the same linestyle (but not the contour object) to be used in legends
    """
    if logscale:
      data=np.log10(data)
    l=data.min(axis=1)
    r=data.max(axis=1)
    X, Y = np.mgrid[l[0]:r[0]:nbin*1j, l[1]:r[1]:nbin*1j]
    positions =np.vstack([X.ravel(), Y.ravel()])
    kernel = gaussian_kde(data)
    Z = np.reshape(kernel(positions).T, X.shape)
    lvls=percent2level(percents,Z)
    #if logscale:
      #h0=plt.contour(np.exp10(X),np.exp10(Y),Z,lvls, colors=colors, **kwargs)
      #plt.loglog()
    #else:
    h0=plt.contour(X,Y,Z,lvls, colors=colors, **kwargs)
    h=Ellipse((0,0),0,0,fill=False, color=list(colors)[0], **kwargs)
    return h,h0
  
def plot_cov_ellipse(cov, pos, nstd=1, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 1 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = gca()

    TS={1:2.3, 2:6.18, 3:11.8} #the TS=x'* at the specific nstd, for a 2-d
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * np.sqrt(TS[nstd]*vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, fill=False, **kwargs)
    #ellip.set_facecolor('none')
    ax.add_artist(ellip)
    return ellip
  
def skeleton(x,y,nbin=10,alpha=0.683):
	"""function [xmed,ymed,ylim,xm,ym,ysig,count]=skeleton(x,y,nbin,alpha)
	% to divide x into bins and give estimation of center and variance of y
	% inside each bin
	%input:
	% x,y: column vectors to extract skeleton from
	% nbin: number of bins or bin edges for x
	% alpha: confidence level for boundary estimation
	%"""

	x=np.array(x)
	y=np.array(y)
	
	count,xbin=np.histogram(x,nbin)
	nbin=len(xbin)-1
	bin=np.digitize(x,xbin)-1
	
	xm=np.empty(nbin)
	#ym=xm[:]     #this is wrong! even though id(ym)!=id(xm), and id(ym[0])!=id(xm[0])
	ym=np.empty_like(xm)
	ysig=np.empty_like(xm)
	xmed=np.empty_like(xm)
	ymed=np.empty_like(xm)
	ylim=np.empty([2,nbin]);
	alpha=(1-alpha)/2;
	
	for i in xrange(nbin):
		xm[i]=np.mean(x[bin==i])
		xmed[i]=np.median(x[bin==i])
		ym[i]=np.mean(y[bin==i])
		ymed[i]=np.median(y[bin==i])
		ysig[i]=np.std(y[bin==i])
		tmp=np.sort(y[bin==i])
		if count[i]:
			ylim[:,i]=[tmp[np.ceil(alpha*count[i])],tmp[np.ceil((1-alpha)*count[i])]]
		else:
			ylim[:,i]=[np.NaN,np.NaN]
			
	return {'x':{'median':xmed,'mean':xm,'bin':xbin,'hist':count},'y':{'median':ymed,'mean':ym,'std':ysig,'CI':ylim}}
