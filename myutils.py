""" utility functions """


import sys
import numpy as np
from scipy.stats import gaussian_kde,norm
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse

def ADSurvFunc(AD):
  return 1-norm.cdf(np.log(AD),loc=-0.22,scale=0.66);

def P2Sig(pval):
  """convert pval to sigma"""
  return norm.ppf(1.0-pval/2)

def AD2Sig(AD):
  """convert AndersonDarling TS to sigma"""
  AD=np.array(AD)
  sig=P2Sig(ADSurvFunc(AD))
  sig[AD>5]=(np.log(AD[AD>5])+0.22)/0.66
  return sig
  
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
  
def percentile_contour(data, nbin=100, percents=0.683, color='r', logscale=False, **kwargs):
    """
    plot contour at specific percentile levels
    
    percents can be a list, specify the contour percentile levels
    data should be shape [2,n] array
    logscale: default False; whether to plot in linear or logspace
    **kwargs specify linestyles
    return a handle artist of the same linestyle (but not the contour object) to be used in legends
    """
    if logscale:
      data=np.log(data)
    l=data.min(axis=1)
    r=data.max(axis=1)
    X, Y = np.mgrid[l[0]:r[0]:nbin*1j, l[1]:r[1]:nbin*1j]
    positions =np.vstack([X.ravel(), Y.ravel()])
    kernel = gaussian_kde(data)
    Z = np.reshape(kernel(positions).T, X.shape)
    lvls=percent2level(percents,Z)
    if logscale:
      h0=contour(np.exp(X),np.exp(Y),Z,lvls, colors=color, **kwargs)
      loglog()
    else:
      h0=contour(X,Y,Z,lvls, colors=color, **kwargs)
    h=Ellipse((0,0),0,0,fill=False, color=color, **kwargs)
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
	
	[count,xbin]=np.histogram(x,nbin)
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
			ylim[:,i]=[NaN,NaN]
			
	return {'x':{'median':xmed,'mean':xm,'bin':xbin,'hist':count},'y':{'median':ymed,'mean':ym,'std':ysig,'CI':ylim}}
