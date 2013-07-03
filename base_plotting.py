
if __name__ != '__main__':
    import eyeSpec
    from eyeSpec.dependencies import np, os, time, deepcopy, pdb, pickle, plt, FormatStrFormatter, wx
    from eyeSpec.interactive_classes import Cursor
    from eyeSpec.base_functions import inv_var_2_var, alt_order_colors
    from eyeSpec.base_classes import query_fits_header

def plot_spec (spec, fig_num = None, ax=None, fig_size = (10.0,8.0), adderror=False,
               alt_color=True, addcursor=False, num_orders=False,
               ylabel='Data', zorder=0, **mpl_kwargs):
    """ 
    This uses matplotlib to create a basic figure with the data of the input spectrum
    
    INPUTS
    ===========   ==============================================================
    keyword       (class) Description
    ===========   ==============================================================
    spec          (eyeSpec_spec) This is the spectrum class for eyeSpec
    fig_num       (int) Number for matplotlib figure. Use the same number in 
                        multiple function calls to overplot.
                        
    ax            (matplotlib.axes.AxesSubplot or None) if given then it will put
                    the data onto this particular plot axes
    fig_size      (tuple) Gives the width and height of the plot, (width, height)
    adderror      (bool) If True will plot the error from the spectrum class
    alt_color     (bool) If True will use alternate order colors based on eyeSpec 
                         internal sceme
    addcursor     (bool) If True it will add a long and height wise cursor
    num_orders    (bool) If True will annotate the plot with the order number
    ylabel        (str) Give the label for the y axis
    zorder        (int) Gives the horizontal stack of matplotlib objects
    **mpl_kwargs  () Can give all the matplotlib kwargs (e.g. color = 'r',
                     linestyle='--') see help(matplotlib.pylab.plot) for more 
    ===========   ==============================================================


    OUTPUTS:
    ===   =====================================================================
    ax    () returns the subplot for the current plot
             ax = plt.figure(fig_num).add_subplot(111)
    ===   =====================================================================


    EXAMPLES:
    For plotting up a single object
    >>> x = readin("foo.fits")
    >>> plot_spec(x)
    >>> plt.show()

    For plotting multiple objects
    >>> x = readin("foo1.fits")
    >>> y = readin("foo2.fits")
    >>> plot_spec(x,10,color='r')
    >>> plot_spec(y,10,color='b')
    >>> plt.show()

    For adding various things to the plot
    >>> x = readin("foo1.fits")
    >>> ax = plot_spec(x)
    >>> ax.axhline()
    >>> ax.axvline(5555.0)
    >>> plt.show()
    
    """
    num_orders = bool(num_orders)
    # !! check spec to make sure it's of the correct class

    #if plt.matplotlib.rcParams['backend'] == 'WXAgg':
    #    wx.App(False)
    if ax is None: ax = plt.figure(fig_num,figsize = fig_size).add_subplot(111)
    elif repr(ax).find("matplotlib.axes.AxesSubplot") == -1: raise ValueError("expected ax which was a matplotlib subplot")
    
    if addcursor:
        try:
            cursor = Cursor(ax)
            cursor.connect()
        except: pass

    if 'color' in mpl_kwargs.keys(): alt_color=False

    ax.set_xlim(spec.get_wlbounds())
    ax.set_ylim(spec.get_databounds())
    ax.xaxis.set_major_formatter(FormatStrFormatter('%10.1f'))


    ax.set_xlabel("Wavelength")
    ax.set_ylabel(str(ylabel))

    title_spec = str(spec).split("\n")
    bandi = spec.get_band()
    bandid = query_fits_header(spec.header,'BANDID'+str(bandi+1))
    current_band = "Using Band: "+str(bandi)
    if bandid.found:
        current_band += "       '"+bandid.val+"'"

    title = "\n".join([title_spec[1],
                       title_spec[0],
                       current_band])
        
    ax.set_title(title)

    c_i = 0
    i = 0
    for order in spec:
        altcolor,c_i = alt_order_colors(c_i)
        c_i += 1

        if alt_color: ax.plot(order[0],order[1],zorder=zorder,color=altcolor,**mpl_kwargs)
        else: ax.plot(order[0],order[1],zorder=zorder,**mpl_kwargs)
            
        if num_orders:
            ax.set_ylim(0,spec.get_max()[1])
            ymean = np.mean(order[1])
            ymax = np.max(order[1])
            xmid = (np.max(order[0])+np.min(order[0]))/2.0
            yoff = np.abs(ymean)*0.1
            ax.plot((xmid,xmid),(0,abs(ymax)*2.0),color='k',alpha=.4)
            ax.annotate("ORDER:"+str(i),
                        (xmid,yoff),
                        (xmid,0),
                        ha='center',va='top',fontweight='bold',fontsize='small')
            i+=1
                        
        if adderror:
            inv_var = deepcopy(order[2])
            var = inv_var_2_var(inv_var)
            yerr = np.sqrt(var)
  
            ignore_ones = (inv_var != 1.0)

            X = order[0][ignore_ones]
            Y = order[1][ignore_ones]
            yerr = yerr[ignore_ones]
            if len(X) > 0:
                if alt_color: ax.errorbar(X,Y,yerr=yerr,zorder=zorder-1,color=altcolor,**mpl_kwargs)
                else: ax.errorbar(X,Y,yerr=yerr,zorder=zorder-1,**mpl_kwargs)

    return ax

plotspec = plot_spec
