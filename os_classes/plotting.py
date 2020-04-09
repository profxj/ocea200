""" Methods and the like for plotting; mainly in Bokeh"""

def bokeh_set_fontsize(p, fsz):
    """
    Parameters
    ----------
    ax : Bokeh plot class
    fsz : float
      Font size
    """
    p.xaxis.axis_label_text_font_size = '{:d}pt'.format(fsz)
    p.xaxis.major_label_text_font_size = "{:d}pt".format(fsz)
    #
    p.yaxis.axis_label_text_font_size = '{:d}pt'.format(fsz)
    p.yaxis.major_label_text_font_size = "{:d}pt".format(fsz)
