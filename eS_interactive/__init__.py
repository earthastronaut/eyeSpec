# 

from eyeSpec import __path__ as _eyeSpec_path #@UnresolvedImport

from eyeSpec.dependencies import _check_for_wx #@UnresolvedImport
_check_for_wx()

import wx #@UnusedImport
from wx.lib.newevent import NewEvent as _NewEvent #@UnusedImport
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as _FigureCanvasWxAgg #@UnusedImport
_FigureCanvas = _FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as _NavigationToolbar2Wx#@UnusedImport

# import edit applications
