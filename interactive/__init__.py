# 

from eyeSpec.dependencies import check_for_wx #@UnresolvedImport
check_for_wx()

import wx #@UnusedImport
from wx.lib.newevent import NewEvent #@UnusedImport
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg #@UnusedImport
FigureCanvas = FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx #@UnusedImport
