import matplotlib
#Set the backend to be used for the rest of ACHP
matplotlib.use('WxAgg')

import wx,os
import time,sys
import ACHPMainFrame

class MySplashScreen(wx.SplashScreen):
    """
    Create a splash screen widget.
    """
    def __init__(self, parent=None):
        # This is a recipe to a the screen.
        # Modify the following variables as necessary.
        aBitmap = wx.Image(name = os.path.join("imgs","Splash.png")).ConvertToBitmap()
        splashStyle = wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT
        splashDuration = 2 # milliseconds
        # Call the constructor with the above arguments in exactly the
        # following order.
        wx.SplashScreen.__init__(self, aBitmap, splashStyle,
                                 splashDuration, parent)
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        wx.Yield()

    def OnExit(self, evt):
        self.Hide()
        evt.Skip()  # Make sure the default handler runs too...

if __name__ == "__main__":
    app = wx.App(False)
    
    if '--nosplash' not in sys.argv:
        Splash=MySplashScreen()
        Splash.Show()
        time.sleep(1.0)
    
    MainFrame = ACHPMainFrame.ACHPMainFrame(None, -1)
    app.SetTopWindow(MainFrame)
    MainFrame.Maximize()
    MainFrame.Show()

    app.MainLoop()