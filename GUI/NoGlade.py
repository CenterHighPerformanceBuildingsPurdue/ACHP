import wx

from LoadGUI import LoadGUI
    

class Example(wx.Frame):
  
    def __init__(self, parent, title):
        super(Example, self).__init__(parent, title=title, 
            size=(800, 600))
        
        LoadGUI(self)
        
        self.Bind(wx.EVT_KEY_DOWN,self.OnKeyDown)
        
        self.Centre()
#        self.Maximize()
        self.Show()
    
    def OnKeyDown(self, event):
        keycode = event.GetKeyCode()
        print keycode
#            if keycode == wx.WXK_ESCAPE:
#                ret  = wx.MessageBox('Are you sure to quit?', 'Question', 
#            wx.YES_NO | wx.NO_DEFAULT, self)
#                if ret == wx.YES:
#                    self.Close()
        event.Skip()

if __name__ == '__main__':
  
    app = wx.App(0) # will redirect stdout to console
    Example(None, title='Review')
    app.MainLoop()