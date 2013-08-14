#!/usr/bin/env python

from time import clock, sleep

import wx
from numpy import *
import os

from wx.lib.floatcanvas import FloatCanvas, NavCanvas


#import hazmat, TAP_mod

ID_ABOUT_MENU = wx.NewId()           
ID_EXIT_MENU  = wx.NewId()   
ID_ZOOM_IN_MENU = wx.NewId() 
ID_ZOOM_OUT_MENU = wx.NewId() 
ID_ZOOM_TO_FIT_MENU = wx.NewId()
ID_DRAWTEST_MENU = wx.NewId()
ID_DRAWMAP_MENU = wx.NewId()
ID_CLEAR_MENU = wx.NewId()
ID_SET_FRAMERATE_MENU = wx.NewId()

ID_OPEN = wx.NewId()
ID_RUN_MOVIE = wx.NewId()
ID_RUNONTOP_MOVIE = wx.NewId()
ID_RERUN_MOVIE = wx.NewId()
ID_PAUSE_BUTTON = wx.NewId()

colorlist = ["BLACK", "RED", "CYAN", "GREEN", "SALMON", "VIOLET"]
CurrentColor = [0]

def GetColor():
    color = colorlist[CurrentColor[0]]
    CurrentColor[0] += 1
    if CurrentColor[0] > len(colorlist):
        CurrentColor[0] = 0
    return color

def EVT_NEW_FRAME_EVENT( window, function ): 
    window.Connect( -1, -1, NEW_FRAME_EVENT, function ) 


class FrameEvent(wx.PyEvent):
    def __init__(self):
        wx.PyEvent.__init__(self)
        self.SetEventType(NEW_FRAME_EVENT)

class DrawFrame(wx.Frame):
    def __init__(self, *args, **kwargs):

        wx.Frame.__init__(self, *args, **kwargs)
        
        ## Set up the MenuBar
        
        MenuBar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(ID_OPEN, "&Open map","Open a bna file")
        wx.EVT_MENU(self, ID_OPEN,      self.Open_bna)
        file_menu.AppendSeparator()
        file_menu.Append(ID_RUN_MOVIE, "Run &Movie","Run a movie of the trajectory")
        wx.EVT_MENU(self, ID_RUN_MOVIE, self.Run_Movie)
        file_menu.Append(ID_RUNONTOP_MOVIE, "Run On Top &Movie","Run a movie of the trajectory on top of existing")
        wx.EVT_MENU(self, ID_RUNONTOP_MOVIE, self.RunOnTop_Movie)
        file_menu.Append(ID_RERUN_MOVIE, "Re Run &Movie","Re-Run the existing movie of the trajectory")
        wx.EVT_MENU(self, ID_RERUN_MOVIE, self.ReRun_Movie)

        file_menu.AppendSeparator()        
        file_menu.Append(ID_EXIT_MENU, "E&xit","Terminate the program")
        wx.EVT_MENU(self, ID_EXIT_MENU,       self.OnQuit)

        wx.EVT_MENU(self, ID_RUN_MOVIE, self.Run_Movie)


        MenuBar.Append(file_menu, "&File")
        
        
        view_menu = wx.Menu()
        view_menu.Append(ID_ZOOM_TO_FIT_MENU, "Zoom to &Fit","Zoom to fit the window")
        wx.EVT_MENU(self, ID_ZOOM_TO_FIT_MENU,self.ZoomToFit)
        view_menu.Append(ID_SET_FRAMERATE_MENU, "Set Frame &Rate","Set the Frame Rate for Movie playback")
        wx.EVT_MENU(self, ID_SET_FRAMERATE_MENU,self.SetFrameRate)


        MenuBar.Append(view_menu, "&View")
        
        help_menu = wx.Menu()
        help_menu.Append(ID_ABOUT_MENU, "&About",
                                "More information About this program")
        wx.EVT_MENU(self, ID_ABOUT_MENU,      self.OnAbout)
        MenuBar.Append(help_menu, "&Help")
        
        self.SetMenuBar(MenuBar)
        
                
        self.CreateStatusBar()
        self.SetStatusText("")
        
        wx.EVT_CLOSE(self, self.OnCloseWindow)
        
        # Add the Canvas

        self.NavCanvas = NavCanvas.NavCanvas(self,-1,(500,500),
                                  ProjectionFun = 'FlatEarth',
                                  Debug = 0,
                                  #BackgroundColor = "DARK SLATE BLUE")
                                  BackgroundColor = "WHITE",
                                  #UseBackground = 1,
                                  ).Canvas
        self.Canvas = NavCanvas.Canvas
        self.Canvas.NumBetweenBlits = 20

        tb = self.NavCanvas.ToolBar
        tb.AddSeparator()

        RewindButton = wx.Button(tb, -1, "Rewind")
        tb.AddControl(RewindButton)
        wx.EVT_BUTTON(self, RewindButton.GetId() , self.Rewind)

        StopButton = wx.Button(tb, -1, "Stop")
        tb.AddControl(StopButton)
        wx.EVT_BUTTON(self, StopButton.GetId() , self.Stop)

        PlayButton = wx.Button(tb, -1, "Play")
        tb.AddControl(PlayButton)
        wx.EVT_BUTTON(self, PlayButton.GetId() ,self.Play)

        tb.Realize()

        self.Show(True)
        
        self.LE_movie = None
        self.LEsObjects = []

        self.TimeStep = 0

        self.FrameDelay = 10 # milliseconds
        
        self.FileDialog = wx.FileDialog(self, "Pick a  file",".","","*",wx.OPEN)

        self.Timer = wx.PyTimer(self.ShowFrame)

        return None

    def Open_bna(self, event):
        dlg = self.FileDialog
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.LoadMap(filename)

    def LoadMap(self, filename):
        self.Canvas.Clear()
        try:
            shorelines = hazmat.read_bna(filename,polytype = "PolygonSet")
            for shoreline in shorelines:
                self.Canvas.AddPolygon(shoreline,
                                       LineWidth = 1,
                                       LineColor = "Black",
                                       FillColor = "Brown",
                                       FillStyle = 'Solid',
                                       Foreground = 0)
            self.Canvas.ZoomToBB()
        except:
            dlg = wx.MessageDialog(self, 'There was something wrong with the selected map file',
                                  'View Trajectories', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

    def Load_Movie(self, event):
        import glob
        dlg = self.FileDialog
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            
            (self.LE_movie,(NumTimesteps,NumLEs),HeaderData,flags) = TAP_mod.ReadTrajectory(filename)
            wx.GetApp().Yield()
            return True
        else:
            return None
                
    def Run_Movie(self, event):
        if self.Load_Movie(None):
            if self.LEsObjects:
                self.Canvas.RemoveObjects(self.LEsObjects)
                self.LEsObjects = []
            self.LEsObjects.append(self.Canvas.AddPointSet(self.LE_movie[0], Color = "Black", Diameter = 1.5,Foreground = 1))
            CurrentColor[0] = 1
            self.ReRun_Movie(None)

    def RunOnTop_Movie(self, event):
        if self.Load_Movie(None):
            for object in self.LEsObjects:
                object.PutInBackground()
            self.LEsObjects.append(self.Canvas.AddPointSet(self.LE_movie[0], Color = GetColor(), Diameter = 1.5,Foreground = 1) )
            self.ReRun_Movie(None)

    def ReRun_Movie(self, event):
        if not self.LE_movie:
            self.Run_Movie(None)
        else:
            self.Play(None)

##    def UpdateThread(self):
##        try:
##            while hasattr(self, 'event') and not self.event.isSet():
##                wx.PostEvent(self, FrameEvent())
##                self.event.wait(self.FrameDelay)
##        except wx.PyDeadObjectError: # BUG: we were destroyed
##            return

    def Running(self):
        """Returns true if the animation is running"""
        return self.Timer.IsRunning()

    def Play(self,event):
        """Start the animation"""
        
        if not self.Running():
            if self.LE_movie:
                #self.event.clear()
                #thread = threading.Thread(target = self.UpdateThread)
                #thread.start()
                self.Timer.Start(self.FrameDelay)
            else:
                self.Run_Movie(None)

    def Stop(self,event):
        self.Timer.Stop()

    def ShowFrame(self):
            if  self.TimeStep < len(self.LE_movie):
                self.SetStatusText("Timestep # %i of %i"%(self.TimeStep+1,len(self.LE_movie)))
                # this sets the data for the next frame
                self.LEsObjects[-1].SetPoints(self.LE_movie[self.TimeStep])
                self.Canvas.Draw()
                self.TimeStep += 1
                wx.GetApp().Yield(True)
            else:
                self.Timer.Stop()


    def Rewind(self,event):
        self.TimeStep = 0
        if self.LE_movie:
            self.LEsObjects[-1].SetPoints(self.LE_movie[self.TimeStep])
            self.SetStatusText("Timestep # %i of %i"%(self.TimeStep+1,len(self.LE_movie)))
            self.Canvas.Draw()


    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "This is a small program to demonstrate\n"
                                                  "the use of the FloatCanvas\n",
                                                  "About Me", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
    def ZoomToFit(self,event):
        self.Canvas.ZoomToBB()
        
    def Clear(self,event = None):
        self.Canvas.Clear()
        self.Canvas.Draw()
        
    def OnQuit(self,event):
        self.Close(True)
        
    def OnCloseWindow(self, event):
        self.Destroy()
        
    def RunMovie(self,event = None):
        import RandomArray
        start = clock()
        shift = RandomArray.randint(0,0,(2,))
        NumFrames = 50
        for i in range(NumFrames):
            points = self.LEs.Points
            shift = RandomArray.randint(-5,5,(2,))
            points += shift
            self.LEs.SetPoints(points)
            self.Canvas.Draw()
        print "running the movie took %f seconds to disply %i frames"%((clock() - start),NumFrames)

    def SetFrameRate(self,event):
        dlg = wx.TextEntryDialog(self,
                                'Please set the time between frames in milliseconds',
                                'ViewTrajectories',
                                "%i"%self.FrameDelay)
        dlg.SetValue("%i"%self.FrameDelay)
        if dlg.ShowModal() == wx.ID_OK:
            try:
                self.FrameDelay = int(dlg.GetValue())
            except:
                pass
        dlg.Destroy()
        
class TrajectoryViewer(wx.App):
    """

    Any bugs, comments, feedback, questions, and especially code are welcome:
    
    -Chris Barker
    
    Chris.Barker@noaa.gov
    
    """
    
    def OnInit(self):
        frame = DrawFrame(None, title="Trajectory Viewer", size=(700,700))

        self.SetTopWindow(frame)

        return True
            


    
if __name__ == "__main__":

    app = TrajectoryViewer(0)

    app.MainLoop()
    
    
    
    
    
    
    
    
    








