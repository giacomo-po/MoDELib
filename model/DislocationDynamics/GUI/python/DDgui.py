#simple GUI

#import everything from Tkinter
from Tkinter import *
from ttk import *

#def write_help_text(lab,text):
#    lab["text"] = text

class Simulation_Tab(Frame):
    
    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()
        self.name='Simulation'
        
        self.stringFields=['dx','Nsteps']
#        self.helpFields=['maximum distance that a dislocation node can travel in one simulation time step',
#                         'number of time steps to be performed']
        self.labelList=[]
        self.entryList=[]
        self.buttons=[]
        
        self.helpbox=Label(self)
        self.helpbox["text"] = 'help'
        self.helpbox.grid(row=0,column=3,rowspan=len(self.stringFields))

        r=0
        for x in self.stringFields:
            self.labelList.append(Label(self))
            self.labelList[r]["text"] = x
            self.labelList[r].grid(row=r,column=0,columnspan=1)
            self.entryList.append(Entry(self))
            #self.entryList[r]["text"] = x
            self.entryList[r].grid(row=r,column=1,columnspan=1)
            self.buttons.append(Button(self))
            self.buttons[r]["text"] = '?'
            #            self.buttons[r]["command" ] = write_help_text(self.helpbox,self.stringFields[r])
            self.buttons[r].grid(row=r,column=2,columnspan=1)
            r +=1


class Material_Tab(Frame):

    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()
        self.name='Material'
        
        self.stringFields=['atomic number']
        self.labelList=[]
        self.entryList=[]
        
        r=0
        for x in self.stringFields:
            self.labelList.append(Label(self))
            self.labelList[r]["text"] = x
            self.labelList[r].grid(row=r,column=0,columnspan=1)
            self.entryList.append(Entry(self))
            #self.entryList[r]["text"] = x
            self.entryList[r].grid(row=r,column=1,columnspan=1)
            r +=1

class Write_Tab(Frame):

    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()
        self.name='DDinput'

        self.button1 = Button(self)
        self.button1["text"] = "Write DDinput"
        self.button1["command" ] = self.write_file
        self.button1.grid(row=0,column=1,columnspan=1, sticky=W)

    def write_file(self):
        file = open('DDinput.txt', 'w')
        r=0
#        for x in self.stringFields:
#            file.write(x + "=" + self.entryList[r]["text"] + ";\n")
#            r+=1

class DDgui(Frame):

    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()

        self.nb = Notebook(self)
        self.nb.pack(fill='both', expand='yes')

        self.tabs=[Simulation_Tab(self.nb), Material_Tab(self.nb), Write_Tab(self.nb)]
        for x in self.tabs:
            self.nb.add(x, text=x.name)


#create the window
root = Tk()
root.title("DDgui")
root.geometry("400x600")

ddgui = DDgui(root)
ddgui.grid()
ddgui.mainloop()
