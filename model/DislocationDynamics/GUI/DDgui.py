#simple GUI

#import everything from Tkinter
from Tkinter import *

class DDgui(Frame):

    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()
        #        self.button_clicks = 0
    
        self.stringFields=['dx','Nsteps','asfsdf','Giacomo','Ben']
        self.labelList=[]
        self.entryList=[]
        self.create_widgets()
    

    def create_widgets(self):
        #create first button
        r=0
        for x in self.stringFields:
            self.labelList.append(Label(self))
            self.labelList[r]["text"] = x
            self.labelList[r].grid(row=r,column=0,columnspan=1)
            self.entryList.append(Entry(self))
            #self.entryList[r]["text"] = x
            self.entryList[r].grid(row=r,column=1,columnspan=1)
            r +=1
        self.button1 = Button(self)
        self.button1["text"] = "Write DDinput"
        self.button1["command" ] = self.write_file
        self.button1.grid(row=r,column=3,columnspan=1, sticky=W)

    def write_file(self):
        file = open('DDinput.txt', 'w')
        r=0
        for x in self.stringFields:
            file.write(x + "=" + self.entryList[r]["text"] + ";\n")
            r+=1

#create the window
root = Tk()
root.title("DDgui")
root.geometry("800x600")

ddgui = DDgui(root)
ddgui.grid()
ddgui.mainloop()