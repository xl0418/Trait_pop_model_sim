from tkinter import *

class trait_player:
    def __init__(self, master):
        self.master = master
        master.title("Parameter estimation for baleen whales")

        self.data = Label(master, text='Data directory').grid(row=0,sticky ="W")
        self.output = Label(master, text='Output file name').grid(row=1,sticky ="W")
        self.enterdata = Entry(master)
        self.enteroutput = Entry(master)
        self.enterdata.grid(row=0, column=1,ipadx ="10")
        self.enteroutput.grid(row=1, column=1,ipadx ="10")

        # choose the number of iterations and particles
        self.iter = Label(master, text='Iterations').grid(row=2,sticky ="W")
        self.iterations = Spinbox(master, from_=1, to=100)
        self.iterations.grid(row=2, column=1,ipadx ="9")
        self.par = Label(master, text='Particles').grid(row=3,sticky ="W")
        self.particles = Spinbox(master, from_=100, to=100000)
        self.particles.grid(row=3, column=1,ipadx ="9")

        # choose the number of threads
        self.thr = Label(master, text='Threads').grid(row=4,sticky ="W")
        self.thread = Spinbox(master, from_=1, to=8)
        self.thread.grid(row=4, column=1,ipadx ="9")

        # choose one summary stats
        self.sstats = Label(master, text='Summary stats').grid(row=5,sticky ="W")

        self.choice = IntVar()
        self.choice.set(1)
        self.smtdbutton = Radiobutton(master, text='smtd', variable=self.choice, value=1).grid(row=5,
                                                                                     column=1)
        self.umtdbutton = Radiobutton(master, text='umtd', variable=self.choice, value=2).grid(
            row=6, column=1)
        self.picsbutton = Radiobutton(master, text='pics', variable=self.choice, value=3).grid(
            row=7, column=1)

    def greet(self):
        print("Greetings!")

root = Tk()
my_gui = trait_player(root)
root.mainloop()

