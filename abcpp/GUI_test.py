from tkinter import *
from tkinter import filedialog
sys.path.append('C:/Liang/Trait_pop_model_sim/abcpp')
from sim_argu_test import simtest
class trait_player:
    def __init__(self, master):
        self.master = master
        master.title("Parameter estimation for baleen whales")
        master.minsize(640, 400)
        # Label frame for the input and output info
        self.labelFrame1 = LabelFrame(master, text="Input and Output")
        self.labelFrame1.grid(column=0, row=0, padx=20, pady=20,sticky="W",columnspan=4)
        self.browsebutton = Button(self.labelFrame1, text="Browse A File", command=self.fileDialog)
        self.browsebutton.grid(column=1, row=0)
        self.data = Label(self.labelFrame1, text='Data directory').grid(row=0,sticky ="W")
        self.output = Label(self.labelFrame1, text='Output file name').grid(row=2,sticky ="W")
        self.enteroutput = Entry(self.labelFrame1)
        self.enteroutput.grid(row=2, column=1,ipadx ="10")
        self.output_set = Button(self.labelFrame1,text="Set",command=self.update_output)
        self.output_set.grid(column=2,row=2)

        # choose the number of iterations and particles
        self.labelFrame2 = LabelFrame(master, text="Structure of the algorithm")
        self.labelFrame2.grid(column=0, row=1, padx=20, pady=20,sticky="W")
        default_iterations = DoubleVar(value=10)   # default value for the iterations
        self.iter = Label(self.labelFrame2, text='Iterations').grid(row=2,sticky ="W")
        self.iterations = Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations,textvariable=default_iterations)
        self.iterations.grid(row=2, column=1,ipadx ="9")
        default_particles = DoubleVar(value=1000)   # default value for the particles
        self.par = Label(self.labelFrame2, text='Particles').grid(row=3,sticky ="W")
        self.particles = Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles,textvariable=default_particles)
        self.particles.grid(row=3, column=1,ipadx ="9")

        # choose the number of threads
        self.labelFrame3 = LabelFrame(master, text="Threads settings")
        self.labelFrame3.grid(column=2, row=1, padx=20, pady=20,sticky="N")
        default_threads = DoubleVar(value=1)   # default value for the threads
        self.thr = Label(self.labelFrame3, text='Threads').grid(row=4,sticky ="W")
        self.thread = Spinbox(self.labelFrame3, from_=1, to=8,command = self.update_threads,
                              textvariable=default_threads)
        self.thread.grid(row=4, column=1,ipadx ="9")

        # choose one summary stats
        self.labelFrame4 = LabelFrame(master, text="Summary statistics")
        self.labelFrame4.grid(column=3, row=1, padx=20, pady=20,sticky="N")
        self.choice = IntVar()
        self.choice.set(1)
        self.smtdbutton = Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1,command=self.update_sstats).grid(row=0, column=0)
        self.umtdbutton = Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).grid(row=0, column=1)
        self.picsbutton = Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).grid(row=0, column=2)


        # Output info
        self.labelFrame4 = LabelFrame(master, text="Report of the progress")
        self.labelFrame4.grid(column=0, row=2, padx=20, pady=20,sticky="W",columnspan=4)
        # self.report = Text(master)
        # Set default values
        self.output_value = self.enteroutput.get()
        self.sstats_value=self.choice.get()
        self.threads_value=self.thread.get()
        self.iterations_value=self.iterations.get()
        self.particles_value=self.particles.get()


        self.gobutton = Button(self.labelFrame4, text="Go", command=self.trait_simer,height = 1,
                               width = 10)
        self.gobutton.grid(row=0, column=0,sticky="E")
        self.outputtextbox = Text(self.labelFrame4)
        self.outputtextbox.grid(column=0, row=1,sticky="W")

    # Open a file dialog
    def fileDialog(self):
        self.filename = filedialog.askdirectory()
        self.label = Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=1)
        self.label.configure(text=self.filename)
        print('Data file is in the directory: %s' % self.filename)

    def trait_simer(self):
        out_put = simtest(files=self.filename,result=self.output_value,
                        num_threads=self.threads_value ,sstats=self.sstats_value,
                        iterations=self.iterations_value,particles=self.particles_value)
        self.outputtextbox.insert(END, str(out_put) + '\n')

    def update_output(self):
        self.output_value = self.enteroutput.get()
        self.label = Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=3)
        self.label.configure(text=self.output_value)
        print('The result will be stored in: %s' % self.output_value)

    def update_sstats(self):
        self.sstats_value=self.choice.get()
        stats_vec = ['smtd','umtd','pics']
        stats = stats_vec[self.sstats_value-1]
        print('The summary stats is chosed to be: %s' % stats)

    def update_threads(self):
        self.threads_value=self.thread.get()
        print('No. of threads to be used: %s' % str(self.threads_value))

    def update_iterations(self):
        self.iterations_value=self.iterations.get()
        print('No. of iterations to be used: %s' % str(self.iterations_value))

    def update_particles(self):
        self.particles_value=self.particles.get()
        print('No. of particles to be used: %s' % str(self.particles_value))

root = Tk()
my_gui = trait_player(root)
root.mainloop()




