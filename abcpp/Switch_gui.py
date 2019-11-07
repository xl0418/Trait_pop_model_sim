import sys
import tkinter as tk
from tkinter import font  as tkfont
import tkinter.filedialog as fd
import tkinter.ttk as ttk
import threading
import shutil
import os
import time
from PIL import Image


sys.path.append('C:/Liang/Trait_pop_model_sim/abcpp')
# from sim_argu_test import simtest
from Trait_simulator import trait_simulator
from Continue_MS_function import Continue_trait_simulator
from Generate_plots import generate_plots
class SampleApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.title("Trait player")
        self.minsize(730, 750)
        self.title_font = tkfont.Font(family='Helvetica', size=18, weight="bold", slant="italic")

        # the container is where we'll stack a bunch of frames
        # on top of each other, then the one we want visible
        # will be raised above the others
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (ParaInf, ConParaInf, GeneCluScri,Plots):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame

            # put all of the pages in the same location;
            # the one on the top of the stacking order
            # will be the one that is visible.
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame("ParaInf")

    def show_frame(self, page_name):
        '''Show a frame for the given page name'''
        frame = self.frames[page_name]
        frame.tkraise()


class ParaInf(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        button1 = tk.Button(self, text="Parameter Inference",
                            command=lambda: controller.show_frame("ParaInf"))
        button2 = tk.Button(self, text="Continue Parameter Inference",
                            command=lambda: controller.show_frame("ConParaInf"))
        button3 = tk.Button(self, text="Generate Cluster Scripts",
                            command=lambda: controller.show_frame("GeneCluScri"))
        button4 = tk.Button(self, text="Plotting",
                            command=lambda: controller.show_frame("Plots"))
        button1.place(x = 15, y = 15, width=150, height=20)
        button2.place(x = 165, y = 15, width=200, height=20)
        button3.place(x = 365, y = 15, width=150, height=20)
        button4.place(x = 515, y = 15, width=100, height=20)

        label = tk.Label(self, text="Parameter Inference", font=controller.title_font)
        label.place(x = 15, y = 45, width=250, height=20)

        # Frame 1
        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.place(x = 15, y = 75, width=300, height=125)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A Dir",
                                      command=self.fileDialog)
        self.browsebutton.place(x = 150, y = 5, width=140, height=20)
        self.data = tk.Label(self.labelFrame1, text='Data directory').place(x = 5, y = 5,
                                                                           width=140, height=20)
        self.output = tk.Label(self.labelFrame1, text='Output file name').place(x = 5, y = 55,
                                                                           width=140, height=20)
        self.enteroutput = tk.Entry(self.labelFrame1)
        self.enteroutput.place(x = 150, y = 55, width=110, height=20)
        self.enteroutput.bind('<Return>', lambda event:self.update_output())

        self.output_set = tk.Button(self.labelFrame1, text="Set", command=self.update_output)
        self.output_set.place(x = 260, y = 55, width=30, height=20)

        # Frame 2
        # choose the number of iterations and particles
        self.labelFrame2 = tk.LabelFrame(self, text="Structure of the algorithm")
        self.labelFrame2.place(x = 15, y = 200, width=300, height=75)
        default_iterations = tk.DoubleVar(value=10)  # default value for the iterations
        self.iter = tk.Label(self.labelFrame2, text='Iterations').place(x = 5, y = 5,
                                                                        width=140, height=20)
        self.iterations = tk.Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations, textvariable=default_iterations)
        self.iterations.place(x = 150, y = 5, width=140, height=20)
        self.iterations.bind('<Return>', lambda event:self.update_iterations())

        default_particles = tk.DoubleVar(value=1000)  # default value for the particles
        self.par = tk.Label(self.labelFrame2, text='Particles').place(x = 5, y = 30, width=140,
                                                                     height=20)
        self.particles = tk.Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles, textvariable=default_particles)
        self.particles.place(x = 150, y = 30, width=140, height=20)
        self.particles.bind('<Return>', lambda event:self.update_particles())

        # Frame 3
        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads setting")
        self.labelFrame3.place(x = 325, y = 200, width=170, height=50)
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').place(x = 5, y = 5, width=50,
                                                                   height=20)
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.place(x = 60, y = 5, width=90, height=20)
        self.thread.bind('<Return>', lambda event:self.update_threads())
        # Frame 4
        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.place(x = 505, y = 200, width=210, height=50)
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1, command=self.update_sstats).place(x = 5,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).place(x = 70,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).place(x = 135,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)

        # Frame 5
        # Output info
        self.labelFrame5 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame5.place(x = 15, y = 280, width=700, height=450)
        self.progressbar = ttk.Progressbar(self.labelFrame5, mode='indeterminate',length=500)
        self.progressbar.place(x = 5, y = 5, width=400, height=20)

        # self.report = Text(master)
        # Set default values
        self.output_value = self.enteroutput.get()
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.iterations_value = self.iterations.get()
        self.particles_value = self.particles.get()

        self.gobutton = tk.Button(self.labelFrame5, text="Go", command=self.go_process, height=1,
                               width=10)
        self.gobutton.place(x = 630, y = 5, width=50, height=20)
        self.outputtextbox = tk.Text(self.labelFrame5)
        self.outputtextbox.place(x = 5, y = 30, width=680, height=390)

        # Open a file dialog

    def fileDialog(self):
        self.filename = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 30, width=285, height=20)
        self.label.configure(text=self.filename)
        print('Data file is in the directory: %s' % self.filename)

    def trait_simer(self):
        self.sstats_value = self.choice.get()
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        print(self.filename,self.output_value,self.threads_value,
              self.sstats_value,self.iterations_value,self.particles_value)
        print(type(self.filename),type(self.output_value),type(self.threads_value),
              type(self.sstats_value),type(self.iterations_value),type(self.particles_value))
        out_put = trait_simulator(files=self.filename, result=self.output_value,
                          num_threads=int(self.threads_value), sstats=stats,
                          iterations=int(self.iterations_value), particles=int(
                self.particles_value))
        self.outputtextbox.insert(tk.END, str(out_put) + '\n')

    def sim_program_thread(self):
        thread = threading.Thread(None, self.trait_simer, None, (), {})
        thread.start()

    def update_output(self):
        self.output_value = self.enteroutput.get()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 80, width=285, height=20)
        self.label.configure(text=self.output_value)
        print('The result will be stored in: %s' % self.output_value)

    def update_sstats(self):
        self.sstats_value = self.choice.get()
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        print('The summary stats is chosed to be: %s' % stats)

    def update_threads(self):
        self.threads_value = self.thread.get()
        print('No. of threads to be used: %s' % str(self.threads_value))

    def update_iterations(self):
        self.iterations_value = self.iterations.get()
        print('No. of iterations to be used: %s' % str(self.iterations_value))

    def update_particles(self):
        self.particles_value = self.particles.get()
        print('No. of particles to be used: %s' % str(self.particles_value))

        # put the test program in a seperate thread so it doesn't lock up the GUI

    # def log_program_thread(self):
    #     thread = threading.Thread(None, self.progress, None, (), {})
    #     thread.start()

    def progress(self):
        self.outputtextbox.delete('1.0', tk.END)
        # read the data
        with open("ParameterInference_log.txt", "r") as f:
            self.outputtextbox.insert(tk.INSERT, f.read())

        self.master.after(10000, self.progress)
        self.outputtextbox.see(tk.END)

    def go_process(self):
        self.progress()
        self.sim_program_thread()

    def foo(self):
        time.sleep(60)  # simulate some work

    def start_foo_thread(self):
        global foo_thread
        foo_thread = threading.Thread(target=self.foo)
        foo_thread.daemon = True
        self.progressbar.start()
        foo_thread.start()
        self.master.after(20, self.check_foo_thread)

    def check_foo_thread(self):
        if foo_thread.is_alive():
            self.master.after(20, self.check_foo_thread)
        else:
            self.progressbar.stop()


class ConParaInf(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        button1 = tk.Button(self, text="Parameter Inference",
                            command=lambda: controller.show_frame("ParaInf"))
        button2 = tk.Button(self, text="Continue Parameter Inference",
                            command=lambda: controller.show_frame("ConParaInf"))
        button3 = tk.Button(self, text="Generate Cluster Scripts",
                            command=lambda: controller.show_frame("GeneCluScri"))
        button4 = tk.Button(self, text="Plotting",
                            command=lambda: controller.show_frame("Plots"))
        button1.place(x = 15, y = 15, width=150, height=20)
        button2.place(x = 165, y = 15, width=200, height=20)
        button3.place(x = 365, y = 15, width=150, height=20)
        button4.place(x = 515, y = 15, width=100, height=20)

        label = tk.Label(self, text="Continue Parameter Inference", font=controller.title_font)
        label.place(x = 15, y = 45, width=400, height=20)


        # Frame 1
        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.place(x = 15, y = 75, width=300, height=125)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A Dir",
                                      command=self.fileDialog)
        self.browsebutton.place(x=150, y=5, width=140, height=20)
        self.data = tk.Label(self.labelFrame1, text='Data directory').place(x=5, y=5,
                                                                            width=140, height=20)
        self.output = tk.Label(self.labelFrame1, text='Output file name').place(x=5, y=55,
                                                                                width=140,
                                                                                height=20)
        self.enteroutput = tk.Entry(self.labelFrame1)
        self.enteroutput.place(x=150, y=55, width=110, height=20)
        self.output_set = tk.Button(self.labelFrame1, text="Set", command=self.update_output)
        self.output_set.place(x=260, y=55, width=30, height=20)
        self.enteroutput.bind('<Return>', lambda event:self.update_output())

        # Frame 1_1
        self.labelFrame1_1 = tk.LabelFrame(self, text="Continuous setting")
        self.labelFrame1_1.place(x = 325, y = 75, width=300, height=125)
        self.predata = tk.Label(self.labelFrame1_1, text='Pre result directory').place(x=5, y=5,
                                                                            width=140, height=20)
        self.contiter = tk.Label(self.labelFrame1_1, text='Continue iterations').place(x=5, y=55,
                                                                                width=140,
                                                                                height=20)
        self.preresultbutton = tk.Button(self.labelFrame1_1, text="Browse A File",
                                     command=self.pre_result_file)
        self.preresultbutton.place(x=150, y=5, width=140, height=20)
        default_num_continue = tk.DoubleVar(value=10)  # default value for the iterations

        self.num_continue = tk.Spinbox(self.labelFrame1_1, from_=1, to=100,
                                     command=self.update_continue_num,
                                     textvariable=default_num_continue)
        self.num_continue.place(x=150, y=55, width=110, height=20)
        self.num_continue.bind('<Return>', lambda event:self.update_continue_num())

        # Frame 3
        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads settings")
        self.labelFrame3.place(x = 15, y = 210, width=170, height=50)
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').place(x = 5, y = 5, width=50,
                                                                   height=20)
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.place(x = 60, y = 5, width=90, height=20)
        self.thread.bind('<Return>', lambda event:self.update_threads())


        # Frame 4
        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.place(x = 195, y = 210, width=210, height=50)
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1, command=self.update_sstats).place(x = 5,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).place(x = 70,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).place(x = 135,
                                                                                 y = 5,
                                                                                 width=60,
                                                                                 height=20)

        # Output info
        self.labelFrame5 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame5.place(x = 15, y = 280, width=700, height=450)
        self.progressbar = ttk.Progressbar(self.labelFrame5, mode='indeterminate',length=500)
        self.progressbar.place(x = 5, y = 5, width=400, height=20)
        # self.report = Text(master)
        # Set default values
        self.output_value = self.enteroutput.get()
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.continue_num_value = self.num_continue.get()
        self.gobutton = tk.Button(self.labelFrame5, text="Go", command=self.conti_trait_simer,
                                  height=1,
                               width=10)
        self.gobutton.place(x = 630, y = 5, width=50, height=20)

        self.outputtextbox = tk.Text(self.labelFrame5)
        self.outputtextbox.place(x = 5, y = 30, width=680, height=390)

        # Open a file dialog

    def fileDialog(self):
        self.filename = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 30, width=285, height=20)
        self.label.configure(text=self.filename)
        print('Data file is in the directory: %s' % self.filename)

    def pre_result_file(self):
        self.pre_result = fd.askopenfilename()
        self.label = tk.Label(self.labelFrame1_1, text="")
        self.label.place(x = 5, y = 30, width=285, height=20)
        self.label.configure(text=self.pre_result)
        print('Data file is in the directory: %s' % self.pre_result)


    def conti_trait_simer(self):
        self.sstats_value = self.choice.get()
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        print(self.filename, self.output_value, self.threads_value,
              self.sstats_value)
        print(type(self.filename), type(self.output_value), type(self.threads_value),
              type(self.sstats_value))
        out_put = Continue_trait_simulator(files=self.filename, result=self.output_value,
                                  num_threads=int(self.threads_value), sstats=stats,
                                  previous_result = self.pre_result,
                                           continue_num=int(self.continue_num_value))
        self.outputtextbox.insert(tk.END, str(out_put) + '\n')

    def update_output(self):
        self.output_value = self.enteroutput.get()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 80, width=285, height=20)
        self.label.configure(text=self.output_value)
        print('The result will be stored in: %s' % self.output_value)

    def update_sstats(self):
        self.sstats_value = self.choice.get()
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        print('The summary stats is chosed to be: %s' % stats)

    def update_threads(self):
        self.threads_value = self.thread.get()
        print('No. of threads to be used: %s' % str(self.threads_value))

    def update_continue_num(self):
        self.continue_num_value = self.num_continue.get()
        print('No. of continue iterations: %s' % str(self.update_continue_num))

        # put the test program in a seperate thread so it doesn't lock up the GUI

    def test_program_thread(self):
        thread = threading.Thread(None, self.progress, None, (), {})
        thread.start()

    def progress(self):
        self.outputtextbox.delete('1.0', tk.END)
        # read the data
        with open("out.txt", "r") as f:
            self.outputtextbox.insert(tk.INSERT, f.read())

        self.master.after(10000, self.progress)
        self.outputtextbox.see(tk.END)


class GeneCluScri(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        button1 = tk.Button(self, text="Parameter Inference",
                            command=lambda: controller.show_frame("ParaInf"))
        button2 = tk.Button(self, text="Continue Parameter Inference",
                            command=lambda: controller.show_frame("ConParaInf"))
        button3 = tk.Button(self, text="Generate Cluster Scripts",
                            command=lambda: controller.show_frame("GeneCluScri"))
        button4 = tk.Button(self, text="Plotting",
                            command=lambda: controller.show_frame("Plots"))
        button1.place(x = 15, y = 15, width=150, height=20)
        button2.place(x = 165, y = 15, width=200, height=20)
        button3.place(x = 365, y = 15, width=150, height=20)
        button4.place(x = 515, y = 15, width=100, height=20)

        label = tk.Label(self, text="Generate Cluster Scripts", font=controller.title_font)
        label.place(x = 15, y = 45, width=400, height=20)

        # Frame1
        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.place(x = 15, y = 75, width=300, height=125)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A Dir",
                                      command=self.fileDialog)
        self.browsebutton.place(x = 150, y = 5, width=140, height=20)
        self.data = tk.Label(self.labelFrame1, text='Data directory').place(x = 5, y = 5,
                                                                           width=140, height=20)
        self.output = tk.Label(self.labelFrame1, text='Output directory').place(x = 5, y = 55,
                                                                           width=140, height=20)
        self.browsebutton2 = tk.Button(self.labelFrame1, text="Browse A Dir",
                                       command=self.savedir)
        self.browsebutton2.place(x = 150, y = 55, width=140, height=20)


        # Frame 2
        # choose the number of iterations and particles
        self.labelFrame2 = tk.LabelFrame(self, text="Structure of the algorithm")
        self.labelFrame2.place(x=15, y=200, width=300, height=75)
        default_iterations = tk.DoubleVar(value=10)  # default value for the iterations
        self.iter = tk.Label(self.labelFrame2, text='Iterations').place(x=5, y=5,
                                                                        width=140, height=20)
        self.iterations = tk.Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations, textvariable=default_iterations)
        self.iterations.place(x=150, y=5, width=140, height=20)
        self.iterations.bind('<Return>', lambda event:self.update_iterations())
        default_particles = tk.DoubleVar(value=1000)  # default value for the particles
        self.par = tk.Label(self.labelFrame2, text='Particles').place(x = 5, y = 30, width=140,
                                                                     height=20)
        self.particles = tk.Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles, textvariable=default_particles)
        self.particles.place(x=150, y=30, width=140, height=20)
        self.particles.bind('<Return>', lambda event:self.update_particles())


        # Frame 3
        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads setting")
        self.labelFrame3.place(x=325, y=200, width=170, height=50)
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').place(x=5, y=5, width=50,
                                                                    height=20)
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.place(x=60, y=5, width=90, height=20)
        self.thread.bind('<Return>', lambda event:self.update_threads())

        # Frame 4
        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.place(x=505, y=200, width=210, height=50)
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                         value=1, command=self.update_sstats).place(x=5,
                                                                                    y=5,
                                                                                    width=60,
                                                                                    height=20)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice,
                                         value=2,
                                         command=self.update_sstats).place(x=70,
                                                                           y=5,
                                                                           width=60,
                                                                           height=20)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice,
                                         value=3,
                                         command=self.update_sstats).place(x=135,
                                                                           y=5,
                                                                           width=60,
                                                                           height=20)

        # Output info
        self.labelFrame5 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame5.place(x=15, y=280, width=700, height=450)
        # self.report = Text(master)
        # Set default values
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.iterations_value = self.iterations.get()
        self.particles_value = self.particles.get()

        self.gobutton = tk.Button(self.labelFrame5, text="Generate", command=self.generate_all,
                                  height=20,
                               width=80)
        self.gobutton.place(x=600, y=5, width=80, height=20)

        self.outputtextbox = tk.Text(self.labelFrame5)
        self.outputtextbox.place(x=5, y=30, width=680, height=390)

        # Open a file dialog

    def fileDialog(self):
        self.treedatadir = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 30, width=285, height=20)
        self.label.configure(text=self.treedatadir)
        print('Data file is in the directory: %s' % self.treedatadir)


    def savedir(self):
        self.save_dir = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x = 5, y = 80, width=285, height=20)
        self.label.configure(text=self.save_dir)
        print('Saving dir is in the directory: %s' % self.save_dir)


    def update_sstats(self):
        self.sstats_value = self.choice.get()
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        print('The summary stats is chosed to be: %s' % stats)


    def update_threads(self):
        self.threads_value = self.thread.get()
        print('No. of threads to be used: %s' % str(self.threads_value))

    def update_iterations(self):
        self.iterations_value = self.iterations.get()
        print('No. of iterations to be used: %s' % str(self.iterations_value))

    def update_particles(self):
        self.particles_value = self.particles.get()
        print('No. of particles to be used: %s' % str(self.particles_value))

    def copy_treedata(self):
        source_file=self.treedatadir
        src_files = os.listdir(source_file)
        # create treedata folder in the destination
        des_treedata_dir = os.path.join(self.save_dir,"treedata")
        if not os.path.exists(des_treedata_dir):
            os.mkdir(des_treedata_dir)
        for file_name in src_files:
            full_file_name = os.path.join(source_file, file_name)
            des_file_name = os.path.join(des_treedata_dir, file_name)
            if os.path.isdir(full_file_name):
                shutil.copytree(full_file_name, des_file_name, symlinks=False, ignore=None)
            else:
                shutil.copy(full_file_name, des_file_name)
        print("Copy the data file from %s to %s" % (source_file,self.save_dir))


    def copy_algorithm_script(self):
        current_algorithm_script = os.getcwd()+"\\abcpp\\Trait_simulator_cluster.py"
        shutil.copy2(current_algorithm_script,self.save_dir)
        piccompute_script = os.getcwd()+"\\abcpp\\pic_compute.py"
        shutil.copy2(piccompute_script,self.save_dir)
        tpupdate_script = os.getcwd()+"\\abcpp\\tp_update_theta.py"
        shutil.copy2(tpupdate_script,self.save_dir)
        print("Copy the algorithm script to %s" % (self.save_dir))

    def create_bash_script(self):
        stats_vec = ['smtd', 'umtd', 'pics']
        stats = stats_vec[self.sstats_value - 1]
        with open(self.save_dir+'\\run.sh', 'w') as rsh:
            rsh.write('''\
#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --partition=gelifes
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=20
#SBATCH --output=MSumtd.log
#SBATCH --job-name=MSumtd

python3 Trait_simulator_cluster.py --treedata %s --result %s --num_threads %i --sstats %s 
--num_iterations %i --num_particles %i
                ''' % ('treedata\\', 'ParaInf_result_'+stats, int(self.threads_value), stats,
                       int(self.iterations_value), int(self.particles_value)))

    def generate_all(self):
        self.create_bash_script()
        self.copy_treedata()
        self.copy_algorithm_script()
        self.GenCluScript_log()
        self.outputtextbox.delete('1.0', tk.END)
        # read the data
        with open("GenerateClusterScript_log.txt", "r") as f:
            self.outputtextbox.insert(tk.INSERT, f.read())
        self.outputtextbox.see(tk.END)
        with open(self.save_dir+"/run.sh", "r") as f:
            self.outputtextbox.insert(tk.INSERT, f.read())
        self.outputtextbox.see(tk.END)

    def GenCluScript_log(self):
        output_log = sys.stdout
        f = open('GenerateClusterScript_log.txt', 'w')
        sys.stdout = f
        print("Copying the data from %s to %s ... \n" % (self.treedatadir, self.save_dir))
        print("Generating run.sh script is done!" )
        print("run.sh is ..." )
        sys.stdout = output_log
        f.close()

class Plots(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        button1 = tk.Button(self, text="Parameter Inference",
                            command=lambda: controller.show_frame("ParaInf"))
        button2 = tk.Button(self, text="Continue Parameter Inference",
                            command=lambda: controller.show_frame("ConParaInf"))
        button3 = tk.Button(self, text="Generate Cluster Scripts",
                            command=lambda: controller.show_frame("GeneCluScri"))
        button4 = tk.Button(self, text="Plotting",
                            command=lambda: controller.show_frame("Plots"))
        button1.place(x=15, y=15, width=150, height=20)
        button2.place(x=165, y=15, width=200, height=20)
        button3.place(x=365, y=15, width=150, height=20)
        button4.place(x=515, y=15, width=100, height=20)

        label = tk.Label(self, text="Generate plots", font=controller.title_font)
        label.place(x=15, y=45, width=200, height=20)

        # Frame1
        self.labelFrame1 = tk.LabelFrame(self, text="Plot setting")
        self.labelFrame1.place(x=15, y=75, width=300, height=125)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A File",
                                      command=self.fileDialog)
        self.browsebutton.place(x=150, y=5, width=140, height=20)
        self.data = tk.Label(self.labelFrame1, text='Result file').place(x=5, y=5,
                                                                            width=140, height=20)
        self.output = tk.Label(self.labelFrame1, text='Plots saving directory').place(x=5, y=55,
                                                                                width=140,
                                                                                height=20)
        self.browsebutton2 = tk.Button(self.labelFrame1, text="Browse A Dir",
                                       command=self.savedir)
        self.browsebutton2.place(x=150, y=55, width=140, height=20)


        # Output info
        self.labelFrame5 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame5.place(x=15, y=280, width=700, height=450)
        self.progressbar = ttk.Progressbar(self.labelFrame5, mode='indeterminate',length=500)
        self.progressbar.place(x = 5, y = 5, width=400, height=20)
        # self.report = Text(master)
        # Set default values


        self.gobutton = tk.Button(self.labelFrame5, text="Generate", command=self.threading_generateplots,
                                  height=20,
                                  width=80)
        self.gobutton.place(x=600, y=5, width=80, height=20)

        self.outputtextbox = tk.Text(self.labelFrame5)
        self.outputtextbox.place(x=5, y=30, width=680, height=390)

        # Open a file dialog

    def fileDialog(self):
        self.resultdir = fd.askopenfilename()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x=5, y=30, width=285, height=20)
        self.label.configure(text=self.resultdir)
        print('Data file is in the directory: %s' % self.resultdir)

    def savedir(self):
        self.save_dir = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.place(x=5, y=80, width=285, height=20)
        self.label.configure(text=self.save_dir)
        print('Saving dir is : %s' % self.save_dir)

    def generate_plots(self):
        generate_plots(result=self.resultdir,plots_dir=self.save_dir)
        # read the data
        with open("Plotting_log.txt", "r") as f:
            self.outputtextbox.insert(tk.INSERT, f.read())
        self.outputtextbox.see(tk.END)




    def threading_generateplots(self):
        self.show_progress(True)

        # Start thread to process file.
        self.thread = threading.Thread(target=self.generate_plots)
        self.thread.daemon = True # Allow the program to terminate without waiting for the thread to finish.
        self.thread.start()

        # Start checking the thread.
        self.process_file_check()


    def process_file_check(self):
        if self.thread.is_alive():
            # Thread is still running, check thread again in 10 milliseconds.
            self.after(10, self.process_file_check)

        else:
            # Thread finished, handle processed results.
            # Do something with `self.fdata`.
            self.show_progress(False)


    def show_progress(self, start):
        if start:
            self.progressbar.start()
        else:
            self.progressbar.stop()
            with Image.open(self.save_dir+'/ms_plots.png') as img:
                img.show()

if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
