import sys
import tkinter as tk                # python 3
from tkinter import font  as tkfont # python 3
import tkinter.filedialog as fd
import tkinter.ttk as ttk
import threading
import shutil
import os
import time

sys.path.append('C:/Liang/Trait_pop_model_sim/abcpp')
# from sim_argu_test import simtest
from Trait_simulator import trait_simulator
class SampleApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.minsize(640, 400)
        self.title_font = tkfont.Font(family='Helvetica', size=18, weight="bold", slant="italic")

        # the container is where we'll stack a bunch of frames
        # on top of each other, then the one we want visible
        # will be raised above the others
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (ParaInf, ConParaInf, GeneCluScri):
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
        button1.grid(row=0,column=0,sticky="W")
        button2.grid(row=0,column=2,sticky="W")
        button3.grid(row=0,column=3,sticky="W")

        label = tk.Label(self, text="Parameter Inference", font=controller.title_font)
        label.grid(row=1,column=0,sticky="W", columnspan=5)

        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.grid(column=0, row=2, padx=20, pady=20, sticky="W", columnspan=5)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A File", command=self.fileDialog)
        self.browsebutton.grid(column=1, row=0)
        self.data = tk.Label(self.labelFrame1, text='Data directory').grid(row=0, sticky="W")
        self.output = tk.Label(self.labelFrame1, text='Output file name').grid(row=2, sticky="W")
        self.enteroutput = tk.Entry(self.labelFrame1)
        self.enteroutput.grid(row=2, column=1, ipadx="10")
        self.output_set = tk.Button(self.labelFrame1, text="Set", command=self.update_output)
        self.output_set.grid(column=2, row=2)

        # choose the number of iterations and particles
        self.labelFrame2 = tk.LabelFrame(self, text="Structure of the algorithm")
        self.labelFrame2.grid(column=0, row=3, padx=20, pady=20, sticky="W")
        default_iterations = tk.DoubleVar(value=10)  # default value for the iterations
        self.iter = tk.Label(self.labelFrame2, text='Iterations').grid(row=2, sticky="W")
        self.iterations = tk.Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations, textvariable=default_iterations)
        self.iterations.grid(row=2, column=1, ipadx="9")
        default_particles = tk.DoubleVar(value=1000)  # default value for the particles
        self.par = tk.Label(self.labelFrame2, text='Particles').grid(row=3, sticky="W")
        self.particles = tk.Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles, textvariable=default_particles)
        self.particles.grid(row=3, column=1, ipadx="9")

        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads settings")
        self.labelFrame3.grid(column=2, row=3, padx=20, pady=20, sticky="N")
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').grid(row=4, sticky="W")
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.grid(row=4, column=1, ipadx="9")

        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.grid(column=3, row=3, padx=20, pady=20, sticky="N")
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1, command=self.update_sstats).grid(row=0, column=0)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).grid(row=0, column=1)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).grid(row=0, column=2)

        # Output info
        self.labelFrame4 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame4.grid(column=0, row=4, padx=20, pady=20, sticky="W", columnspan=5)
        self.progressbar = ttk.Progressbar(self.labelFrame4, mode='indeterminate',length=500)
        self.progressbar.grid(column=0, row=0, sticky="W")

        # self.report = Text(master)
        # Set default values
        self.output_value = self.enteroutput.get()
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.iterations_value = self.iterations.get()
        self.particles_value = self.particles.get()

        self.gobutton = tk.Button(self.labelFrame4, text="Go", command=self.go_process, height=1,
                               width=10)
        self.gobutton.grid(row=0, column=1, sticky="E")
        self.progressbutton = tk.Button(self.labelFrame4, text="Progress",
                                     command=self.progress, height=1, width=10)
        self.progressbutton.grid(row=1, column=1, sticky="E")
        self.outputtextbox = tk.Text(self.labelFrame4)
        self.outputtextbox.grid(column=0, row=2, sticky="W")

        # Open a file dialog

    def fileDialog(self):
        self.filename = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=1)
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
        self.label.grid(column=1, row=3)
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
        button1.grid(row=0, column=0, sticky="W")
        button2.grid(row=0, column=2, sticky="W")
        button3.grid(row=0, column=3, sticky="W")

        label = tk.Label(self, text="Continue Parameter Inference", font=controller.title_font)
        label.grid(row=1,column=0,sticky="W", columnspan=5)


        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.grid(column=0, row=2, padx=20, pady=20, sticky="W", columnspan=5)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A File", command=self.fileDialog)
        self.browsebutton.grid(column=1, row=0)
        self.data = tk.Label(self.labelFrame1, text='Data directory').grid(row=0, sticky="W")
        self.output = tk.Label(self.labelFrame1, text='Output file name').grid(row=2, sticky="W")
        self.enteroutput = tk.Entry(self.labelFrame1)
        self.enteroutput.grid(row=2, column=1, ipadx="10")
        self.output_set = tk.Button(self.labelFrame1, text="Set", command=self.update_output)
        self.output_set.grid(column=2, row=2)

        # choose the number of iterations and particles
        self.labelFrame2 = tk.LabelFrame(self, text="Structure of the algorithm")
        self.labelFrame2.grid(column=0, row=3, padx=20, pady=20, sticky="W")
        default_iterations = tk.DoubleVar(value=10)  # default value for the iterations
        self.iter = tk.Label(self.labelFrame2, text='Iterations').grid(row=2, sticky="W")
        self.iterations = tk.Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations, textvariable=default_iterations)
        self.iterations.grid(row=2, column=1, ipadx="9")
        default_particles = tk.DoubleVar(value=1000)  # default value for the particles
        self.par = tk.Label(self.labelFrame2, text='Particles').grid(row=3, sticky="W")
        self.particles = tk.Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles, textvariable=default_particles)
        self.particles.grid(row=3, column=1, ipadx="9")

        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads settings")
        self.labelFrame3.grid(column=2, row=3, padx=20, pady=20, sticky="N")
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').grid(row=4, sticky="W")
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.grid(row=4, column=1, ipadx="9")

        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.grid(column=3, row=3, padx=20, pady=20, sticky="N")
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1, command=self.update_sstats).grid(row=0, column=0)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).grid(row=0, column=1)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).grid(row=0, column=2)

        # Output info
        self.labelFrame4 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame4.grid(column=0, row=4, padx=20, pady=20, sticky="W", columnspan=5)
        # self.report = Text(master)
        # Set default values
        self.output_value = self.enteroutput.get()
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.iterations_value = self.iterations.get()
        self.particles_value = self.particles.get()

        self.gobutton = tk.Button(self.labelFrame4, text="Go", command=self.trait_simer, height=1,
                               width=10)
        self.gobutton.grid(row=0, column=0, sticky="E")
        self.progressbutton = tk.Button(self.labelFrame4, text="Progress",
                                     command=self.test_program_thread, height=1, width=10)
        self.progressbutton.grid(row=1, column=0, sticky="E")
        self.outputtextbox = tk.Text(self.labelFrame4)
        self.outputtextbox.grid(column=0, row=2, sticky="W")

        # Open a file dialog

    def fileDialog(self):
        self.filename = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=1)
        self.label.configure(text=self.filename)
        print('Data file is in the directory: %s' % self.filename)

    def trait_simer(self):
        out_put = simtest(files=self.filename, result=self.output_value,
                          num_threads=self.threads_value, sstats=self.sstats_value,
                          iterations=self.iterations_value, particles=self.particles_value)
        self.outputtextbox.insert(tk.END, str(out_put) + '\n')

    def update_output(self):
        self.output_value = self.enteroutput.get()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=3)
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
        button1.grid(row=0, column=0, sticky="W")
        button2.grid(row=0, column=2, sticky="W")
        button3.grid(row=0, column=3, sticky="W")

        label = tk.Label(self, text="Generate Cluster Scripts", font=controller.title_font)
        label.grid(row=1,column=0,sticky="W", columnspan=5)


        self.labelFrame1 = tk.LabelFrame(self, text="Input and Output")
        self.labelFrame1.grid(column=0, row=2, padx=20, pady=20, sticky="W", columnspan=5)
        self.browsebutton = tk.Button(self.labelFrame1, text="Browse A File", command=self.fileDialog)
        self.browsebutton.grid(column=1, row=0)
        self.data = tk.Label(self.labelFrame1, text='Data directory').grid(row=0, sticky="W")
        self.output = tk.Label(self.labelFrame1, text='Output directory').grid(row=2, sticky="W")
        self.browsebutton2 = tk.Button(self.labelFrame1, text="Browse A File",
                                       command=self.savedir)
        self.browsebutton2.grid(column=1, row=2)


        # choose the number of iterations and particles
        self.labelFrame2 = tk.LabelFrame(self, text="Structure of the algorithm")
        self.labelFrame2.grid(column=0, row=3, padx=20, pady=20, sticky="W")
        default_iterations = tk.DoubleVar(value=10)  # default value for the iterations
        self.iter = tk.Label(self.labelFrame2, text='Iterations').grid(row=2, sticky="W")
        self.iterations = tk.Spinbox(self.labelFrame2, from_=1, to=100,
                                  command=self.update_iterations, textvariable=default_iterations)
        self.iterations.grid(row=2, column=1, ipadx="9")
        default_particles = tk.DoubleVar(value=1000)  # default value for the particles
        self.par = tk.Label(self.labelFrame2, text='Particles').grid(row=3, sticky="W")
        self.particles = tk.Spinbox(self.labelFrame2, from_=100, to=100000,
                                 command=self.update_particles, textvariable=default_particles)
        self.particles.grid(row=3, column=1, ipadx="9")

        # choose the number of threads
        self.labelFrame3 = tk.LabelFrame(self, text="Threads settings")
        self.labelFrame3.grid(column=2, row=3, padx=20, pady=20, sticky="N")
        default_threads = tk.DoubleVar(value=1)  # default value for the threads
        self.thr = tk.Label(self.labelFrame3, text='Threads').grid(row=4, sticky="W")
        self.thread = tk.Spinbox(self.labelFrame3, from_=1, to=8, command=self.update_threads,
                              textvariable=default_threads)
        self.thread.grid(row=4, column=1, ipadx="9")

        # choose one summary stats
        self.labelFrame4 = tk.LabelFrame(self, text="Summary statistics")
        self.labelFrame4.grid(column=3, row=3, padx=20, pady=20, sticky="N")
        self.choice = tk.IntVar()
        self.choice.set(1)
        self.smtdbutton = tk.Radiobutton(self.labelFrame4, text='smtd', variable=self.choice,
                                      value=1, command=self.update_sstats).grid(row=0, column=0)
        self.umtdbutton = tk.Radiobutton(self.labelFrame4, text='umtd', variable=self.choice, value=2,
                                      command=self.update_sstats).grid(row=0, column=1)
        self.picsbutton = tk.Radiobutton(self.labelFrame4, text='pics', variable=self.choice, value=3,
                                      command=self.update_sstats).grid(row=0, column=2)

        # Output info
        self.labelFrame4 = tk.LabelFrame(self, text="Report of the progress")
        self.labelFrame4.grid(column=0, row=4, padx=20, pady=20, sticky="W", columnspan=5)
        # self.report = Text(master)
        # Set default values
        self.sstats_value = self.choice.get()
        self.threads_value = self.thread.get()
        self.iterations_value = self.iterations.get()
        self.particles_value = self.particles.get()

        self.gobutton = tk.Button(self.labelFrame4, text="Generate", command=self.generate_all,
                                  height=1,
                               width=10)
        self.gobutton.grid(row=0, column=0, sticky="E")

        self.outputtextbox = tk.Text(self.labelFrame4)
        self.outputtextbox.grid(column=0, row=2, sticky="W")

        # Open a file dialog

    def fileDialog(self):
        self.treedatadir = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=1)
        self.label.configure(text=self.treedatadir)
        print('Data file is in the directory: %s' % self.treedatadir)

    def savedir(self):
        self.save_dir = fd.askdirectory()
        self.label = tk.Label(self.labelFrame1, text="")
        self.label.grid(column=1, row=3)
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

if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
