import os
import sys
import threading
import time
import subprocess
import platform
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkCommonDialog
import ScrolledText
import tkMessageBox
import ttk

from math import ceil
from sklearn.manifold import MDS
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as col
from numpy import arange, sin, pi, random
import numpy as np
import networkx as nx

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor,_DistanceMatrix

APP_NAME = "CAFE: aCcelerated Alignment-FrEe sequence analysis"
GLOBAL_HASH_DIR = ""
GLOBAL_JELLYFISH_PATH = ""
GLOBAL_QUANTILE = 10

class ToolTip(object):
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, _cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 27
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))

        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def createToolTip( widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

class PreferencesWindow():    
    def __init__(self, view):
        self.parent = view
        self.openFile_icon = PhotoImage(file='image/openfile.gif')
        self.jellyfishVal = StringVar()
        self.hashVal = StringVar()
        self.quantVal = IntVar()
        self.create_prefereces_window()
        
    def create_prefereces_window(self):
        self.pref_window = Toplevel(self.parent)
        self.pref_window.title("Settings")
        self.create_prefereces_list()
        self.pref_window.transient(self.parent)
        
    def create_prefereces_list(self):
        global GLOBAL_HASH_DIR
        global GLOBAL_JELLYFISH_PATH
        global GLOBAL_QUANTILE
        
        self.jellyfishVal.set(GLOBAL_JELLYFISH_PATH)
        Label(self.pref_window, text="Jellyfish executable path:").grid(row=1, column=0, sticky=W, padx=5, pady=5)
        ttk.Button(self.pref_window, image=self.openFile_icon, command=self.on_open_file_jellyfish).grid(row=1, column=1, columnspan=2, sticky=E, padx=5, pady=5)
        self.jellyfish_Entry = Entry(self.pref_window, textvariable=self.jellyfishVal).grid(row=1, column=4, columnspan=5, pady=4)
        
        self.hashVal.set(GLOBAL_HASH_DIR)
        Label(self.pref_window, text="Saved hash directory:").grid(row=2, column=0, sticky=W, padx=5, pady=5)
        ttk.Button(self.pref_window, image=self.openFile_icon, command=self.on_open_file_hash).grid(row=2, column=1, columnspan=2, sticky=E, padx=5, pady=5)
        self.hash_Entry = Entry(self.pref_window, textvariable=self.hashVal).grid(row=2, column=4, columnspan=5, pady=4)
        
        self.quantVal.set(GLOBAL_QUANTILE)
        Label(self.pref_window, text="Quantile for Network:").grid(row=3, column=0, sticky=W, padx=5, pady=5)
        Spinbox(self.pref_window, from_=1, to=100, width=4, textvariable=self.quantVal, increment=1).grid(row=3, column=1, columnspan=2, pady=4)
        
        ttk.Button(self.pref_window, text="Save", command=self.on_save_button_clicked).grid(
            row=5, column=2, sticky=E, padx=5, pady=5)
        ttk.Button(self.pref_window, text="Cancel", command=self.on_cancel_button_clicked).grid(
            row=5, column=1, sticky=E, padx=5, pady=5)
    

    def on_open_file_jellyfish(self):
        jellyfish_path = tkFileDialog.askopenfilename(filetypes=[('All files', '*.*')])
        if jellyfish_path:
            self.jellyfishVal.set(jellyfish_path)
    
    def on_open_file_hash(self):
        hash_dir = tkFileDialog.askdirectory()
        if hash_dir:
            self.hashVal.set(hash_dir)
    
    def on_save_button_clicked(self):
        global GLOBAL_HASH_DIR
        global GLOBAL_JELLYFISH_PATH
        global GLOBAL_QUANTILE

        GLOBAL_HASH_DIR = self.hashVal.get()
        GLOBAL_JELLYFISH_PATH = self.jellyfishVal.get()
        GLOBAL_QUANTILE = self.quantVal.get()
        
        print(GLOBAL_HASH_DIR)
        print(GLOBAL_JELLYFISH_PATH)
        print(GLOBAL_QUANTILE)
        self.pref_window.destroy()
    
    def on_cancel_button_clicked(self):
        self.pref_window.destroy()
    
class GUIApp:
       
    def __init__(self, root):
        self.root = root
        self.distVal = StringVar()
        self.kVal = IntVar()
        self.mVal = IntVar()
        self.thresVal = IntVar()
        self.revComplVal = BooleanVar()
        self.consoleVal = StringVar()
        
        self.distVal.set('Ma')
        self.kVal.set(8)
        self.mVal.set(0)
        self.thresVal.set(0)
        self.revComplVal.set(False)
        self.logURL = 'log.txt'
        self.logURL_replicate = 'log1.txt'
        
        self.input_list = []
        self.label_list = []
        
        self.addFile_icon = PhotoImage(file='image/addFile.gif')
        self.addDir_icon = PhotoImage(file='image/addDir.gif')
        self.remove_icon = PhotoImage(file='image/remove.gif')
        self.clear_icon = PhotoImage(file='image/clear.gif')
        #self.save_icon = PhotoImage(file='image/save.gif')
        self.setting_icon = PhotoImage(file='image/setting.gif')
        
        self.root.protocol('WM_DELETE_WINDOW', self.exit_app)
        self.create_gui()

    def create_gui(self):
        self.root.title(APP_NAME)
        self.create_menu()
        self.create_top_bar()
        self.create_left_bar()
        self.create_right_frame_tab()
    
    def create_menu(self):
        menu_bar = Menu(self.root)
        file_menu = Menu(menu_bar, tearoff=0)
        file_menu.add_command(label='Add File', compound='left', image=self.addFile_icon, command=self.on_addFile_button_clicked)
        file_menu.add_command(label='Add Dir', compound='left', image=self.addDir_icon, command=self.on_addDir_button_clicked)
        file_menu.add_command(label='Remove Selected', compound='left', image=self.remove_icon, command=self.on_remove_button_clicked)
        file_menu.add_command(label='Remove All', compound='left', image=self.clear_icon, command=self.on_clear_button_clicked)
        #file_menu.add_command(label='Save', state='disabled', compound='left', image=self.save_icon, command=self.on_save_button_clicked)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.exit_app)
        menu_bar.add_cascade(label='File', menu=file_menu)
        
        about_menu = Menu(menu_bar, tearoff=0)
        about_menu.add_command(label='Setting', compound='left', image=self.setting_icon, command=self.on_setting_button_clicked)
        about_menu.add_separator()
        about_menu.add_command(label="About", command=self.about_app)
        about_menu.add_command(label="Help", command=self.help_app)
        menu_bar.add_cascade(label='About', menu=about_menu)
        self.root.config(menu=menu_bar)
    
    def create_top_bar(self):
        topbar_frame = Frame(self.root, height=25)
        topbar_frame.grid(row=0, columnspan=20, rowspan=10, pady=5)
        
        file_button = ttk.Button(topbar_frame, image=self.addFile_icon, command=self.on_addFile_button_clicked)
        file_button.grid(row=0, column=1)
        createToolTip(file_button, 'Add File')
        
        dir_button = ttk.Button(topbar_frame, image=self.addDir_icon, command=self.on_addDir_button_clicked)
        dir_button.grid(row=0, column=2)
        createToolTip(dir_button, 'Add Dir')
        
        remove_button = ttk.Button(topbar_frame, image=self.remove_icon, command=self.on_remove_button_clicked)
        remove_button.grid(row=0, column=3)
        createToolTip(remove_button, 'Remove Selected')
        
        clear_button = ttk.Button(topbar_frame, image=self.clear_icon, command=self.on_clear_button_clicked)
        clear_button.grid(row=0, column=4)
        createToolTip(clear_button, 'Remove all')
        
        ttk.Separator(topbar_frame, orient='vertical').grid(row=0, column=5, sticky="ns", padx=5)
        
        Label(topbar_frame, text='Distance:').grid(row=0, column=6)
        dist_combobox = ttk.Combobox(topbar_frame, width=10, textvariable=self.distVal)
        dist_combobox.grid(row=0, column=7)
        dist_combobox['values'] = ('Anderberg','Antidice','Canberra','Ch','Chisq','Cosine','CVtree','Dice','D2','D2star','D2shepp','Eu','Gower','Hamman','Hamming','Jaccard','JS','Kulczynski','Ma','Matching','Ochiai','Pearson','Phi','Russel','Sneath','Tanimoto','Yule')
        dist_combobox.set('Ma')
        dist_combobox.bind('<<ComboboxSelected>>', self.on_dist_changed)
        
        Label(topbar_frame, text='K:').grid(row=0, column=8, padx=3)
        self.kVal.set(8)
        self.kVal_spinbox = Spinbox(topbar_frame, from_=2, to=20, width=4, textvariable=self.kVal, increment=1, command=self.on_k_changed)
        self.kVal_spinbox.grid(row=0, column=9)
        
        Label(topbar_frame, text='Markov Order:').grid(row=0, column=10, padx=3)
        self.mVal.set(0)
        self.mVal_spinbox = Spinbox(topbar_frame, state='disabled', from_=-1, to=20, width=4, textvariable=self.mVal, increment=1)
        self.mVal_spinbox.grid(row=0, column=11)
        
        Label(topbar_frame, text='Threshold:').grid(row=0, column=12, padx=3)
        self.thresVal.set(0)
        Spinbox(topbar_frame, from_=0, to=100, width=4, textvariable=self.thresVal, increment=1).grid(row=0, column=13)
        
        self.revComplVal.set(False)
        ttk.Checkbutton(topbar_frame, text='Reverse Complmentary', variable=self.revComplVal).grid(row=0, column=14, padx=5)
        ttk.Separator(topbar_frame, orient='vertical').grid(row=0, column=15, sticky="ns", padx=5)
        
        self.run_button = ttk.Button(topbar_frame, text='Run',  command=self.on_run_button_clicked)
        self.run_button.grid(row=0, column=18, padx=5)
    
    def create_left_bar(self):
        Label(self.root, text='Input:').grid(row=11,column=1, sticky='w')
        
        file_frame = Frame(self.root)
        file_frame.grid(row=12, column=1, columnspan=4, rowspan = 10, sticky='wse')
        
        self.file_box = Listbox(file_frame, activestyle='none', cursor='hand2', selectmode=EXTENDED)
        self.file_box.pack(side=LEFT, fill=BOTH, expand=1)
        self.file_box.bind("<Button-3>", self.show_context_menu)
        
        file_scroll_bar = Scrollbar(file_frame)
        file_scroll_bar.pack(side=RIGHT, fill=BOTH)
        self.file_box.config(yscrollcommand=file_scroll_bar.set)
        file_scroll_bar.config(command=self.file_box.yview)
        
        self.context_menu = Menu(self.file_box, tearoff=0)
        self.context_menu.add_command( label="Delete", command=self.on_remove_selected_context_menu_clicked)
        
        Label(self.root, text='Console:').grid(row=22, column=1, sticky='w')
        
        console_frame = Frame(self.root)
        console_frame.grid(row=23, column=1, columnspan=4, rowspan = 10, sticky='wse')
        
        self.console_box = Text(console_frame, state=DISABLED, width=40)
        self.console_box.pack(side=LEFT, fill=BOTH, expand=1)
        
        console_scroll_bar = Scrollbar(console_frame)
        console_scroll_bar.pack(side=RIGHT, fill=BOTH)
        self.console_box.config(yscrollcommand=console_scroll_bar.set)
        console_scroll_bar.config(command=self.console_box.yview)
    
    def create_right_frame_tab(self):
        tabControl = ttk.Notebook(self.root)
        
        self.tab4 = ttk.Frame(tabControl)
        tabControl.add(self.tab4, text='Dendrogram')
        self.tab3 = ttk.Frame(tabControl)
        tabControl.add(self.tab3, text='Principal Coordinates Analysis')
        self.tab2 = ttk.Frame(tabControl)
        tabControl.add(self.tab2, text='Heatmap')
        self.tab1 = ttk.Frame(tabControl) 
        tabControl.add(self.tab1, text='Network')
        
        tabControl.grid(row=14, column=6, columnspan=15, rowspan = 10, padx=2, sticky='nw')
        
        

        
    def show_context_menu(self, event):
        self.context_menu.tk_popup(event.x_root, event.y_root)
    
    def on_addFile_button_clicked(self):
        input_file = tkFileDialog.askopenfilename(filetypes=[('All supported', '.fasta .fa .fna'), ('.fasta files', '.fasta'), ('.fa files', '.fa'), ('.fna files', '.fna')])
        if input_file:
            if not input_file in self.input_list:   
                input_file_path, input_file_name = os.path.split(input_file)
                self.label_list.append(input_file_name)
                self.input_list.append(input_file)
                self.file_box.insert(END, input_file_name)
            
    def on_remove_button_clicked(self):
        try:
            selected_indexes = self.file_box.curselection()
            for index in reversed(selected_indexes):
                self.file_box.delete(index)
                del self.input_list[index]
                del self.label_list[index]
        except IndexError:
            pass
    
    def on_remove_selected_context_menu_clicked(self):
        self.on_remove_button_clicked()
    
    def on_addDir_button_clicked(self):
        directory_path = tkFileDialog.askdirectory()
        if not directory_path:
            return
        input_files_in_directory = []
        for (dirpath, dirnames, filenames) in os.walk(directory_path):
            for input_file in filenames:
                if input_file.endswith(".fasta") or input_file.endswith(".fa") or input_file.endswith(".fna"):
                    input_files_in_directory.append(dirpath + "/" + input_file)
        for input_file in input_files_in_directory:
            if not input_file in self.input_list: 
                input_file_path, input_file_name = os.path.split(input_file)
                self.label_list.append(input_file_name)
                self.input_list.append(input_file)
                self.file_box.insert(END, input_file_name)   
    
    def on_clear_button_clicked(self):
        del self.input_list[:]
        self.input_list = []
        del self.label_list[:]
        self.label_list = []
        self.file_box.delete(0, END)
    
#    def on_save_button_clicked(self):
#        pass
    
    def on_setting_button_clicked(self):
        PreferencesWindow(self.root)
    
    def display_console(self):
        self.check_console = True
        while self.check_console:
            self.console_box.config(state='normal')
            self.console_box.delete(1.0, 'end')
            
            if os.path.isfile(self.logURL):
                os.system("cp " + self.logURL + " " + self.logURL_replicate)  
                
                log_lines = [line.rstrip('\n') for line in open(self.logURL_replicate)]
                for log_line in log_lines:
                    print(log_line)
                    message = log_line.encode('utf-8')
                    self.console_box.insert('end', message.decode('utf-8') + '\n')
                    self.console_box.yview(END)
                    
                if log_lines[-1] == "Done":
                    self.check_console = False
                    self.run_button.config(state='normal')

            self.console_box.config(state='disabled')
            time.sleep(1)
    
    def callViz(self):
        vizURL = "result."+self.distVal.get()+".phylip"
        print("vizURL: "+vizURL)
        if not os.path.isfile(vizURL):
            tkMessageBox.showerror(APP_NAME,vizURL + "  Not Exist!!")
            return
        
        for widget in self.tab4.winfo_children():
            widget.destroy()
        for widget in self.tab3.winfo_children():
            widget.destroy()
        for widget in self.tab2.winfo_children():
            widget.destroy()
        for widget in self.tab1.winfo_children():
            widget.destroy()
        
        #phylogenetic tree
        canvas_frame_tab4 = Frame(self.tab4)
        canvas_frame_tab4.grid(row=0, column=0, columnspan=12, rowspan = 8, sticky='nw')
        canvas_tab4 = Canvas(canvas_frame_tab4, background="white", width=500, height=500)
        
        x_scroll_tab4 = Scrollbar(canvas_frame_tab4, orient="horizontal")
        x_scroll_tab4.pack(side="bottom", fill="x")
        x_scroll_tab4.config(command=canvas_tab4.xview)
        y_scroll_tab4 = Scrollbar(canvas_frame_tab4, orient="vertical")
        y_scroll_tab4.pack(side="right", fill="y")
        y_scroll_tab4.config(command=canvas_tab4.yview)
        
        fig_tab4, ax_tab4 = plt.subplots()
        #tree = Phylo.read(StringIO('(A,(B,C),(D,E));'), 'newick')
        
        firstCol = np.loadtxt(vizURL,usecols=range(1),dtype='str')
        genomeSize = int(firstCol[0])
        genomeDist = np.loadtxt(vizURL,usecols=range(1,len(firstCol)),skiprows=1)
        genome_min = np.amin(genomeDist)
        genomeDist = genomeDist - genome_min;
        genome_min = np.amin(genomeDist)
        genome_max = np.amax(genomeDist)
        genomeDist_normed = genomeDist / genome_max
        genome_cmap = LinearSegmentedColormap.from_list('mycmap', [(0, 'blue'),(1, 'white')])
        
        genomeDist = genomeDist.tolist()
        genomeList = firstCol[1:].tolist()
        
        genomeTriDist = []
        for row in range(genomeSize):
            genomeTriDist.append(genomeDist[row][0:row+1])
        dm = _DistanceMatrix(genomeList, matrix=genomeTriDist)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        
        tree.ladderize()
        Phylo.draw(tree,do_show=False,axes=ax_tab4)
        
        canvas_tab4 = FigureCanvasTkAgg(fig_tab4, master=canvas_frame_tab4)
        canvas_tab4.get_tk_widget().config(xscrollcommand=x_scroll_tab4.set, yscrollcommand=y_scroll_tab4.set)
        canvas_tab4.get_tk_widget().pack(side=RIGHT, fill=BOTH, expand=YES)
        
        #network
        canvas_frame_tab1 = Frame(self.tab1)
        canvas_frame_tab1.grid(row=0, column=0, columnspan=12, rowspan = 8, sticky='nw')
        canvas_tab1 = Canvas(canvas_frame_tab1, background="white", width=500, height=500)
        
        x_scroll_tab1 = Scrollbar(canvas_frame_tab1, orient="horizontal")
        x_scroll_tab1.pack(side="bottom", fill="x")
        x_scroll_tab1.config(command=canvas_tab1.xview)
        y_scroll_tab1 = Scrollbar(canvas_frame_tab1, orient="vertical")
        y_scroll_tab1.pack(side="right", fill="y")
        y_scroll_tab1.config(command=canvas_tab1.yview)
        
        fig_tab1, ax_tab1 = plt.subplots()
        edgeVal = []
        for fr in range(1,genomeSize):
            for to in range(fr+1,genomeSize):
                edgeVal.append(genomeDist[fr][to])
        edgeVal.sort();
        global GLOBAL_QUANTILE
        idx = int(ceil(GLOBAL_QUANTILE*0.01*len(edgeVal)))
        if idx >= len(edgeVal)-1:
            idx = len(edgeVal)-1
        thres = edgeVal[ idx ]
        
        G = nx.Graph()
        for row in range(genomeSize):
            G.add_node(row)
        for fr in range(1,genomeSize):
            for to in range(fr+1,genomeSize):
                if genomeDist[fr][to] <= thres:
                    G.add_edge(fr, to, weight=genomeDist[fr][to])      
                              
        labelsNet = {}
        for node in G.nodes():
            labelsNet[node] = genomeList[node]        
        
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos) 
        nx.draw_networkx_labels(G, pos,labelsNet, font_size=8)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        plt.axis('off')
        
        canvas_tab1 = FigureCanvasTkAgg(fig_tab1, master=canvas_frame_tab1)
        canvas_tab1.get_tk_widget().config(xscrollcommand=x_scroll_tab1.set, yscrollcommand=y_scroll_tab1.set)
        canvas_tab1.get_tk_widget().pack(side=RIGHT, fill=BOTH, expand=YES)
        
        #heatmap
        canvas_frame_tab2 = Frame(self.tab2)
        canvas_frame_tab2.grid(row=0, column=0, columnspan=12, rowspan = 8, sticky='nw')
        canvas_tab2 = Canvas(canvas_frame_tab2, background="white", width=500, height=500)
        
        x_scroll_tab2 = Scrollbar(canvas_frame_tab2, orient="horizontal")
        x_scroll_tab2.pack(side="bottom", fill="x")
        x_scroll_tab2.config(command=canvas_tab2.xview)
        y_scroll_tab2 = Scrollbar(canvas_frame_tab2, orient="vertical")
        y_scroll_tab2.pack(side="right", fill="y")
        y_scroll_tab2.config(command=canvas_tab2.yview)
        
        fig_tab2, ax_tab2 = plt.subplots()
        im = ax_tab2.pcolor(genomeDist_normed, cmap=genome_cmap, vmin=genome_min, vmax=genome_max, edgecolors='white')
        cbar = fig_tab2.colorbar(im)
        ax_tab2.set_frame_on(False)
        ax_tab2.set_aspect('equal')
        ax_tab2.invert_yaxis()
        ax_tab2.set_xticks(np.arange(genomeDist_normed.shape[1]) + 0.5, minor=False)
        ax_tab2.set_xticklabels(genomeList,rotation=60,fontsize=8)
        ax_tab2.set_yticklabels([])
        
        canvas_tab2 = FigureCanvasTkAgg(fig_tab2, master=canvas_frame_tab2)
        canvas_tab2.get_tk_widget().config(xscrollcommand=x_scroll_tab2.set, yscrollcommand=y_scroll_tab2.set)
        canvas_tab2.get_tk_widget().pack(side=RIGHT, fill=BOTH, expand=YES)
        
        #PCOA
        canvas_frame_tab3 = Frame(self.tab3)
        canvas_frame_tab3.grid(row=0, column=0, columnspan=12, rowspan = 8, sticky='nw')
        canvas_tab3 = Canvas(canvas_frame_tab3, background="white", width=500, height=500)
        
        x_scroll_tab3 = Scrollbar(canvas_frame_tab3, orient="horizontal")
        x_scroll_tab3.pack(side="bottom", fill="x")
        x_scroll_tab3.config(command=canvas_tab3.xview)
        y_scroll_tab3 = Scrollbar(canvas_frame_tab3, orient="vertical")
        y_scroll_tab3.pack(side="right", fill="y")
        y_scroll_tab3.config(command=canvas_tab3.yview)
        
        fig_tab3, ax_tab3 = plt.subplots()

        mds = MDS(n_components=2, dissimilarity='precomputed')
        mdsDim = mds.fit_transform(genomeDist)
        ax_tab3.scatter(mdsDim[:,0], mdsDim[:,1], s = 100.)
        for i, txt in enumerate(genomeList):
            ax_tab3.annotate(txt, (mdsDim[i,0], mdsDim[i,1]),fontsize=8)
        plt.grid(True)
        
        canvas_tab3 = FigureCanvasTkAgg(fig_tab3, master=canvas_frame_tab3)
        canvas_tab3.get_tk_widget().config(xscrollcommand=x_scroll_tab3.set, yscrollcommand=y_scroll_tab3.set)
        canvas_tab3.get_tk_widget().pack(side=RIGHT, fill=BOTH, expand=YES)
    
    def callCAFE(self):
        optSys = platform.system()
        if optSys == 'Linux' :
            exePath = "./cafe_linux"
            os.system("chmod 777 " + exePath)
        elif optSys == 'Darwin' :
            exePath = "./cafe_mac"
            os.system("chmod 777 " + exePath)
        elif optSys == 'Windows' :
            exePath = "cafe_win.exe"
            import ctypes
            SEM_NOGPFAULTERRORBOX = 0x0002 
            ctypes.windll.kernel32.SetErrorMode(SEM_NOGPFAULTERRORBOX);
            subprocess_flags = 0x8000000

        global GLOBAL_HASH_DIR
        global GLOBAL_JELLYFISH_PATH

        paramArr = [exePath]
        paramArr.append("-I")
        paramArr.append(",".join(self.input_list))
        if GLOBAL_HASH_DIR:
            paramArr.append("-S")
            paramArr.append(GLOBAL_HASH_DIR)
        if GLOBAL_JELLYFISH_PATH:
            paramArr.append("-J")
            paramArr.append(GLOBAL_JELLYFISH_PATH)
        if self.thresVal.get():
            paramArr.append("-L")
            paramArr.append(str(self.thresVal.get()))
        
        paramArr.append("-K")
        paramArr.append(str(self.kVal.get()))
        paramArr.append("-M")
        paramArr.append(str(self.mVal.get()))
        paramArr.append("-D")
        paramArr.append(self.distVal.get())
        paramArr.append("-O")
        paramArr.append("result")
        if self.revComplVal.get():
            paramArr.append("-R")
        paramArr.append(">")
        paramArr.append(self.logURL)
        
        #absFileName = os.path.abspath('cafe.exe')
        dirName = os.path.dirname(os.path.abspath(__file__))
        p = subprocess.Popen(paramArr, cwd=dirName, shell=True)
        p.communicate()
        
        self.callViz()
    
    def on_run_button_clicked(self,event=None):
        optSys = platform.system()
        if optSys == 'Linux' :
            exePath = "cafe_linux"
        elif optSys == 'Darwin' :
            exePath = "cafe_mac"
        elif optSys == 'Windows' :
            exePath = "cafe_win.exe"
        
        if not len(self.input_list) or len(self.input_list)==1:
            tkMessageBox.showerror(APP_NAME,'Program cannot proceed when the input size <= 1')
            return
        
        os.system("rm " + self.logURL)
        os.system("rm " + self.logURL_replicate)  
        
        self.checkConsoleT = threading.Thread(target=self.display_console)
        self.checkConsoleT.setDaemon(True)    
        self.checkConsoleT.start()
        
        self.callT = threading.Thread(target=self.callCAFE)
        self.callT.setDaemon(True)    
        self.callT.start()
        
        self.run_button.config(state='disabled')
        
    def exit_app(self):
        self.check_console = False
        if tkMessageBox.askokcancel("Quit", "Do you really want to quit?"):
            self.root.destroy()
    
    def about_app(self):
        tkMessageBox.showinfo(APP_NAME,'Copyright (C) 2016 University of Southern California\nAuthors: \nYang Lu: ylu465@usc.edu\nProf. Fengzhu Sun: fsun@usc.edu')
    
    def help_app(self):
        tkMessageBox.showinfo(APP_NAME,'Manual is available at https://github.com/younglululu/CAFE')
    
    def on_dist_changed(self,event=None):
        if self.distVal.get() == 'Ch' or self.distVal.get() == 'Eu' or self.distVal.get() == 'Ma' or self.distVal.get() == 'D2' \
        or self.distVal.get() == 'Chisq' or self.distVal.get() == 'Cosine' or self.distVal.get() == 'Canberra' or self.distVal.get() == 'Pearson' \
        or self.distVal.get() == 'Hamming' or self.distVal.get() == 'Anderberg' or self.distVal.get() == 'Antidice' or self.distVal.get() == 'Dice' \
        or self.distVal.get() == 'Gower' or self.distVal.get() == 'Hamman' or self.distVal.get() == 'Jaccard' or self.distVal.get() == 'Kulczynski' \
        or self.distVal.get() == 'Matching' or self.distVal.get() == 'Ochiai' or self.distVal.get() == 'Pearson' or self.distVal.get() == 'Phi' \
        or self.distVal.get() == 'Russel' or self.distVal.get() == 'Sneath' or self.distVal.get() == 'Tanimoto' or self.distVal.get() == 'Yule':
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='disabled')
        elif self.distVal.get() == 'D2star' or self.distVal.get() == 'D2shepp':
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='normal')
        elif self.distVal.get() == 'CVtree':
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='disabled')
        elif self.distVal.get() == 'JS':
            self.kVal_spinbox.config(state='disabled')
            self.mVal_spinbox.config(state='normal')
        else:
            print("No suitable distance!")
    
    def on_k_changed(self,event=None):
        max_mVal = int(self.kVal.get())-1
        self.mVal.set(min(0, max_mVal))
        self.mVal_spinbox.config(to=max_mVal)

if __name__ == '__main__':
    root = Tk()
    app = GUIApp(root)
    root.geometry('{}x{}'.format(1000, 630))
    root.resizable(False, False)
    root.iconbitmap('image/logo.ico')
    root.mainloop()