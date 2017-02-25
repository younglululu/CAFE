import os
import sys
import threading
import time
import subprocess
import platform
import ctypes

#from Tkinter import *
from mtTkinter import *
import Tkconstants
import tkFileDialog
import tkMessageBox
import ttk
from ttk import *

from math import ceil
from sklearn.manifold import MDS

if platform.system() == 'Darwin':
    import matplotlib
    matplotlib.use("TkAgg")

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
        self.optSys = platform.system()
        self.currentTab = 4
        self.logURL = 'log.txt'
        self.logURL_replicate = 'log1.txt'
        
        self.input_list = []
        self.label_list = []
        self.log_list = []
        
        self.addFile_icon = PhotoImage(file='image/addFile.gif')
        self.addDir_icon = PhotoImage(file='image/addDir.gif')
        self.remove_icon = PhotoImage(file='image/remove.gif')
        self.clear_icon = PhotoImage(file='image/clear.gif')
        self.load_icon = PhotoImage(file='image/load.gif')
        self.setting_icon = PhotoImage(file='image/setting.gif')
        self.zoomin_icon = PhotoImage(file='image/zoomin.gif')
        self.zoomout_icon = PhotoImage(file='image/zoomout.gif')
        self.save_icon = PhotoImage(file='image/save.gif')
        
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
        file_menu.add_command(label='Setting', compound='left', image=self.setting_icon, command=self.on_setting_button_clicked)
        file_menu.add_separator()
        file_menu.add_command(label="About", command=self.about_app)
        file_menu.add_command(label="Help", command=self.help_app)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.exit_app)
        menu_bar.add_cascade(label='File', menu=file_menu)
        
        self.root.config(menu=menu_bar)
    
    def create_top_bar(self):
        topbar_frame = Frame(self.root, height=25)
        topbar_frame.grid(row=0, columnspan=20, rowspan=10, pady=5, sticky=Tkconstants.NSEW)
        
        load_button = ttk.Button(topbar_frame, image=self.load_icon, command=self.on_load_button_clicked)
        load_button.grid(row=0, column=0)
        createToolTip(load_button, 'Load Existing Results in Phylip format')
        
        file_button = ttk.Button(topbar_frame, image=self.addFile_icon, command=self.on_addFile_button_clicked)
        file_button.grid(row=0, column=1)
        createToolTip(file_button, 'Add one genome sequence to the list')
        
        dir_button = ttk.Button(topbar_frame, image=self.addDir_icon, command=self.on_addDir_button_clicked)
        dir_button.grid(row=0, column=2)
        createToolTip(dir_button, 'Add all genome sequences from directory to the list')
        
        remove_button = ttk.Button(topbar_frame, image=self.remove_icon, command=self.on_remove_button_clicked)
        remove_button.grid(row=0, column=3)
        createToolTip(remove_button, 'Remove Selected genome sequences in the list')
        
        clear_button = ttk.Button(topbar_frame, image=self.clear_icon, command=self.on_clear_button_clicked)
        clear_button.grid(row=0, column=4)
        createToolTip(clear_button, 'Remove all genome sequences in the list')
        
        ttk.Separator(topbar_frame, orient='vertical').grid(row=0, column=5, sticky="ns", padx=3)
        
        Label(topbar_frame, text='Distance:').grid(row=0, column=6)
        dist_combobox = ttk.Combobox(topbar_frame, width=10, textvariable=self.distVal)
        dist_combobox.grid(row=0, column=7)
        dist_combobox['values'] = ('Anderberg','Antidice','Canberra','Ch','Chisq','Cosine','Co-phylog','CVtree','Dice','D2','D2star','D2shepp','Eu','FFP','Gower','Hamman','Hamming','Jaccard','JS','Kulczynski','Ma','Matching','Ochiai','Pearson','Phi','Russel','Sneath','Tanimoto','Yule')
        dist_combobox.set('Ma')
        dist_combobox.bind('<<ComboboxSelected>>', self.on_dist_changed)
        
        Label(topbar_frame, text='K:').grid(row=0, column=8, padx=3)
        self.kVal.set(8)
        self.kVal_spinbox = Spinbox(topbar_frame, from_=2, to=20, width=4, textvariable=self.kVal, increment=1, command=self.on_k_changed)
        self.kVal_spinbox.grid(row=0, column=9)
        
        Label(topbar_frame, text='Markov Order:').grid(row=0, column=10, padx=3)
        self.mVal.set(0)
        self.mVal_spinbox = Spinbox(topbar_frame, state='disabled', from_=0, to=0, width=4, textvariable=self.mVal, increment=1)
        self.mVal_spinbox.grid(row=0, column=11)
        
        Label(topbar_frame, text='Cutoff:').grid(row=0, column=12, padx=3)
        self.thresVal.set(0)
        Spinbox(topbar_frame, from_=0, to=100, width=4, textvariable=self.thresVal, increment=1).grid(row=0, column=13)
        
        self.revComplVal.set(False)
        ttk.Checkbutton(topbar_frame, text='Rev.Complement', variable=self.revComplVal).grid(row=0, column=14, padx=3)
        ttk.Separator(topbar_frame, orient='vertical').grid(row=0, column=15, sticky="ns", padx=3)
        
        self.run_button = ttk.Button(topbar_frame, text='Run', state='disabled', command=self.on_run_button_clicked)
        self.run_button.grid(row=0, column=16, padx=3)
        
        ttk.Separator(topbar_frame, orient='vertical').grid(row=0, column=17, sticky="ns", padx=3)
        
        self.zoomin_button = ttk.Button(topbar_frame, image=self.zoomin_icon, state='disabled', command=lambda:self.on_zoom_clicked(1))
        self.zoomin_button.grid(row=0, column=18)
        createToolTip(self.zoomin_button, 'Zoom in the current figure')
        
        self.zoomout_button = ttk.Button(topbar_frame, image=self.zoomout_icon, state='disabled', command=lambda:self.on_zoom_clicked(-1))
        self.zoomout_button.grid(row=0, column=19)
        createToolTip(self.zoomout_button, 'Zoom out the current figure')
        
        self.save_button = ttk.Button(topbar_frame, image=self.save_icon, state='disabled', command=self.on_saveFig_clicked)
        self.save_button.grid(row=0, column=20)
        createToolTip(self.save_button, 'Save the current figure')
    
    def create_left_bar(self):
        Label(self.root, text='Input:').grid(row=11,column=0, sticky='w')
        
        file_frame = Frame(self.root)
        file_frame.grid(row=12, column=0, columnspan=4, rowspan = 10, sticky=Tkconstants.NSEW)
        
        self.file_box = Listbox(file_frame, activestyle='none', cursor='hand2', selectmode=EXTENDED)
        self.file_box.pack(side=LEFT, fill=BOTH, expand=1)
        
        if self.optSys == 'Darwin':
            self.file_box.bind("<Button-2>", self.show_context_menu)
        else:
            self.file_box.bind("<Button-3>", self.show_context_menu)
        
        file_scroll_bar = Scrollbar(file_frame)
        file_scroll_bar.pack(side=RIGHT, fill=BOTH)
        self.file_box.config(yscrollcommand=file_scroll_bar.set)
        file_scroll_bar.config(command=self.file_box.yview)
        
        self.context_menu = Menu(self.file_box, tearoff=0)
        self.context_menu.add_command( label="Delete", command=self.on_remove_selected_context_menu_clicked)
        
        Label(self.root, text='Console:').grid(row=22, column=0, sticky='w')
        
        console_frame = Frame(self.root)
        console_frame.grid(row=23, column=0, columnspan=4, rowspan = 10, sticky=Tkconstants.NSEW)
        
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
        
        tabControl.grid(row=14, column=5, columnspan=15, rowspan = 10, padx=2, sticky='nw')
        
    def show_context_menu(self, event):
        self.context_menu.tk_popup(event.x_root, event.y_root)
    
    def on_load_button_clicked(self):
        vizURL = tkFileDialog.askopenfilename(filetypes=[('All supported', '.phylip'), ('.phylip files', '.phylip')])
        if not vizURL:
            return
        self.callViz(vizURL)
    
    def on_addFile_button_clicked(self):
        input_file = tkFileDialog.askopenfilename(filetypes=[('All supported', '.fasta .fa .fna'), ('.fasta files', '.fasta'), ('.fa files', '.fa'), ('.fna files', '.fna')])
        if input_file:
            if not input_file in self.input_list:   
                input_file_path, input_file_name = os.path.split(input_file)
                self.label_list.append(input_file_name)
                self.input_list.append(input_file)
                self.file_box.insert(END, input_file_name)    
        if len(self.input_list) > 1:
            self.run_button.config(state='normal')
            
    def on_remove_button_clicked(self):
        try:
            selected_indexes = self.file_box.curselection()
            for index in reversed(selected_indexes):
                self.file_box.delete(index)
                del self.input_list[index]
                del self.label_list[index]
        except IndexError:
            pass
        if len(self.input_list) <= 1:
            self.run_button.config(state='disabled')
    
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
        if len(self.input_list) > 1:
            self.run_button.config(state='normal')
                
    def on_clear_button_clicked(self):
        del self.input_list[:]
        self.input_list = []
        del self.label_list[:]
        self.label_list = []
        self.file_box.delete(0, END)
        self.run_button.config(state='disabled')
    
    def on_setting_button_clicked(self):
        PreferencesWindow(self.root)
    
    def display_console(self):
        self.check_console = True
        while self.check_console:
            self.console_box.config(state='normal')
            #self.console_box.delete(1.0, 'end')
            
            if os.path.isfile(self.logURL):
                if self.optSys == 'Windows' :
                    os.system("copy " + self.logURL + " " + self.logURL_replicate)
                else :
                    os.system("cp " + self.logURL + " " + self.logURL_replicate)  
                
                log_lines = [line.rstrip('\n') for line in open(self.logURL_replicate)]
                for lineCnt in range(len(log_lines)):
                    if lineCnt < len(self.log_list):
                        continue
                    message = log_lines[lineCnt].encode('utf-8')
                    self.console_box.insert('end', message.decode('utf-8') + '\n')
                    self.console_box.yview(END)
                    self.log_list.append(message);        

                if len(log_lines) > 1:
                    if log_lines[-1] == "Done":
                        self.check_console = False
                        self.run_button.config(state='normal')

            self.console_box.config(state='disabled')
            time.sleep(1)
            
        vizURL = "result."+self.distVal.get()+".phylip"
        self.callViz(vizURL)
    
    def _phyloLabel_callback(self,clade):
        if clade.name.startswith('Inner'): return None
        return clade
    
    def _bound_to_mousewheel(self, event):
        if event.widget == self.canvas_tab4_container:
            self.currentTab = 4
            self.canvas_tab4.get_tk_widget().bind_all("<MouseWheel>", self._plot_on_mousewheel_tab4)
            if self.optSys == 'Darwin':
                self.canvas_tab4.get_tk_widget().bind_all("<Button-2>", self.show_context_menu_tab)
            else:
                self.canvas_tab4.get_tk_widget().bind_all("<Button-3>", self.show_context_menu_tab)
        elif event.widget == self.canvas_tab3_container:
            self.currentTab = 3
            self.canvas_tab3.get_tk_widget().bind_all("<MouseWheel>", self._plot_on_mousewheel_tab3)
            if self.optSys == 'Darwin':
                self.canvas_tab3.get_tk_widget().bind_all("<Button-2>", self.show_context_menu_tab)
            else:
                self.canvas_tab3.get_tk_widget().bind_all("<Button-3>", self.show_context_menu_tab)
        elif event.widget == self.canvas_tab2_container:
            self.currentTab = 2
            self.canvas_tab2.get_tk_widget().bind_all("<MouseWheel>", self._plot_on_mousewheel_tab2)
            if self.optSys == 'Darwin':
                self.canvas_tab2.get_tk_widget().bind_all("<Button-2>", self.show_context_menu_tab)
            else:
                self.canvas_tab2.get_tk_widget().bind_all("<Button-3>", self.show_context_menu_tab)
        elif event.widget == self.canvas_tab1_container:
            self.currentTab = 1
            self.canvas_tab1.get_tk_widget().bind_all("<MouseWheel>", self._plot_on_mousewheel_tab1)
            if self.optSys == 'Darwin':
                self.canvas_tab1.get_tk_widget().bind_all("<Button-2>", self.show_context_menu_tab)
            else:
                self.canvas_tab1.get_tk_widget().bind_all("<Button-3>", self.show_context_menu_tab)
        else:
            pass 

    def _unbound_to_mousewheel(self, event):
        if event.widget == self.canvas_tab4_container:
            #self.canvas_tab4.get_tk_widget().unbind_all("<MouseWheel>")
            self.canvas_tab4.get_tk_widget().unbind_all("<Button-2>")
            self.canvas_tab4.get_tk_widget().unbind_all("<Button-3>")
        elif event.widget == self.canvas_tab3_container:
            #self.canvas_tab3.get_tk_widget().unbind_all("<MouseWheel>")
            self.canvas_tab3.get_tk_widget().unbind_all("<Button-2>")
            self.canvas_tab3.get_tk_widget().unbind_all("<Button-3>")
        elif event.widget == self.canvas_tab2_container:
            #self.canvas_tab2.get_tk_widget().unbind_all("<MouseWheel>")
            self.canvas_tab2.get_tk_widget().unbind_all("<Button-2>")
            self.canvas_tab2.get_tk_widget().unbind_all("<Button-3>")
        elif event.widget == self.canvas_tab1_container:
            #self.canvas_tab1.get_tk_widget().unbind_all("<MouseWheel>")
            self.canvas_tab1.get_tk_widget().unbind_all("<Button-2>")
            self.canvas_tab1.get_tk_widget().unbind_all("<Button-3>")
        else:
            pass 
    
    def show_context_menu_tab(self, event):
        if event.widget == self.canvas_tab4.get_tk_widget():
            self.context_menu_tab4.tk_popup(event.x_root, event.y_root)
        elif event.widget == self.canvas_tab3.get_tk_widget():
            self.context_menu_tab3.tk_popup(event.x_root, event.y_root)
        elif event.widget == self.canvas_tab2.get_tk_widget():
            self.context_menu_tab2.tk_popup(event.x_root, event.y_root)
        elif event.widget == self.canvas_tab1.get_tk_widget():
            self.context_menu_tab1.tk_popup(event.x_root, event.y_root)
        else:
            pass 
        
    def on_save_context_menu_clicked(self, tabID):
        filename = tkFileDialog.asksaveasfilename(filetypes=[('All supported', '.jpg .png'), ('.jpg files', '.jpg'), ('.png files', '.png')])
        
        if tabID == 4:
            self.fig_tab4.savefig(filename, bbox_inches='tight')
        elif tabID == 3:
            self.fig_tab3.savefig(filename, bbox_inches='tight')
        elif tabID == 2:
            self.fig_tab2.savefig(filename, bbox_inches='tight')
        elif tabID == 1:
            self.fig_tab1.savefig(filename, bbox_inches='tight')
        else:
            pass 
    
    def on_saveFig_clicked(self):
        if self.currentTab  == 4:
            self.on_save_context_menu_clicked(4)
        elif self.currentTab  == 3:
            self.on_save_context_menu_clicked(3)
        elif self.currentTab  == 2:
            self.on_save_context_menu_clicked(2)
        elif self.currentTab  == 1:
            self.on_save_context_menu_clicked(1)
        else:
            pass 
    
    def on_zoom_clicked(self, delta):
        factor = 1
        if delta > 0:
            factor = factor * 1.1;
        else:
            factor = factor / 1.1
        
        if self.currentTab  == 4:
            self.dealWithMouseWheel(factor, self.fig_tab4, self.canvas_tab4.get_tk_widget(), self.canvas_tab4_container, self.cwid_tab4)
        elif self.currentTab  == 3:
            self.dealWithMouseWheel(factor, self.fig_tab3, self.canvas_tab3.get_tk_widget(), self.canvas_tab3_container, self.cwid_tab3)
        elif self.currentTab  == 2:
            self.dealWithMouseWheel(factor, self.fig_tab2, self.canvas_tab2.get_tk_widget(), self.canvas_tab2_container, self.cwid_tab2)
        elif self.currentTab  == 1:
            self.dealWithMouseWheel(factor, self.fig_tab1, self.canvas_tab1.get_tk_widget(), self.canvas_tab1_container, self.cwid_tab1)
        else:
            pass 
      
    def _plot_on_mousewheel_tab4(self,event):
        factor = 1
        if event.delta/120 > 0:
            factor = factor * 1.1;
        else:
            factor = factor / 1.1
        self.dealWithMouseWheel(factor, self.fig_tab4, self.canvas_tab4.get_tk_widget(), self.canvas_tab4_container, self.cwid_tab4)

    def _plot_on_mousewheel_tab3(self,event):
        factor = 1
        if event.delta/120 > 0:
            factor = factor * 1.1;
        else:
            factor = factor / 1.1
        self.dealWithMouseWheel(factor, self.fig_tab3, self.canvas_tab3.get_tk_widget(), self.canvas_tab3_container, self.cwid_tab3)
        
    def _plot_on_mousewheel_tab2(self,event):
        factor = 1
        if event.delta/120 > 0:
            factor = factor * 1.1;
        else:
            factor = factor / 1.1
        self.dealWithMouseWheel(factor, self.fig_tab2, self.canvas_tab2.get_tk_widget(), self.canvas_tab2_container, self.cwid_tab2)
        
    def _plot_on_mousewheel_tab1(self,event):
        factor = 1
        if event.delta/120 > 0:
            factor = factor * 1.1;
        else:
            factor = factor / 1.1
        self.dealWithMouseWheel(factor, self.fig_tab1, self.canvas_tab1.get_tk_widget(), self.canvas_tab1_container, self.cwid_tab1)
    
    def dealWithMouseWheel(self, factor, figure, canvas, canvasContainer, window):
        oldSize = figure.get_size_inches()
        figure.set_size_inches([factor * s for s in oldSize])
        wi,hi = [i*figure.dpi for i in figure.get_size_inches()]
        canvas.config(width=wi, height=hi)
        canvasContainer.itemconfigure(window, width=wi, height=hi)
        canvasContainer.config(scrollregion=canvasContainer.bbox(Tkconstants.ALL),width=800,height=600)
        figure.canvas.draw()
    
    def callViz(self, vizURL):
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
        canvas_frame_tab4.grid(row=1, column=1, sticky=Tkconstants.NS)
        canvas_frame_tab4.rowconfigure(1, weight=1)
        canvas_frame_tab4.columnconfigure(1, weight=1)
        
        self.fig_tab4, ax_tab4 = plt.subplots()
        ax_tab4.autoscale(enable=True)
        
        self.canvas_tab4_container = Canvas(canvas_frame_tab4)
        self.canvas_tab4_container.grid(row=1, column=1, sticky=Tkconstants.NSEW)
        x_scroll_tab4 = Scrollbar(canvas_frame_tab4, orient=HORIZONTAL)
        y_scroll_tab4 = Scrollbar(canvas_frame_tab4, orient=VERTICAL)
        x_scroll_tab4.grid(row=2, column=1, sticky=Tkconstants.EW)
        y_scroll_tab4.grid(row=1,column=2, sticky=Tkconstants.NS)
        self.canvas_tab4_container.config(xscrollcommand=x_scroll_tab4.set)
        x_scroll_tab4.config(command=self.canvas_tab4_container.xview)
        self.canvas_tab4_container.config(yscrollcommand=y_scroll_tab4.set)
        y_scroll_tab4.config(command=self.canvas_tab4_container.yview)
        
        self.context_menu_tab4 = Menu(self.canvas_tab4_container, tearoff=0)
        self.context_menu_tab4.add_command( label="Save to Image", command=lambda:self.on_save_context_menu_clicked(4))
        
        self.canvas_tab4_container.bind('<Enter>', self._bound_to_mousewheel)
        self.canvas_tab4_container.bind('<Leave>', self._unbound_to_mousewheel)
        
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
        
        Phylo.draw(tree,label_func=self._phyloLabel_callback, do_show=False, axes=ax_tab4)
        
        self.canvas_tab4 = FigureCanvasTkAgg(self.fig_tab4, master=self.canvas_tab4_container)
        self.cwid_tab4 = self.canvas_tab4_container.create_window(0, 0, window=self.canvas_tab4.get_tk_widget(), anchor=Tkconstants.NW)
        self.canvas_tab4_container.config(scrollregion=self.canvas_tab4_container.bbox(Tkconstants.ALL),width=800,height=600)

        #PCOA
        canvas_frame_tab3 = Frame(self.tab3)
        canvas_frame_tab3.grid(row=1, column=1, sticky=Tkconstants.NS)
        canvas_frame_tab3.rowconfigure(1, weight=1)
        canvas_frame_tab3.columnconfigure(1, weight=1)
         
        self.fig_tab3, ax_tab3 = plt.subplots()
        ax_tab3.autoscale(enable=True)
         
        self.canvas_tab3_container = Canvas(canvas_frame_tab3)
        self.canvas_tab3_container.grid(row=1, column=1, sticky=Tkconstants.NSEW)
        x_scroll_tab3 = Scrollbar(canvas_frame_tab3, orient=HORIZONTAL)
        y_scroll_tab3 = Scrollbar(canvas_frame_tab3, orient=VERTICAL)
        x_scroll_tab3.grid(row=2, column=1, sticky=Tkconstants.EW)
        y_scroll_tab3.grid(row=1,column=2, sticky=Tkconstants.NS)
        self.canvas_tab3_container.config(xscrollcommand=x_scroll_tab3.set)
        x_scroll_tab3.config(command=self.canvas_tab3_container.xview)
        self.canvas_tab3_container.config(yscrollcommand=y_scroll_tab3.set)
        y_scroll_tab3.config(command=self.canvas_tab3_container.yview)
        
        self.context_menu_tab3 = Menu(self.canvas_tab3_container, tearoff=0)
        self.context_menu_tab3.add_command( label="Save to Image", command=lambda:self.on_save_context_menu_clicked(3))
        
        self.canvas_tab3_container.bind('<Enter>', self._bound_to_mousewheel)
        self.canvas_tab3_container.bind('<Leave>', self._unbound_to_mousewheel)
        
        mds = MDS(n_components=2, dissimilarity='precomputed')
        mdsDim = mds.fit_transform(genomeDist)
        ax_tab3.scatter(mdsDim[:,0], mdsDim[:,1], s = 50., color='blue')
        for i, txt in enumerate(genomeList):
            ax_tab3.annotate(txt, (mdsDim[i,0], mdsDim[i,1]),fontsize=8)
        plt.grid(True)
        
        self.canvas_tab3 = FigureCanvasTkAgg(self.fig_tab3, master=self.canvas_tab3_container)
        self.cwid_tab3 = self.canvas_tab3_container.create_window(0, 0, window=self.canvas_tab3.get_tk_widget(), anchor=Tkconstants.NW)
        self.canvas_tab3_container.config(scrollregion=self.canvas_tab3_container.bbox(Tkconstants.ALL),width=800,height=600)

        #heatmap
        canvas_frame_tab2 = Frame(self.tab2)
        canvas_frame_tab2.grid(row=1, column=1, sticky=Tkconstants.NS)
        canvas_frame_tab2.rowconfigure(1, weight=1)
        canvas_frame_tab2.columnconfigure(1, weight=1)
         
        self.fig_tab2, ax_tab2 = plt.subplots()
        ax_tab2.autoscale(enable=True)
         
        self.canvas_tab2_container = Canvas(canvas_frame_tab2)
        self.canvas_tab2_container.grid(row=1, column=1, sticky=Tkconstants.NSEW)
        x_scroll_tab2 = Scrollbar(canvas_frame_tab2, orient=HORIZONTAL)
        y_scroll_tab2 = Scrollbar(canvas_frame_tab2, orient=VERTICAL)
        x_scroll_tab2.grid(row=2, column=1, sticky=Tkconstants.EW)
        y_scroll_tab2.grid(row=1,column=2, sticky=Tkconstants.NS)
        self.canvas_tab2_container.config(xscrollcommand=x_scroll_tab2.set)
        x_scroll_tab2.config(command=self.canvas_tab2_container.xview)
        self.canvas_tab2_container.config(yscrollcommand=y_scroll_tab2.set)
        y_scroll_tab2.config(command=self.canvas_tab2_container.yview)
        
        self.context_menu_tab2 = Menu(self.canvas_tab2_container, tearoff=0)
        self.context_menu_tab2.add_command( label="Save to Image", command=lambda:self.on_save_context_menu_clicked(2))
        
        self.canvas_tab2_container.bind('<Enter>', self._bound_to_mousewheel)
        self.canvas_tab2_container.bind('<Leave>', self._unbound_to_mousewheel)
        
        im = ax_tab2.pcolor(genomeDist_normed, cmap=genome_cmap, vmin=genome_min, vmax=genome_max, edgecolors='white')
        cbar = self.fig_tab2.colorbar(im)
        ax_tab2.set_frame_on(False)
        ax_tab2.set_aspect('equal')
        ax_tab2.invert_yaxis()
        ax_tab2.set_xticks(np.arange(genomeDist_normed.shape[1]) + 0.5, minor=False)
        ax_tab2.set_xticklabels(genomeList,rotation=60,fontsize=8)
        ax_tab2.set_yticks(np.arange(genomeDist_normed.shape[1]) + 0.5, minor=False)
        ax_tab2.set_yticklabels(genomeList)
        
        self.canvas_tab2 = FigureCanvasTkAgg(self.fig_tab2, master=self.canvas_tab2_container)
        self.cwid_tab2 = self.canvas_tab2_container.create_window(0, 0, window=self.canvas_tab2.get_tk_widget(), anchor=Tkconstants.NW)
        self.canvas_tab2_container.config(scrollregion=self.canvas_tab2_container.bbox(Tkconstants.ALL),width=800,height=600)
    
        #network
        canvas_frame_tab1 = Frame(self.tab1)
        canvas_frame_tab1.grid(row=1, column=1, sticky=Tkconstants.NS)
        canvas_frame_tab1.rowconfigure(1, weight=1)
        canvas_frame_tab1.columnconfigure(1, weight=1)
         
        self.fig_tab1, ax_tab1 = plt.subplots()
        ax_tab1.autoscale(enable=True)
         
        self.canvas_tab1_container = Canvas(canvas_frame_tab1)
        self.canvas_tab1_container.grid(row=1, column=1, sticky=Tkconstants.NSEW)
        x_scroll_tab1 = Scrollbar(canvas_frame_tab1, orient=HORIZONTAL)
        y_scroll_tab1 = Scrollbar(canvas_frame_tab1, orient=VERTICAL)
        x_scroll_tab1.grid(row=2, column=1, sticky=Tkconstants.EW)
        y_scroll_tab1.grid(row=1,column=2, sticky=Tkconstants.NS)
        self.canvas_tab1_container.config(xscrollcommand=x_scroll_tab1.set)
        x_scroll_tab1.config(command=self.canvas_tab1_container.xview)
        self.canvas_tab1_container.config(yscrollcommand=y_scroll_tab1.set)
        y_scroll_tab1.config(command=self.canvas_tab1_container.yview)
        
        self.context_menu_tab1 = Menu(self.canvas_tab1_container, tearoff=0)
        self.context_menu_tab1.add_command( label="Save to Image", command=lambda:self.on_save_context_menu_clicked(1))
        
        self.canvas_tab1_container.bind('<Enter>', self._bound_to_mousewheel)
        self.canvas_tab1_container.bind('<Leave>', self._unbound_to_mousewheel)
        
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
        nx.draw_networkx_nodes(G, pos, node_color='b') 
        nx.draw_networkx_labels(G, pos,labelsNet, font_size=8)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        plt.axis('off')
        
        self.canvas_tab1 = FigureCanvasTkAgg(self.fig_tab1, master=self.canvas_tab1_container)
        self.cwid_tab1 = self.canvas_tab1_container.create_window(0, 0, window=self.canvas_tab1.get_tk_widget(), anchor=Tkconstants.NW)
        self.canvas_tab1_container.config(scrollregion=self.canvas_tab1_container.bbox(Tkconstants.ALL),width=800,height=600)
        
        self.zoomin_button.config(state='normal')
        self.zoomout_button.config(state='normal')
        self.save_button.config(state='normal')

    def callCAFE(self):
        if self.optSys == 'Linux' :
            exePath = "./cafe_linux"
            os.system("chmod 777 " + exePath)
        elif self.optSys == 'Darwin' :
            exePath = "./cafe_mac"
            os.system("chmod 777 " + exePath)
        elif self.optSys == 'Windows' :
            #exePath = os.path.join(os.getcwd(),"cafe_win.exe")
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
        
        if self.optSys == 'Darwin':
            p = subprocess.Popen([" ".join(paramArr)], shell=True)
            p.communicate()
        else:
            #dirName = os.path.dirname(os.path.abspath(__file__))
            #dirName = os.getcwd()
            #p = subprocess.Popen(paramArr, cwd=dirName, shell=True)
            
            p = subprocess.Popen(paramArr, shell=True)
            p.communicate()
            
            #f = open('run.cmd', 'w+')
            #f.write(" ".join(paramArr))
            #f.close()
            
            #shell32 = ctypes.windll.shell32
            #ret = shell32.ShellExecuteW(None, u"runas", u"cmd", u"/user:Administrator", None, 1)
            #ret = shell32.ShellExecuteW(None, u"open", u"cmd.exe", u"/C run.cmd", None, 1)
            #print("return ShellExecuteW is: " + str(ret))
        
        #vizURL = "result."+self.distVal.get()+".phylip"
        #self.callViz(vizURL)
    
    def on_run_button_clicked(self,event=None):
        if not len(self.input_list) or len(self.input_list)==1:
            tkMessageBox.showerror(APP_NAME,'Program cannot proceed when the input size <= 1')
            return
        
        if len(self.log_list) > 0:
            self.console_box.delete(1.0, 'end')
            self.log_list = []
        
        if self.optSys == 'Windows' :
            os.system("del " + self.logURL)
            os.system("del " + self.logURL_replicate)  
        else:
            os.system("rm " + self.logURL)
            os.system("rm " + self.logURL_replicate)  
        
        self.checkConsoleT = threading.Thread(target=self.display_console)
        self.checkConsoleT.setDaemon(True)    
        self.checkConsoleT.start()
        
        #self.callT = threading.Thread(target=self.callCAFE)
        #self.callT.setDaemon(True)    
        #self.callT.start()
        
        self.run_button.config(state='disabled')
        self.callCAFE()
        
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
        or self.distVal.get() == 'Russel' or self.distVal.get() == 'Sneath' or self.distVal.get() == 'Tanimoto' or self.distVal.get() == 'Yule' \
        or self.distVal.get() == 'FFP' or self.distVal.get() == 'Co-phylog' :
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='disabled')
            self.mVal_spinbox.config(from_=0)
            self.mVal_spinbox.config(to=0)
            self.mVal.set(0)
        elif self.distVal.get() == 'D2star' or self.distVal.get() == 'D2shepp':
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='normal')
            max_mVal = int(self.kVal.get())-1
            self.mVal_spinbox.config(from_=-1)
            self.mVal_spinbox.config(to=max_mVal)
            self.mVal.set(-1)
        elif self.distVal.get() == 'CVtree':
            self.kVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(state='disabled')
            max_mVal = int(self.kVal.get())-2
            self.mVal_spinbox.config(to=max_mVal)
            self.mVal_spinbox.config(from_=max_mVal)
            self.mVal.set(max_mVal)
        elif self.distVal.get() == 'JS':
            self.kVal_spinbox.config(state='disabled')
            self.mVal_spinbox.config(state='normal')
            self.mVal_spinbox.config(to=20)
            self.mVal_spinbox.config(from_=1)
            self.mVal.set(8)
        else:
            print("No suitable distance!")
    
    def on_k_changed(self,event=None):
        if self.distVal.get() == 'D2star' or self.distVal.get() == 'D2shepp':
            max_mVal = int(self.kVal.get())-1
            self.mVal_spinbox.config(from_=-1)            
            self.mVal_spinbox.config(to=max_mVal)
        elif self.distVal.get() == 'CVtree':
            max_mVal = int(self.kVal.get())-2
            self.mVal_spinbox.config(from_=max_mVal)
            self.mVal_spinbox.config(to=max_mVal)
        else:
            pass

if __name__ == '__main__':  
    
#     if not os.environ["TCL_LIBRARY"]:
#         os.environ["TCL_LIBRARY"] = os.path.join(os.getcwd(), 'tcl')
#         print(os.environ["TCL_LIBRARY"])
#     if not os.environ["TK_LIBRARY"]:
#         os.environ["TK_LIBRARY"] = os.path.join(os.getcwd(), 'tk')
#         print(os.environ["TK_LIBRARY"])

    print(os.getcwd())
    os.chdir(os.getcwd())
    print(os.environ["TCL_LIBRARY"])
    
    root = Tk()
    root.style = Style()
    root.style.theme_use("clam")
    
    app = GUIApp(root)
    root.geometry('{}x{}'.format(1280, 800))
    root.resizable(False, False)
    root.iconbitmap('image/logo.ico')
    root.mainloop()