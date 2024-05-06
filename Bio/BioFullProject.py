# Importing Libs
import pandas as pd
import numpy as np
import customtkinter as ctk
import tkinter as tk
from tkinter import ttk
import bisect
from itertools import permutations,zip_longest, islice


# Algorithm Functions
def Complement(seq):
    dic = {"G": "C", "C": "G", "A": "T", "T": "A"}

    # Filter out non-nucleotide characters
    seq = ''.join(char for char in seq if char in dic)

    s = list(seq)
    for i in range(len(seq)):
        s[i] = str(dic[s[i]])
    s = "".join(s)
    return s

def Reverse(seq):
    s = list(seq)
    s = reversed(s)
    s = "".join(s)
    return s

def Reverse_Complement(seq):
    seq = Reverse(seq)
    seq = Complement(seq)
    return seq
#################################################
# Translate
def Translation_Table(seq):
    dic = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
    s=""
    sf=""
    flag=0
    for i in range(0,len(seq)-2,3):
        if dic[seq[i:i+3]]=="M":
            flag=1
        elif (dic[seq[i:i+3]]=="*"):
            flag=0
        if flag==1:
            s+=dic[seq[i:i+3]]
        sf+=dic[seq[i:i+3]]
    return sf,s

# Matching
def match(seq,sub_seq):
    x=-1
    for i in range(len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
    return x

#################################################

# Bad Characters
def Badchars(seq,sub_seq):
    table=np.zeros([4,len(sub_seq)])     
    row=["C","G","A","T"]
    sumk = 0
    for i in range (4):
        num=-1
        for j in range (len(sub_seq)):
            if row[i]==sub_seq[j]:
                table[i,j]=-1
                num=-1
            else:
                num+=1
                table[i,j]=num
    x=-1
    i=0
    while(i<len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
        else:
            for j in range(i+len(sub_seq)-1,i-1,-1):
                if seq[j] != sub_seq[int(j-i)]:
                    k=row.index(seq[j])
                    i+=table[k,j-i]
                    sumk=sumk+k
                    break
        i=int(i+1)
    return x,table,sumk

#################################################

# Python program for KMP Algorithm
def computeLPSArray(pat, M, lps):
    len = 0 # length of the previous longest prefix suffix
 
    lps[0] # lps[0] is always 0
    i = 1
 
    # the loop calculates lps[i] for i = 1 to M-1
    while i < M:
        if pat[i]== pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            # This is tricky. Consider the example.
            # AAACAAAA and i = 7. The idea is similar 
            # to search step.
            if len != 0:
                len = lps[len-1]
 
                # Also, note that we do not increment i here
            else:
                lps[i] = 0
                i += 1
                
#################################################                

def KMPSearch(txt, pat):
    M = len(pat)
    N = len(txt)
 
    # create lps[] that will hold the longest prefix suffix 
    # values for pattern
    lps = [0]*M
    j = 0 # index for pat[]
 
    # Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps)
 
    i = 0 # index for txt[]
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1
 
        if j == M:
            answer = str(i-j)
            j = lps[j-1]
 
        # mismatch after j matches
        elif i < N and pat[j] != txt[i]:
            # Do not match lps[0..lps[j-1]] characters,
            # they will match anyway
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return answer

#################################################
 
def calcHamming(seqA, seqB):
    
    #Calculate the Hamming distance between two sequences (as strings)of equal length.
    distH = 0
    if len(seqA) == len(seqB):
        for index, value in enumerate(seqA):
            if value != seqB[index]:
                distH += 1
    elif len(seqA) != len(seqB):
        tk.messagebox.showerror("Error", "Error: Sequences must be the same length!")
    return distH

from Bio import SeqIO
from numpy import zeros

#################################################

def indel(sequence, i):
    return sequence[:i] + '-' + sequence[i:]


def edit_alignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)

    edit_distance_matrix = zeros((m + 1, n + 1), dtype=int)
    trace_back_matrix = zeros((m + 1, n + 1), dtype=int)

    for i in range(0, m + 1):
        for j in range(0, n + 1):
            if j == 0:
                edit_distance_matrix[i][0] = i
            elif i == 0:
                edit_distance_matrix[0][j] = j
            else:
                if (not (seq1[i - 1] is seq2[j - 1])):
                    edit_distance_matrix[i - 1][j - 1] += 1
                scores = [edit_distance_matrix[i - 1][j - 1], edit_distance_matrix[i - 1][j] + 1,
                          edit_distance_matrix[i][j - 1] + 1]
                edit_distance_matrix[i][j] = min(scores)
                trace_back_matrix[i][j] = scores.index(edit_distance_matrix[i][j])

    edited_seq1 = seq1
    edited_seq2 = seq2
    i = m
    j = n

    while i > 0 and j > 0:
        if trace_back_matrix[i][j] == 1:
            i -= 1
            edited_seq2 = indel(edited_seq2, j)
        elif trace_back_matrix[i][j] == 2:
            j -= 1
            edited_seq1 = indel(edited_seq1, i)
        else:
            i -= 1
            j -= 1
            
    return edited_seq1, edited_seq2,str(edit_distance_matrix[m][n])

#################################################

def dynamicEditDistance(s, t):
    rows = len(s)+1
    cols = len(t)+1
    dist = [[0 for x in range(cols)] for x in range(rows)]

    # source prefixes can be transformed into empty strings 
    # by deletions:
    for i in range(1, rows):
        dist[i][0] = i

    # target prefixes can be created from an empty source string
    # by inserting the characters
    for i in range(1, cols):
        dist[0][i] = i
        
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0
            else:
                cost = 1
            dist[row][col] = min(dist[row-1][col] + 1,      # deletion
                                 dist[row][col-1] + 1,      # insertion
                                 dist[row-1][col-1] + cost) # substitution
            
    return dist[row][col]

#################################################

def IndexSorted(seq,ln):
    index = []
    for i in range(len(seq)-ln+1):
        index.append((seq[i:i+ln], i))
    index.sort() 
    return index

def query(t,p,index):
    keys = [r[0] for r in index]
    st = bisect.bisect_left(keys,p[:len(keys[0])])
    
    en = bisect.bisect(keys,p[:len(keys[0])])
    hits = index[st:en] 
    l=[h[1] for h in hits ]
    offsets=[]
    for i in l:
        if t[i:i+len(p)]==p:
            offsets.append(i)
    return offsets

#################################################

def overlap(a,b,min_length=3):
    start=0
    while True:
        start=a.find(b[:min_length],start)
        if start==-1:
            return 0
        if b.startswith(a[start:]):
            return len(a[start:])
        start+=1
        

###########

def native_overlap(reads,k):
    olap={}
    for a,b in permutations(reads,2):
        olen=overlap(a, b,k)
        if olen>0:
            olap[(a,b)]=olen
    return olap

#################################################

def Distance(seq,seq2):
    tag=3
    dic={}
    for i in range(0,len(seq)-tag):
        dic[seq[i:i+tag]]=dic.get(seq[i:i+tag],0)+1
    
    dic2={}
    for i in range(0,len(seq2)-tag):
        dic2[seq2[i:i+tag]]=dic2.get(seq2[i:i+tag],0)+1
        
    k=list(dic.keys())
    for i in range(len(k)):
        dic2[k[i]]=(dic2.get(k[i],0)-dic[k[i]])
    d=list(dic2.values())
    Sum=0
    for i in range(len(d)):
        Sum+=d[i]**2
    distance=Sum**(0.5)
    return distance

#################################################

def Suffix_Array(T):
    l=[]
    for i in range(len(T)):
        l.append(T[i:])
    l2=l.copy()
    l.sort()
    dec={}
    for i in range(len(l)):
        dec[l[i]]=i
    table=[]
    for i in range(len(l)):
        table.append([l2[i],i,dec[l2[i]]])
    return table    


def to_int_keys(l):
    index = {v: i for i, v in enumerate(sorted(set(l)))}
    return [index[v] for v in l]

#suffix Matrix construction
def suffix_matrix(s):
    n = len(s)
    k = 1
    line = to_int_keys(s)
    ans = [line]
    while max(line) < n - 1:
        line = to_int_keys(
            list(zip_longest(line, islice(line, k, None),
                             fillvalue=-1)))
        ans.append(line)
        k <<= 1
    return ans

#################################################

# Applying The Selected Algorithm
def apply_algorithm(selected):
    text = input_text.get("1.0", "end-1c")  # Retrieve the text
    pattern = pattern_text.get("1.0","end-1c")      # Retrieve the pattern
    if pattern == '' or pattern == 'Pattern Here':               # Ensuring No empty pattern
        pattern = '1'
    output_text.configure(state="normal")  # Enable modification to set result
    output_text.delete("1.0", "end")  # Clear previous result
    if text == '':
        tk.messagebox.showerror("Error", "Please Enter a String")     # Ensuring that there is no empty text inputted
    elif selected == "DNA complementary":                             # Applying DNA Complementary
        complement_val = selected_radio.get()
        if complement_val == "complement":
            result = Complement(text)
            output_text.insert("end", result)
        elif complement_val == "reverse":
            result = Reverse(text)
            output_text.insert("end", result)
        elif complement_val == "both":
            result = Reverse_Complement(text)
            output_text.insert("end", result)
            
    elif selected == "Translate":                                     # Applying Translation
        result1,result2 = Translation_Table(text)
        output_text.insert("end", result1+'\n\n'+result2)
        
    elif selected == "Bad character algorithm":                       # Applying BadChar with its table
        result,table,skips = Badchars(text,pattern)
        table = pd.DataFrame(table)
        tree["column"] = list(table.columns)
        tree["show"] = "headings"
        for column in tree["columns"]:
            tree.heading(column, text=column) # let the column heading = column name

        df_rows = table.to_numpy().tolist() # turns the dataframe into a list of lists
        
        for row in df_rows:
            tree.insert("", "end", values=row)
        output_text.insert("end", str(result) + '\n' + str(skips))
        
        
    elif selected == "Query and Indexing":
        length = int(len_text.get("1.0","end-1c"))
        result = query(text,pattern,IndexSorted(text,length))
        output_text.insert("end", result)
        
    elif selected == "KMP":
        result = KMPSearch(text,pattern)
        output_text.insert("end", result)
        
    elif selected == "Hamming":
        result = calcHamming(text,pattern)
        output_text.insert("end", result)
        
    elif selected == "EditAlignment":
        result1,result2,result3 = edit_alignment(text,pattern)
        output_text.insert("end",result1+'\n\n'+result2+'\n\n'+result3)
        
    elif selected == "Dynamic Programming":
        result = dynamicEditDistance(text,pattern)
        output_text.insert("end",str(result))
        
    elif selected == "Overlap":
        length = int(len_text.get("1.0","end-1c"))
        result1 = overlap(text,pattern,length)
        reads = [text,pattern]
        result2 = str(native_overlap(reads,length))
        output_text.insert("end", str(result1) + '\n' + result2)
        
    elif selected == "Distance":
        result = Distance(text,pattern)
        output_text.insert("end", result)
    elif selected == "Suffix Array":
        result = Suffix_Array(text)
        table = suffix_matrix(text)
        table = pd.DataFrame(table)
        tree["column"] = list(table.columns)
        tree["show"] = "headings"
        for column in tree["columns"]:
            tree.heading(column, text=column) # let the column heading = column name
            
        df_rows = table.to_numpy().tolist() # turns the dataframe into a list of lists
        
        for row in df_rows:
            tree.insert("", "end", values=row)
            
        output_text.insert("end", result)
    
    output_text.configure(state="disabled")  # Disable modification after setting result


####################################################################################################################################

# GUI Functions
def show_table(checked):
    if checked =="on":
        pre_window.deiconify()  # Show the window
        tree_frame.pack(padx=10,pady=10)
        tree.pack(padx=10,pady=10)
    elif checked =="off":    # Checkbox is unmarked
        pre_window.withdraw()  # Hide the window
    
        
def toggle_visibility(selected):
    
    if selected == "DNA complementary":
        complement.grid(row=4, column=2, pady=10)
        reverse.grid(row=5, column=2, pady=10)
        both.grid(row=6, column=2, pady=10)
        table.grid_forget()
        len_text.grid_forget()
        
    elif selected =="Bad character algorithm":
        table.grid(row=7, column=2, padx=20, pady=(0, 20), sticky="we")
        complement.grid_forget()
        reverse.grid_forget()
        both.grid_forget()
        len_text.grid_forget()
    elif selected =="Suffix Array":
        table.grid(row=7, column=2, padx=20, pady=(0, 20), sticky="we")
        complement.grid_forget()
        reverse.grid_forget()
        both.grid_forget()
        len_text.grid_forget()
        
    elif selected == "Query and Indexing":
        len_text.grid(row=1, column=2, padx=20, pady=20, sticky="nsew")
        complement.grid_forget()
        reverse.grid_forget()
        both.grid_forget()
        table.grid_forget()
    elif selected == "Overlap":
        len_text.grid(row=1, column=2, padx=20, pady=20, sticky="nsew")
        complement.grid_forget()
        reverse.grid_forget()
        both.grid_forget()
        table.grid_forget()
    else:
        complement.grid_forget()
        reverse.grid_forget()
        both.grid_forget()
        table.grid_forget()
        len_text.grid_forget()

def clear_text():
    input_text.delete(1.0, "end-1c")
    pattern_text.delete(1.0, "end-1c")
    len_text.delete(1.0, "end-1c")
    output_text.configure(state="normal") 
    output_text.delete(1.0, "end-1c")
    output_text.configure(state="disabled")  
  
    
def browse_fasta_file():
    # Open a file dialog to select a FASTA file
    global file_path
    file_path = tk.filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta;*.fa")])

    # Check if a file was selected
    if file_path:
        # Read the contents of the FASTA file
        with open(file_path, 'r') as file:
            fasta_content = file.read()

        # Update the text in the text box or do whatever you want with the content
        input_text.delete(1.0, "end")  # Clear previous content
        input_text.insert("end", fasta_content)
        
        
def display_dict_in_textbox(my_dict):
    for key, value in my_dict.items():
        output_text.insert(tk.END, f"{key}: {value}\n")


####################################################################################################################################


# Starting Gui
app = ctk.CTk()
app.title("DNA Processing")
app.after(201, lambda :app.iconbitmap('F:\AI\Projects\Bio\dna.ico'))
app.geometry("1000x600")
app.resizable(0,0)
app.pack_propagate(False)


# Building GUI

# Defining right frame
right_frame = ctk.CTkFrame(app, width=100, height=300,fg_color=('#242424'))
right_frame.grid(column=1,sticky='ns')

# Input Text
input_text = ctk.CTkTextbox(app, width=80, height=10)   
input_text.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
input_text.insert(1.0,"Text Here")

output_text = ctk.CTkTextbox(app, width=80, height=10, state="disabled")
output_text.grid(row=3, column=0, padx=20, pady=20, sticky="nsew")
output_text.configure(state='normal')
output_text.insert(1.0,"Output Here")
output_text.configure(state='disabled')

pattern_text = ctk.CTkTextbox(right_frame,width=20, height=10)
pattern_text.insert(1.0,"Pattern Here")
pattern_text.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")

len_text = ctk.CTkTextbox(right_frame,width=20, height=10)
len_text.insert(1.0,"Length Here")

# Apply Buttons
button_frame1 = ctk.CTkFrame(app, width=200, height=200,fg_color=('#242424'))
button_frame1.grid(row=1,column=0)

apply_button = ctk.CTkButton(button_frame1, text="Apply Algorithm", command=lambda: apply_algorithm(selected_menu.get()))
apply_button.grid(row=1, column=0, pady=10,padx=20, sticky="ew")

browse_button = ctk.CTkButton(button_frame1, text="Browse FASTA File", command=browse_fasta_file)
browse_button.grid(row=1, column=1, pady=10,padx=20, sticky="ew")

button_frame2 = ctk.CTkFrame(app, width=200, height=200,fg_color=('#242424'))
button_frame2.grid(row=4,column=0)

clear_button = ctk.CTkButton(button_frame2, text="Clear", command=clear_text)
clear_button.grid(row=4, column=0, pady=10,padx=20, sticky="ew")

exit_button = ctk.CTkButton(button_frame2, text="Exit", command= lambda:app.destroy())
exit_button.grid(row=4, column=1, pady=10,padx=20, sticky="ew")



# Apply OptionMenu
selected_menu = ctk.StringVar(value="Select Algorithm")
optionmenu = ctk.CTkOptionMenu(right_frame, values=["DNA complementary","Translate",
                                            "Bad character algorithm","Suffix Array",
                                            "Query and Indexing","KMP",
                                            "Hamming","EditAlignment",
                                            "Dynamic Programming","Overlap",
                                            'Distance'],
                                command=toggle_visibility,
                                variable=selected_menu)
optionmenu.grid(row=2, column=2, padx=20, pady=10, sticky="ew")

# Apply RadioButton
selected_radio = ctk.StringVar(value="complement")
complement = ctk.CTkRadioButton(right_frame, text="Complement", variable=selected_radio, value="complement")
reverse = ctk.CTkRadioButton(right_frame, text="Reverse", variable=selected_radio, value="reverse")
both = ctk.CTkRadioButton(right_frame, text="Both", variable=selected_radio, value="both")

# Apply CheckBox
check_var = ctk.StringVar(value="off")
table = ctk.CTkCheckBox(right_frame, text="Show Preprocessing Table", command=lambda:show_table(check_var.get()),
                                     variable=check_var, onvalue="on", offvalue="off", hover_color=('white'),
                                     text_color=('light Blue'))

# Apply TreeView
pre_window = ctk.CTkToplevel(app)
pre_window.title("Preprocessing Table")
pre_window.geometry("500,500")
pre_window.withdraw()

tree_frame = tk.LabelFrame(pre_window,height=500, width=500)
tree = ttk.Treeview(tree_frame)
treescrolly = tk.Scrollbar(tree_frame, orient="vertical", command=tree.yview) # command means update the yaxis view of the widget
treescrollx = tk.Scrollbar(tree_frame, orient="horizontal", command=tree.xview) # command means update the xaxis view of the widget
tree.configure(xscrollcommand=treescrollx.set, yscrollcommand=treescrolly.set) # assign the scrollbars to the Treeview Widget
treescrollx.pack(side="bottom", fill="x") # make the scrollbar fill the x axis of the Treeview widget
treescrolly.pack(side="right", fill="y") # make the scrollbar fill the y axis of the Treeview widget

# Configure grid weights to make resizing work properly
app.grid_columnconfigure(0, weight=1)
app.grid_rowconfigure(0, weight=1)
app.grid_rowconfigure(3, weight=1)
app.mainloop()