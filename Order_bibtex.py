# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 13:45:30 2018

@author: siirias
"""
import re
import Tkinter, tkFileDialog

root = Tkinter.Tk()
root.withdraw()

file_name = tkFileDialog.askopenfilename()

#file_name="koebib.bib"
#save_file="koe_ready.bib"
if(file_name!=''):
    data=open(file_name,'r').readlines()
    data="\n".join(data)
    #data="@article{Joku2011, jotain {muuta},\n ja vielä}, %jep joo @article{Aasi1992, olio, kaura{kges sgs}, gesr}  @article{Zenon1211, Geage {sguuta}, ja vielä}"
    
    start_at=-1
    first_open=-1
    open_queue=-1
    closing_mark=-1
    found_refs={}
    for i,j in enumerate(data):
        if(data[i]=='@' and open_queue<=0):
            start_at=i
        if(data[i]=="{"):
            if(open_queue==0):
                first_open=i
                open_queue+=1
                closing_mark=-1
            else:
                open_queue+=1
                
        if(data[i]=="}"):
            if(open_queue==0):
                closing_mark=i
                reference=data[start_at:closing_mark+1]
                print reference
                reference=re.sub("^([^,]*)\s*,\s*([^@])","\g<1>,\n\g<2>",reference)
                reference=re.sub("\n\s*","\n    ",reference)
                reference=re.sub("([^\n])\}$","\g<1>\n}",reference)
                
                found_refs[re.search("{\s*([^,]*?)\s*,",reference).group(1).lower()]=reference
                start_at=-1
                first_open=-1
                open_queue=-1
                closing_mark=-1
                           
            else:
                open_queue-=1
                
    output_data=""
    for i in sorted(found_refs.keys()):
        output_data+=found_refs[i]+"\n\n"
    output_data=re.subn("\n\s*\n\s*([^@])","\n\g<1>",output_data)[0]
    save_file = tkFileDialog.asksaveasfilename()
    if(save_file!=''):
        open(save_file,'w').write(output_data)
