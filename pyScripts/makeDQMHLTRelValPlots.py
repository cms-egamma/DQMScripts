#!/usr/bin/env python

import os
import argparse
import ROOT   
import glob
import shutil
import sys

ROOT.gROOT.SetBatch(ROOT.kTRUE)

def find_nth(string, substring, n):
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start+len(substring))
        n -= 1
    return start

class RelValInfo:
    cmssw_version_full = ""
    cmssw_version = ""
    pu_type = ""
    year_stamp = ""
    sample_type = ""
    sample_type_full = ""
    global_tag = ""
    mahi_status = ""
    data_info = ""
    
    def __init__(self,filename):
        filename_split = filename.split("__")
        gt_pu_start = find_nth(filename_split[2],"-",1)
        gt_end = filename_split[2].rfind("-v")
        self.cmssw_version_full = filename_split[2][:gt_pu_start]
        self.cmssw_version = self.cmssw_version_full[6:].replace('_','')

        #now figure out if there is PU which starts at the front of this string if so
        #otherwise its just the GT
        gt_pu_str = filename_split[2][gt_pu_start+1:gt_end]        
        if gt_pu_str[:2]=="PU":
            self.pu_type = gt_pu_str.split('_')[0]
            self.global_tag = gt_pu_str.split('_')[1]
        else:
            self.pu_type = ""
            self.global_tag = gt_pu_str

        self.sample_type_full = filename_split[1]
        self.sample_type = self.sample_type_full[6:-3]

        tmp_str_list = filename_split[2].replace("-","_").split("_")
        for part in tmp_str_list:
            if "mahi" in part:
                self.mahi_status = part
            if "sig" in part:
                self.data_info = part
            if "2017" in part or "2018" in part:
                self.year_stamp = part
    
#        print self.cmssw_version_full,self.cmssw_version,self.pu_type,self.global_tag,self.sample_type_full,self.sample_type


def makeRelValPlots(filename,ref_filename,base_output_dir,update,label1,label2):
    


    sample_info = RelValInfo(filename)
    ref_info = RelValInfo(ref_filename)
    
    subdir_name="EGRelVal_"+sample_info.sample_type+"_"+sample_info.cmssw_version+sample_info.pu_type+sample_info.mahi_status+label1+"Vs"+ref_info.cmssw_version+ref_info.pu_type+ref_info.mahi_status+label2

    output_dir = base_output_dir.rstrip("/")+"/"+subdir_name
    
    leg_entry = sample_info.cmssw_version
    if sample_info.pu_type != "": leg_entry+="-"+sample_info.pu_type
    if sample_info.mahi_status != "": leg_entry+="-"+sample_info.mahi_status
    if label1 != "": leg_entry+="-"+label1
    refleg_entry = ref_info.cmssw_version+ref_info.pu_type
    if ref_info.mahi_status != "": refleg_entry+="-"+ref_info.mahi_status
    if label2 != "": refleg_entry+="-"+ label2

    if os.path.isdir(output_dir):
        if update:
            shutil.rmtree(output_dir)
        else:
            print "\n",output_dir,"exists already, use --update option to delete it and remake it\n"
            sys.exit(0)

    ROOT.gROOT.ProcessLine(".L rootScripts/makeDQMHLTRelValPlots.C+")
    ROOT.printAllPlots(output_dir,filename,ref_filename,leg_entry,refleg_entry)
    
    subdirs = glob.glob(output_dir+"/*")
    for subdir in subdirs:
        files = glob.glob(subdir+"/*.gif")
        with open(subdir+"/index.html","w") as f:
            for filename in files:
                filename = filename.split("/")[-1]
                if subdir.split("/")[-1] == "EGTagAndProbe":
                    html_str = "Path: {} Filter: {} <br>\n".format(filename.split('-')[1],filename.split('-')[2].split(".")[0])
                    f.write(html_str)
                html_str = "<a href=\"{}\"><img class=\"image\" width=\"1000\" src=\"{}\" ALIGH=TOP></a><br><br>\n".format(filename,filename)
                f.write(html_str)
    with open(output_dir+"/index.html","w") as f:
        for subdir in subdirs:
            html_str = "<a href=\"{}\">{}</a><br>\n".format(subdir.split("/")[-1],subdir.split("/")[-1])
            f.write(html_str)
    
       
                

    


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='reads DQM histogram files and produces formated plots for easier validation')
    parser.add_argument('filename',help='filename')
    parser.add_argument('refFilename',help='reference filename')
    parser.add_argument('-o','--outputDir',help='output base directory',required=True)
    parser.add_argument('--update',action='store_true',help='allows overwriting of existing directory')
    parser.add_argument('--l1',help='label of sample 1',default = "")
    parser.add_argument('--l2',help='label of sample 2',default = "")
    args = parser.parse_args()

    makeRelValPlots(args.filename,args.refFilename,args.outputDir,args.update,args.l1,args.l2)
