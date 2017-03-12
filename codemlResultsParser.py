#!/usr/bin/env python

"""
This script takes a codeml results file from either
pairwise or branch models, checks what model is run and
prints:
Branch models: (runmode = 0, NSsites = 0)
model = 0 -> gene name, STOP codons, dN, dS, dN/dS, lnL
model = 2 -> gene name, STOP codons, (sp1_dN, sp1_dS, sp_1dN/dS,...) lnL
                                         for each spp

Branch-site models: (runmode = 0, NSsites = 0)
model = 0 -> gene name, STOP codons, dN, dS, dN/dS, lnL
model = 2 -> gene name, STOP codons, (sp1_dN, sp1_dS, sp_1dN/dS,...) lnL
                                         for each spp

Gene name and model number must be in file name.

"""
import os
import csv
import re
import argparse
import sys

def parse_codeml_results(fil,runmode,om,model,nssites,gene_name):
    read = open(fil, 'rb')
    lines = read.readlines()
    stop = "0"
####################################################################
    if runmode == "0" and om == "est" and model == "2":
        if nssites == "0":
            # if model = 2 and omega is estimated
            # different dN/dS values on some branches,
            # and 1 site class within a sequence
            # output file will be formatted a particular way.
            w_tree_line = "none"
            ds_tree_line = "none"
            dn_tree_line = "none"
            for i in range(0,len(lines)):
                if "dS tree:\n" in lines[i]:
                    # get tree with ds as branch labels
                    ds_tree_line = i+1
                    ds_tree = lines[ds_tree_line]
                    ds_tree = ds_tree.replace(")", "")
                    ds_tree = ds_tree.replace("\n", "")
                    ds_tree = ds_tree.replace(")", "")
                    ds_tree = ds_tree.replace(";", "")
                    ds_tree = ds_tree.replace(" ", "")
                    ds_tree = re.split(",", ds_tree)
                    ds = []
                    for l in ds_tree:
                        ds.append(re.split(":", l)[1])
                if "dN tree:\n" in lines [i]:
                    # get tree with dn as branch labels
                    dn_tree_line = i+1
                    dn_tree = lines[dn_tree_line]
                    dn_tree = dn_tree.replace(")", "")
                    dn_tree = dn_tree.replace("\n", "")
                    dn_tree = dn_tree.replace("(", "")
                    dn_tree = dn_tree.replace(";", "")
                    dn_tree = dn_tree.replace(" ", "")
                    dn_tree = re.split(",", dn_tree)
                    dn = []
                    for m in dn_tree:
                        dn.append(re.split(":", m)[1])
                if "w ratios as labels for TreeView:\n" in lines[i]:
                    # get tree with omega as branch labels
                    w_tree_line = i+1
                    w_tree = lines[w_tree_line]
                    w_tree = w_tree.replace(")", "")
                    w_tree = w_tree.replace("\n", "")
                    w_tree = w_tree.replace("(", "")
                    w_tree = w_tree.replace(";", "")
                    w_tree = w_tree.replace(" ", "")
                    w_tree = w_tree.replace("#", ":")
                    w_tree = re.split(",", w_tree)
                    dnds = []
                    for n in w_tree:
                        dnds.append(re.split(":", n)[1])
                if "lnL(ntime: " in lines[i]:
                    # get lnL estimate
                    lnL_line = lines[i]
                    lnL_values = re.split(":", lnL_line)[-1]
                    np_values = re.split(":",lnL_line)[-2]
                    lnL_values = lnL_values.replace(" ", "")
                    lnL = float(re.split("\+", lnL_values)[0])
                    np = np_values.replace(" ","")
                    np = np_values.replace(")","")
                if "???" in lines[i]:
                    # check for premature stop codons
                    stop = "1"
                if "branch" in lines[i]:
                    NSdata_s = i
                    NSdata_e = i+5
                    lines2 = lines[NSdata_s:NSdata_e]
                    NSdata = re.sub("\s",",", lines2[2])
                    NSdata = re.split(",+",NSdata)
                    N = NSdata[3]
                    S = NSdata[4]
            if w_tree_line == "none" or ds_tree_line == "none"\
               or dn_tree_line == "none":
                # no calculations performed by CODEML. Sequences ambiguous.
                print gene_name,",","no_data"
            else:
                print gene_name,",",stop,",",N,",",S,",",\
                      ",".join(dn),",",",".join(ds),",",\
                      ",".join(dnds),",",lnL,",",np,",",model,
        if nssites == "2":
            #print "BRANCH-SITE MODEL"
            # if model = 2 and omega is estimated
            # different dN/dS values on some branches,
            # and 4 site classes within a sequence
            # only two branch types are allowed: foreground and background
            # output file will be formatted a particular way.
            for i in range(0,len(lines)):
                if "lnL(ntime: " in lines[i]:
                    # get lnL estimate and np (number of parameters)
                    lnL_line = lines[i]
                    lnL_line = re.split(":", lnL_line)
                    lnL_values = lnL_line[-1]
                    lnL_values = lnL_values.replace(" ", "")
                    lnL = float(re.split("\+", lnL_values)[0])
                    np = lnL_line[2]
                    np = np.replace(" ","")
                    np = np.replace(")","")
                if "???" in lines[i]:
                    # check for premature stop codons
                    stop = "1"
                if "site class " in lines[i]:
                    # get the omega estimates for each class
                    # and the proportion of sites that belong
                    # to each class.
                    #print lines[i+1]
                    prop_line = lines[i+1]
                    prop_line = prop_line.replace("\n","")
                    prop_line = prop_line.replace(" ",",")
                    prop_line = re.split(",+",prop_line)
                    #print lines[i+2]
                    bkgrndw_line = lines[i+2]
                    bkgrndw_line = bkgrndw_line.replace("\n","")
                    bkgrndw_line = bkgrndw_line.replace(" ",",")
                    bkgrndw_line = re.split(",+",bkgrndw_line)
                    #print lines[i+3]
                    frgrndw_line = lines[i+3]
                    frgrndw_line = frgrndw_line.replace("\n","")
                    frgrndw_line = frgrndw_line.replace(" ",",")
                    frgrndw_line = re.split(",+",frgrndw_line)
                    #print prop_line,bkgrndw_line,frgrndw_line
                    c0p=prop_line[1]
                    c1p=prop_line[2]
                    c2ap=prop_line[3]
                    c2bp=prop_line[4]
                    c0bw=bkgrndw_line[2]
                    c1bw=bkgrndw_line[3]
                    c2abw=bkgrndw_line[4]
                    c2bbw=bkgrndw_line[5]
                    c0fw=frgrndw_line[2]
                    c1fw=frgrndw_line[3]
                    c2afw=frgrndw_line[4]
                    c2bfw=frgrndw_line[5]

            print gene_name,\
                  ",",c0p,",",c1p,",",c2ap,",",c2bp,\
                  ",",c0bw,",",c1bw,",",c2abw,",",c2bbw,\
                  ",",c0fw,",",c1fw,",",c2afw,",",c2bfw,\
                  lnL,",",np,",",stop,",",model,
####################################################################
    elif runmode == "0" and om == "fixed" and model == "2":
        #print "NOT WRITTEN THIS PART OF PARSER YET"
        if nssites == "0":
            #print "BRANCH MODEL"
            # if model = 2 and omega is estimated
            # different dN/dS values on some branches,
            # and 1 site class within a sequence
            # output file will be formatted a particular way.
            w_tree_line = "none"
            ds_tree_line = "none"
            dn_tree_line = "none"
            for i in range(0,len(lines)):
                if "dS tree:\n" in lines[i]:
                    # get tree with ds as branch labels
                    ds_tree_line = i+1
                    ds_tree = lines[ds_tree_line]
                    ds_tree = ds_tree.replace(")", "")
                    ds_tree = ds_tree.replace(")", "")
                    ds_tree = ds_tree.replace("\n", "")
                    ds_tree = ds_tree.replace(";", "")
                    ds_tree = ds_tree.replace(" ", "")
                    ds_tree = re.split(",", ds_tree)
                    ds = []
                    for l in ds_tree:
                        ds.append(re.split(":", l)[1])
                if "dN tree:\n" in lines [i]:
                    # get tree with dn as branch labels
                    dn_tree_line = i+1
                    dn_tree = lines[dn_tree_line]
                    dn_tree = dn_tree.replace(")", "")
                    dn_tree = dn_tree.replace("(", "")
                    dn_tree = dn_tree.replace("\n", "")
                    dn_tree = dn_tree.replace(";", "")                    
                    dn_tree = dn_tree.replace(" ", "")
                    dn_tree = re.split(",", dn_tree)
                    dn = []
                    for m in dn_tree:
                        dn.append(re.split(":", m)[1])
                if "w ratios as labels for TreeView:\n" in lines[i]:
                    # get tree with omega as branch labels
                    w_tree_line = i+1
                    w_tree = lines[w_tree_line]
                    w_tree = w_tree.replace(")", "")
                    w_tree = w_tree.replace("(", "")
                    w_tree = w_tree.replace("\n", "")
                    w_tree = w_tree.replace(";", "")
                    w_tree = w_tree.replace(" ", "")
                    w_tree = w_tree.replace("#", ":")
                    w_tree = re.split(",", w_tree)
                    dnds = []
                    for n in w_tree:
                        dnds.append(re.split(":", n)[1])
                if "lnL(ntime: " in lines[i]:
                    # get lnL estimate
                    lnL_line = lines[i]
                    lnL_values = re.split(":", lnL_line)[-1]
                    np_values = re.split(":",lnL_line)[-2]
                    lnL_values = lnL_values.replace(" ", "")
                    lnL = float(re.split("\+", lnL_values)[0])
                    np = np_values.replace(" ","")
                    np = np_values.replace(")","")
                if "???" in lines[i]:
                    # check for premature stop codons
                    stop = "1"
                if "branch" in lines[i]:
                    NSdata_s = i
                    NSdata_e = i+5
                    lines2 = lines[NSdata_s:NSdata_e]
                    NSdata = re.sub("\s",",", lines2[2])
                    NSdata = re.split(",+",NSdata)
                    N = NSdata[3]
                    S = NSdata[4]
            if w_tree_line == "none" or ds_tree_line == "none"\
               or dn_tree_line == "none":
                # no calculations performed by CODEML. Sequences ambiguous.
                print gene_name,",","no_data"
            else:
                print gene_name,",",stop,",",N,",",S,",",\
                      ",".join(dn),",",",".join(ds),",",\
                      ",".join(dnds),",",lnL,",",np,",",model,
        if nssites == "2":
            #print "BRANCH-SITE MODEL"
            # if model = 2 and omega is estimated
            # different dN/dS values on some branches,
            # and 4 site classes within a sequence
            # only two branch types are allowed: foreground and background
            # output file will be formatted a particular way.
            for i in range(0,len(lines)):
                if "lnL(ntime: " in lines[i]:
                    # get lnL estimate and np (number of parameters)
                    lnL_line = lines[i]
                    lnL_line = re.split(":", lnL_line)
                    lnL_values = lnL_line[-1]
                    lnL_values = lnL_values.replace(" ", "")
                    lnL = float(re.split("\+", lnL_values)[0])
                    np = lnL_line[2]
                    np = np.replace(" ","")
                    np = np.replace(")","")
                if "???" in lines[i]:
                    # check for premature stop codons
                    stop = "1"
                if "site class " in lines[i]:
                    # get the omega estimates for each class
                    # and the proportion of sites that belong
                    # to each class.
                    #print lines[i+1]
                    prop_line = lines[i+1]
                    prop_line = prop_line.replace("\n","")
                    prop_line = prop_line.replace(" ",",")
                    prop_line = re.split(",+",prop_line)
                    #print lines[i+2]
                    bkgrndw_line = lines[i+2]
                    bkgrndw_line = bkgrndw_line.replace("\n","")
                    bkgrndw_line = bkgrndw_line.replace(" ",",")
                    bkgrndw_line = re.split(",+",bkgrndw_line)
                    #print lines[i+3]
                    frgrndw_line = lines[i+3]
                    frgrndw_line = frgrndw_line.replace("\n","")
                    frgrndw_line = frgrndw_line.replace(" ",",")
                    frgrndw_line = re.split(",+",frgrndw_line)
                    #print prop_line,bkgrndw_line,frgrndw_line
                    c0p=prop_line[1]
                    c1p=prop_line[2]
                    c2ap=prop_line[3]
                    c2bp=prop_line[4]
                    c0bw=bkgrndw_line[2]
                    c1bw=bkgrndw_line[3]
                    c2abw=bkgrndw_line[4]
                    c2bbw=bkgrndw_line[5]
                    c0fw=frgrndw_line[2]
                    c1fw=frgrndw_line[3]
                    c2afw=frgrndw_line[4]
                    c2bfw=frgrndw_line[5]

            print gene_name,",",c0p,",",c1p,",",c2ap,",",c2bp,\
                  ",",c0bw,",",c1bw,",",c2abw,",",c2bbw,\
                  ",",c0fw,",",c1fw,",",c2afw,",",c2bfw,\
                  lnL,",",np,",",stop,",",model,

####################################################################
    elif runmode == "0" and om == "est" and model == "0":
        # one dN/dS across the whole tree but still estimated
        omega = "none"
        if nssites == "0":
            # one site class
            TREE = "none"
            for i in range(0, len(lines)):
                # Find the TREE line to know what the tip numbers are
                if "TREE #" in lines[i]:
                    TREE = lines[i].replace("TREE #  1:  ","")
                    TREE = re.sub(";.*","",TREE)
                    TREE = TREE.replace("(","")
                    TREE = TREE.replace(")","")
                    TREE = TREE.replace(" ","")
                    TREE = TREE.replace("\n","")
                    TREE = TREE.split(",")
                    #print TREE
                # Find the dN,dS,S,N data for each tip in TREE
                if "tree length for dN:" in lines[i]:
                    dn = lines[i]
                    dn = dn.replace(" ", "")
                    dn = dn.replace("\n", "")
                    dn = re.split(":", dn)[1]
                if "tree length for dS:" in lines[i]:
                    ds = lines[i]
                    ds = ds.replace(" ", "")
                    ds = ds.replace("\n", "")
                    ds = re.split(":", ds)[1]
                if "omega (dN/dS)" in lines[i]:
                    omega_line = i
                    omega_dat = lines[omega_line]
                    omega_dat = omega_dat.replace("\n", "")
                    omega_dat = omega_dat.replace(" ", "")
                    omega_dat = re.split("=", omega_dat)
                    omega = float(omega_dat[1])
                if "lnL(ntime: " in lines[i]:
                    lnL_line = lines[i]
                    lnL_values = re.split(":", lnL_line)[-1]
                    np_values = re.split(":",lnL_line)[-2]
                    lnL_values = lnL_values.replace(" ", "")
                    lnL = float(re.split("\+", lnL_values)[0])
                    np = np_values.replace(" ","")
                    np = np_values.replace(")","")
                if "???" in lines[i]:
                    stop = "1"
                # Find the table of dn,ds,N and S data for each tip
                # and node. Parse this table.
                if "branch" in lines[i][:10]:
                    NSdata_s = i
                    NSdata_e = i+5
                    lines2 = lines[NSdata_s:]
                    sp_dn = []
                    sp_ds = []
                    for line2 in lines2:
                        for tip in TREE:
                            if ".."+tip+" " in line2:
                                #print line2
                                line2=re.sub("\s",",", line2)
                                line2=re.split(",+",line2)
                                sp_dn.append(line2[-5])
                                sp_ds.append(line2[-4])
                    NSdata = re.sub("\s",",", lines2[2])
                    NSdata = re.split(",+",NSdata)
                    N = NSdata[3]
                    S = NSdata[4]
            if omega == "none":
                # no calculations performed by CODEML. Sequences ambiguous
                print gene_name,",","no_data"
            else:
                print gene_name,",",stop,",",N,",",S,",",dn,",",ds,","\
                  ,omega,",",lnL,",",np,",",model,",",",".join(sp_dn),",",\
                  ",".join(sp_ds)
        #There are many more possibilities when model=0
        if nssites == "3":
            print "SITES MODEL (NO BRANCHES)"
        if nssites == "4":
            print "SITES MODEL (NO BRANCHES)"
        if nssites == "5":
            print "SITES MODEL (NO BRANCHES)"
        if nssites == "6":
            print "SITES MODEL (NO BRANCHES)"
        if nssites == "7":
            print "SITES MODEL (NO BRANCHES)"

####################################################################
    elif runmode == "0" and om == "fixed" and model == "0":
        # one dN/dS across the whole tree and fixed
        print "NOT WRITTEN THIS PART OF PARSER YET"
        omega = 1
        for i in range(0, len(lines)):
            if "lnL" in lines[i]:
                data = lines[i:]
                data = [d for d in data if d != "\n"]
                lnl = re.split("  ",data[0])[4]
                if lnl == '':
                    lnl = float(re.split(":",re.split("  ",data[0])[3])[1])
                else:
                    lnl = float(re.split("  ",data[0])[4])
            if "???" in lines[i]:
                stop = "1"
        print gene_name,",",stop,",",lnl
####################################################################
    elif runmode == "-3":
        # Pairwise Bayesian estimation of dN/dS
        ml_omega = "none"
        b_omega = "none"
        b_p = "none"
        b = "0"
        for i in range(0,len(lines)):
            if "lnL =" in lines[i]:
                data = lines[i:]
                data = [d for d in data if d != "\n"]
                ml_data = [data[0],data[2]]
                ml_data_o = re.sub("\s\s+"," ",ml_data[1])
                ml_data_o = re.split("[ ]",ml_data_o)
                ml_omega = ml_data_o[ml_data_o.index("dN/dS=")+1]
                N = ml_data_o[ml_data_o.index("N=")+1]
                S = ml_data_o[ml_data_o.index("S=")+1]
                dn = ml_data_o[ml_data_o.index("dN")+1]
                if dn == "=":
                    dn = ml_data_o[ml_data_o.index("dN")+2]
                    dn = re.sub("\n","",dn)
                    dn = re.sub("=","",dn)
                else:
                    dn = re.sub("\n","",dn)
                    dn = re.sub("=","",dn)
                try:
                    ds = ml_data_o[ml_data_o.index("dS")+2]
                    ds = re.sub("\n","",ds)
                except IndexError:
                    ds = ml_data_o[ml_data_o.index("dS")+1]
                    ds = re.sub("\n","",ds)
                    ds = re.sub("=","",ds)
                lnL = float(re.split("=",ml_data[0])[1].replace("\n",""))
                b = "0"
                for line in data:
                    if "Bayesian" in line:
                        # bayesian estimates have not been calculated
                        # print "NA" for values
                        b = "1"
                if b == "1":
                    b_data = re.sub("\s+"," ",data[-1])
                    b_data = re.split("[ ]",b_data)
                    b_omega = float(b_data[2])
                    b_p = float(b_data[-1])
                    
            if "???" in lines[i]:
                stop = "1" 
        if ml_omega == "none":
                # no calculations performed by CODEML.
                # Sequences ambiguous
            print gene_name,",","no_data"            
        else:
            print gene_name,",",stop,",",N,",",S,",",dn,",",ds,","\
                  ,ml_omega,",",lnL,",",b_omega,",",b_p,",",model
####################################################################
    elif runmode == "-2":
        # Pairwise ML estimation of dN/dS
        omega = "none"
        for i in range(0,len(lines)):
            if "lnL =" in lines[i]:
                data = lines[i:]
                data = [d for d in data if d != "\n"]
                data = [data[0],data[2]]
                data_o = re.sub("\s\s+"," ",data[1])
                data_o = re.split("[ ]",data_o)
                omega = data_o[data_o.index("dN/dS=")+1]
                N = data_o[data_o.index("N=")+1]
                S = data_o[data_o.index("S=")+1]
                dn = data_o[data_o.index("dN")+2]
                ds = data_o[data_o.index("dS")+2]
                ds = re.sub("\n","",ds)
                lnL = float(re.split("=",data[0])[1].replace("\n",""))
            if "???" in lines[i]:
                stop = "1" 
        if omega == "none":
                # no calculations performed by CODEML.
                # Sequences ambiguous
            print gene_name,",","no_data"
        else:
            print gene_name,",",stop,",",N,",",S,",",dn,",",ds,","\
                  ,omega,",",lnL,",",model

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    parser.add_argument('-fil',
                       help = 'paml results file')
    
    args = vars(parser.parse_args())

    try:
        filename = re.split("_",args['fil'][:-12])
        #print filename #script tester line
        RM = str([rm for rm in filename if "rm" in rm])
        NS = str([ns for ns in filename if "NS" in ns])
        M = str([m for m in filename if "Model" in m])
        #print RM, NS, M# Script tester line
        if "0" in RM:
            # User specified tree
            runmode = "0"
            if "2" in M:
                # Model=2 multiple w estimated
                model = "2"
                if "2" in NS:
                    nssites = "2"
                if "0" in NS:
                    nssites = "0"
                if "fo" in filename:
                    om = "fixed"
                if "eo" in filename:
                    om = "est"
                    
            elif "0" in M:
                # A single w estimated
                model = "0"
                if "2" in NS:
                    nssites = "2"
                if "0" in NS:
                    nssites = "0"
                if "fo" in filename:
                    om = "fixed"
                if "eo" in filename:
                    om = "est"

        if "-3" in RM:
            # Bayesian pairwise dN/dS
            runmode = "-3"
            if "fo" in filename:
                om = "fixed"
            if "eo" in filename:
                om = "est"
            
        if "-2" in RM:
            # Pairwise dN/dS
            runmode = "-2"
            if "fo" in filename:
                om = "fixed"
            if "eo" in filename:
                om = "est"
                
        gene_name = filename[0] # get gene name from results file name
        #print filename, om,model,nssites,gene_name # Script tester line
        parse_codeml_results(args['fil'],runmode,om,model,nssites,gene_name)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
 
