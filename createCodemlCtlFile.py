#!/usr/bin/env python3
"""
This script takes a phylip format alignment file and creates a basic
control file for input to codeml. Also allows user to specify some
option parameters. See paml manual for details of options.
Many will default and need not be provided. All options not implemented
yet.

Input seqfile name in the format (GENENAME.best.nuc_clean.phy).
Output from PRANK (Loytynoja and Goldman 2005) is properly formatted

If an option is not needed pass: '*' to that option.
"""
import csv
import argparse
import re
def codemlfile(args):
    gene = re.split("\.", args['seqfile'])[0]
    if int(args['runmode']) < 0:
        # runmode -2 or -3 runs pairwise models
        # no "model" will be specified
        if int(args['fix_omega']) == 1:
            
            ctlfilename = gene+"_fo_pw"\
                          +"_rm"+str(args['runmode'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
            ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
            args['outfile'] = gene+"_fo_pw"\
                              +"_rm"+str(args['runmode'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
        else:
            ctlfilename = gene+"_eo_pw"\
                          +"_rm"+str(args['runmode'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
            
            ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
            args['outfile'] = gene+"_eo_pw"\
                              +"_rm"+str(args['runmode'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
    elif int(args['runmode']) == 0:
        # any value for runmode >= 0 runs phylogenetic models
        # with either user defined trees or automatic trees
        # runmode = 0: the user specifies a phylogenetic tree topology
        if int(args['fix_omega']) == 1:
            ctlfilename = gene+"_rm"+args['runmode']\
                          +"_fo_Model"+str(args['model'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
            ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
            args['outfile'] = gene+"_rm"+args['runmode']\
                              +"_fo_Model"+str(args['model'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
        else:
            ctlfilename = gene+"_rm"+args['runmode']\
                          +"_eo_Model"+str(args['model'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
            ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
            args['outfile'] = gene+"_rm"+args['runmode']\
                              +"_eo_Model"+str(args['model'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
    else:
        # runmode any other value > 2
        if int(args['fix_omega']) == 1:
            ctlfilename = gene+"_rm"+args['runmode']\
                          +"_fo_Model"+str(args['model'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
            ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
            args['outfile'] = gene+"_rm"+args['runmode']\
                              +"_fo_Model"+str(args['model'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
        else:
           ctlfilename = gene+"_rm"+args['runmode']\
                          +"_eo_Model"+str(args['model'])\
                          +"_NSs"+str(args['NSsites'])\
                          +args['handle']\
                          +".ctl"
           ctlfile = csv.writer(open(ctlfilename,'w'),delimiter = " ")
           args['outfile'] = gene+"_rm"+args['runmode']\
                              +"_eo_Model"+str(args['model'])\
                              +"_NSs"+str(args['NSsites'])\
                              +args['handle']\
                              +".results.txt"
    
    for opt in args.keys():
        if args[opt] == '*':
            pass
        else:
            if opt != 'handle':
                ctlfile.writerow([opt, "=", args[opt]])

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

#MAIN CONTROLS
    parser.add_argument('-seqfile',
                       help = 'sequence file name')

    parser.add_argument('-treefile', default = '*',
                        help = 'phylogenetic tree file name')

    parser.add_argument('-handle',default = "",
                        help = 'handle to add to beginning of output\n\
                        files and .ctl files in case there are several \n\
                        similar models being run')
 #KAPPA CONTROLS    
    parser.add_argument('-fix_kappa', nargs = '?', const = 1,
                        default = 0,
                        help = 'fix kappa? (boolean 0[no]/1[yes])')

    parser.add_argument('-kappa',default = 2,
                        help = 'value for kappa/initial value '+
                        'for kappa')
 #OMEGA CONTROLS
    parser.add_argument('-fix_omega', nargs = '?', const = 1,
                        default = 0,
                       help = 'fix omega? (boolean)')

    parser.add_argument('-omega',default = 0.001,
                       help = 'value for omega/initial value for '+
                       'omega')
 #OUTPUT CONTROLS
    parser.add_argument('-noisy', default = 0,
                        help = 'how much rubbish printed to screen'+
                        '(0-3,9')

    parser.add_argument('-verbose', default = 0,
                        help = 'how much rubbish in result file (0-1)')

    parser.add_argument('-runmode', default = -2,
                        help = 'runmode (-2,-3[...],0-5): (see paml manual)')
 #MOLECULAR EVOLUTION MODEL CONTROLS
    parser.add_argument('-seqtype', default = 1,
                        help = 'seqtype: 1=codons, 2=AAs, '+
                        '3=codons-->AAs')

    parser.add_argument('-CodonFreq', default = 0,
                        help = 'equilibrium codon frequency (0 - 3)')

    parser.add_argument('-model', default = 0,
                        help = 'branch models: 0 = one omega for all'+
                        ' branches')

 #   parser.add_argument('-aaDist', default = '*',
 #                       help = 'AA distances equal (0) or
#                        Granthams matrix (1)')

    parser.add_argument('-NSsites', default = 0,
                        help = 'specifies site model')

    parser.add_argument('-icode', default = 0,
                        help = 'which genetic code (0-11)')

    parser.add_argument('-Mgene', default = 0,
                        help = 'partition models for codon' +
                        ' substitution')
 #ALPHA CONTROLS
 #   parser.add_argument('-fix_alpha', nargs = '?', const = 0,
 #                       default = 1,
 #                       help = 'fix alpha? (boolean) (substitution
 #                       rate variation)')

 #   parser.add_argument('-alpha', default = 0
 #                       help = 'alpha, 0 = infinity, constant rate')
    
 #RHO CONTROLS
    parser.add_argument('-fix_rho', nargs = '?', const = 1,
                        default = 0,
                        help = 'fix correlation of substitution rates'+
                         ' at adjacent sites? (boolean)')

    parser.add_argument('-rho', default = 0,
                        help = 'correlation of substitution rates at'+
                        ' adjacent sites')
#OTHER CONTROLS
    parser.add_argument('-cleandata', nargs = '?', const = 1,
                        default = 0,
                        help = 'remove sites with ambiguity data? '+
                        '(boolean)')
 #   parser.add_argument('-RateAncestor', default = 0,
 #                       help = '0 or 1, 1 forces more analysis')
 
    parser.add_argument('-getSE', nargs = '?', const = 1,
                        default = 0,
                        help = 'estimate SE? (boolean)')
    
 #   parser.add_argument('-Small_Diff', default = 1e-6
 #                       help = 'small value used in differentials')

 #   parser.add_argument('-fix_blength', default = 0,
 #                       help = 'branch length behaviour')

 #   parser.add_argument('-method', default = 0,
 #                       help = 'iteration algorithm to use')

    parser.add_argument('-clock', default = 0,
                        help = 'what kind of molecular clock model')
                        
    args = vars(parser.parse_args())

    try:
        codemlfile(args)        
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass


