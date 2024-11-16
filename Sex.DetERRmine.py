#!/usr/bin/env python3
import sys, argparse, json
from math import sqrt
from collections import OrderedDict

VERSION="1.1.3"

def determine_sex(rate_x, rate_y):
    """
    Determine biological sex based on X and Y chromosome coverage ratios
    
    Parameters:
    rate_x: X/Autosome coverage ratio
    rate_y: Y/Autosome coverage ratio
    
    Returns:
    str: "Female", "Male", or "Undetermined"
    """
    # Thresholds based on expected ratios
    # Males: X/A ≈ 0.5, Y/A > 0.1
    # Females: X/A ≈ 1.0, Y/A ≈ 0
    
    if rate_x >= 0.8 and rate_y < 0.05:
        return "Female"
    elif 0.35 <= rate_x <= 0.65 and rate_y >= 0.1:
        return "Male"
    else:
        return "Undetermined"

def CalcErrors(AutSnps, XSnps, YSnps, NrAut, NrX, NrY):
    SNPs=[AutSnps, XSnps, YSnps]
    Reads=[NrAut, NrX, NrY]
    p={}
    ErrNr={}
    dp={}
    Errdp={}
    rate={}
    rateErr={}
    
    Total=sum(Reads)
    
    ## If no sites are covered at all, set error and rate to 0. 
    if Total == 0:
        (rate, rateErr) = ({"Aut":0.0,"X":0.0,"Y":0.0}, {"Aut":0.0,"X":0.0,"Y":0.0})
    else :
        for Bin,Idx in zip(["Aut", "X", "Y"], range(3)):
            p[Bin]=Reads[Idx]/Total
            ErrNr[Bin]=sqrt(Total*p[Bin])
            dp[Bin]=Reads[Idx]/SNPs[Idx]
            Errdp[Bin]=ErrNr[Bin]/SNPs[Idx]
    
        for Bin in ["X","Y"]:
            rate[Bin]=dp[Bin]/dp["Aut"]
            rateErr[Bin]=sqrt((Errdp[Bin]/dp["Aut"])**2 + (Errdp["Aut"]*dp[Bin]/(dp["Aut"]**2))**2)
    
    return (rate, rateErr)

#### MAIN ####

parser = argparse.ArgumentParser(description="Calculate the relative X- and Y-chromosome coverage of data, as well as the associated error bars for each.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input samtools depth file. Omit to read from stdin.", required=False)
parser.add_argument("-f", "--SampleList", type=argparse.FileType('r'), help="A list of samples/bams that were in the depth file. One per line. Should be in the order of the samtools depth output.")
parser.add_argument("-v", "--version", action="store_true", help="Print the version of the script and exit.")
args = parser.parse_args()

if args.version:
    print(VERSION, file=sys.stderr)
    sys.exit(0)

if args.Input == None:
    args.Input = sys.stdin

Names=OrderedDict()
if args.SampleList != None:
    Samples = [line.strip() for line in args.SampleList]
    for idx,Sample in enumerate(Samples):
        Names.update({Sample:idx})
    NrAut  = [0 for x in range(len(Names))]
    NrX    = [0 for x in range(len(Names))]
    NrY    = [0 for x in range(len(Names))]
    
Reads={}
AutSnps=0
YSnps=0
XSnps=0
for line in args.Input:
    fields=line.strip().split()
    if fields [0][0:3]=="chr":
        Chrom=fields[0][3:]
    else:
        Chrom=fields[0]
    if fields[0][0]=="#":
        if args.SampleList==None:
            Zip    = zip(fields[2:],range(len(fields[2:])))
            for Sample,Index in Zip:
                Names.update({Sample:Index})
            NrAut  = [0 for x in range(len(Names))]
            NrX    = [0 for x in range(len(Names))]
            NrY    = [0 for x in range(len(Names))]
            continue
        else:
            continue
    depths=[int(x) for x in fields[2:]]
    if Chrom != "Y" and Chrom != "X":
        AutSnps+=1
    if Chrom == "Y":
        YSnps+=1
    if Chrom == "X":
        XSnps+=1
    for x in Names:
        if Chrom != "Y" and Chrom != "X":
            NrAut[Names[x]]+=depths[Names[x]]
        if Chrom == "Y":
            NrY[Names[x]]+=depths[Names[x]]
        if Chrom == "X":
            NrX[Names[x]]+=depths[Names[x]]

if sum([AutSnps, XSnps, YSnps]) == 0:
    raise IOError("""
Input depth file is empty. Stopping execution.

Please ensure you are using a bed file compatible with your reference genome coordinates.""")

print ("#Sample", "#SnpsAut", "#SNPsX", "#SnpsY", "NrAut", "NrX", "NrY", "x-rate", "y-rate", "Err(x-rate)", "Err(y-rate)", "Sex", sep="\t", file=sys.stdout)
data=OrderedDict()
data['Metadata'] = {'tool_name' : "Sex.DetERRmine", "version" : VERSION}
for Ind in Names:
    rate,rateErr=CalcErrors(AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]])
    sex_estimate = determine_sex(rate["X"], rate["Y"])
    data[Ind] = {
        "Snps Autosomal" : AutSnps, 
        "XSnps" : XSnps, 
        "YSnps": YSnps, 
        "NR Aut" : NrAut[Names[Ind]],
        "NrX": NrX[Names[Ind]], 
        "NrY": NrY[Names[Ind]], 
        "RateX" : rate["X"], 
        "RateY": rate["Y"], 
        "RateErrX" : rateErr["X"], 
        "RateErrY" : rateErr["Y"],
        "Sex": sex_estimate
    }
    print (Ind, AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]], 
           rate["X"], rate["Y"], rateErr["X"], rateErr["Y"], sex_estimate, sep="\t", file=sys.stdout)
    
with open('sexdeterrmine.json', 'w') as outfile:
    json.dump(data, outfile)
