#!/usr/bin/env python
#! coding:utf-8
"""
Usage:
 gff2topsl.py <queryChromSizes> <targetChromSizes> <inGff> <out.psl>

convert a GFF2 file to a PSL file.
Input file must conform to the GFF2 specification:
         http://gmod.org/wiki/GFF2

Duo Xie (BGI,xieduo@genomics.cn) 2018-12-22

Arguments:
 <queryChromSizes>    file with query (main coordinates) CDS sizes,such as ENSGALP00000000002<tab>2022.
 <targetChromSizes>   file with target (Target attribute)  chromosome sizes,such as chr1<tab>195276750.
 <inGff>              GFF2 formatted file.
 <out.psl>            PSL formatted output.
Options:
 -h,--help

"""
import re
from docopt import docopt

def readSize(ChromSizes):
        f=open(ChromSizes,"r")
        lines=f.readlines()
        size={}
        for l in lines:
                l=l.rstrip("\n")
                arr=l.split("\t")
                size[arr[0]]=arr[1]
        return size
        f.close()
def gffToPsl(gff,qChromSize,tChromSize,outPsl):
        f=open(gff,"r")
        psl=open(outPsl,'w')
        lines=f.readlines()
        count=0
        for l in lines:
                l=l.rstrip("\n")
                arr=l.split("\t")
                count+=1
                if l.find("mRNA")!=-1:
                        blockSizes=list()
                        blockCount=0
                        tBaseInsert=0
                        tStarts=list()
                        qStarts=list()
                        blockEnd=list()
                        tName=arr[0]
                        qName=re.search(r'ID=([^;]+);',arr[-1]).group(1)
                        strand=arr[6]
                        tStart=int(arr[3])-1
                        tEnd=arr[4]
                else:
                        tStarts.append(str(int(arr[3])-1)+",")
                        blockSizes.append(int(arr[4])-int(arr[3]))
                        blockEnd.append(arr[4])
                        blockCount+=1
                        if len(blockEnd)>1:
                                tBaseInsert+=int(arr[3])-int(blockEnd[-2])+1
                        if len(blockSizes)>1:
                                qStarts.append(str(sum(blockSizes[:-1]))+",")
                        else:
                                qStarts.append("0,")
                        if count==len(lines) or lines[count].find("mRNA")!=-1:
                                blockSizes = [str(x) for x in blockSizes]
                                output=str(qChromSize[qName])+"\t"+str(0)+"\t"  \
                                        +str(0)+"\t"+str(0)+"\t"+str(0)+"\t"+str(0)+"\t" \
                                        +str(blockCount-1)+"\t" \
                                       +str(tBaseInsert)+"\t"+strand+"\t"+qName+"\t"+str(qChromSize[qName])+"\t"   \
                                       +str(0)+"\t"   \
                                       +str(qChromSize[qName])+"\t"+tName+"\t"+str(tChromSize[tName])+"\t" \
                                       +str(tStart)+"\t"+str(tEnd)+"\t"+str(blockCount)+"\t"+",".join(blockSizes)+","+ \
                                       "\t"+"".join(qStarts)+"\t"+"".join(tStarts)+"\n"
                                psl.write(output)
                                       
        f.close()
        psl.close()
if __name__ == '__main__':
        arguments = docopt(__doc__)
        gff=arguments["<inGff>"]
        psl=arguments["<out.psl>"]
        qChromSizes=readSize(arguments["<queryChromSizes>"])
        tChromSizes=readSize(arguments["<targetChromSizes>"])
        gffToPsl(gff,qChromSizes,tChromSizes,psl)
