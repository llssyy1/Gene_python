import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import re
import warnings
# from tkinter import *
import optparse

#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='input gtf files',help='enter your transcript gtf)')
parse.add_option('-t','--txt',dest='txt',action='store',metavar='txt file name',help='please enter your txt files')
parse.add_option('-i','--id',dest='id',action='store',metavar='',help='input gene ensg id')
(options,args) = parse.parse_args()
#outPutFileName = options.outfile
GTF = options.gtf
TXT = options.txt
GE_ID = options.id

class showGene(object):

    def __init__(self, gtfText_path, Gene_name, idtext_path):
        self.genename = ''
        self.genedata = ''
        self.genenum = ''
        # print(idtext_path)
        self.genedata = gtfText_path
        self.genename = Gene_name
        self.show(idtext_path)

    def show(self, idtext_path):
        global genetxt_path
        global ax
        traNum = -1
        transcript_num = 0
        gene_exist = 0
        transcript_list = []
        # trNumList = []
        line_width = 5
        if self.genedata == '' or self.genename == '':
            return

        data_path = str(os.getcwd()) + '/' + self.genedata
        genename = self.genename
        result_path = os.getcwd() + '/My result/' + genename
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        genetxt_path = result_path + '/' + genename + '.txt'

        f = open(data_path, 'r')
        fp = open(genetxt_path, 'w')
        allLines = f.readlines()

        for eachLine in allLines:
            if not eachLine.startswith('#'):
                eachLine_list = eachLine.split('\t')
                name = re.search('gene_id "(.*?)";', eachLine_list[-1]).group(1)
                if name == genename:
                    gene_exist = 1
                    fp.write(eachLine)
                    if eachLine_list[2] == 'exon':
                        transcript_name = re.search('transcript_id "(.*?)";' , eachLine_list[-1]).group(1)
                        if transcript_name not in transcript_list:
                            transcript_list.append(transcript_name)
                            transcript_num += 1
        f.close()
        fp.close()
        if transcript_num >= 21:
            line_width = 4
        lista = []
        for eachnum in transcript_list:
            try:
                trLabel = '(' + str(self.read_transcriptNum(idtext_path, eachnum)) + ')'
                lista.append(trLabel)
            except:
                print("cannot find the %s file",idtext_path)
                trLabel = '(1)'
                lista.append(trLabel)
        # print(lista)
        if gene_exist == 0:
            print("genename no exist,please check and retry")
        else:
            self.paint(result_path,genename,transcript_num,transcript_list,lista,line_width,traNum)

    def paint(self,result_path,genename,transcript_num,transcript_list,lista,line_width,traNum):
        fig = plt.figure(1)
        num = 0
        pdf_path = result_path + '/' + genename + '.pdf'

        fp = open(genetxt_path, 'r')
        allLines = fp.readlines()

        startList = []
        endList = []

        trStart = []
        trEnd = []
        for eachLine in allLines:
            eachLine_list = eachLine.split('\t')
            startList.append(eachLine_list[3])
            endList.append(eachLine_list[4])
        transcript_start = min(startList)
        transcript_end = max(endList)

        tssList = self.readHG38('hg38.cage_peak_phase1and2combined_fair_ann.txt.gz.extract.tsv'
                      , str(eachLine_list[0]), eachLine_list[6], transcript_start, transcript_end)

        for eachLine in allLines:
            traNum = traNum + 1
            eachLine_list = eachLine.split('\t')

            trStart.append(eachLine_list[3])
            trEnd.append(eachLine_list[4])

            tName = re.search('transcript_id "(.*?)";', allLines[traNum]).group(1)
            if (traNum + 1) < len(allLines):
                tName1 = re.search('transcript_id "(.*?)";', allLines[traNum + 1]).group(1)
            else:
                tName1 = ""
            if eachLine_list[2] == 'exon':

                if eachLine_list[6] == '+':
                    arr = '4'
                else:
                    arr = '3'

                ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])

                ax.set_xlim(int(transcript_start) - 10, int(transcript_end) + 10)
                ax.set_ylim(-0.5, transcript_num + 1)
                ax.set_xticks(np.linspace(int(transcript_start) - 10, int(transcript_end) + 10, 2))
                ax.set_yticks([0.1] + list(range(1, transcript_num + 1)))
                ax.set_yticklabels(['tss'] + transcript_list)
                ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
                (int(transcript_start) - 10, int(transcript_end) + 10, 2)]])

                ax.spines['top'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.get_xaxis().set_tick_params(direction='out')
                ax.tick_params(axis=u'y', which=u'both', length=0)
                ax1 = ax.twinx()
                ax1.set_ylim(-0.5, transcript_num + 1)
                ax1.set_yticks([0.1] + list(range(1, transcript_num + 1)))
                ax1.set_yticklabels([' '] + lista)
                ax1.spines['top'].set_visible(False)
                ax1.spines['left'].set_visible(False)
                ax1.spines['right'].set_visible(False)
                ax1.tick_params(axis=u'y', which=u'both', length=0)
                for label in ax1.get_yticklabels():
                    label.set_fontsize(6)

                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(6)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(6)

                if num == 0:
                    # print("start = %s,end = %s"%(transcript_start,transcript_end))
                    line1 = [int(transcript_start), num + 0.1], [int(transcript_end), num + 0.1]
                    (line1_xs, line1_ys) = zip(*line1)
                    ax.add_line(lines.Line2D(line1_xs, line1_ys, linewidth=0.5, c='black'))
                    for eachtss in tssList:
                        eachtss_array = eachtss.split(',')
                        print("tssList: %s, %s " % (eachtss_array[0], eachtss_array[1]))
                        line = [(int(eachtss_array[0]) - 1, num+0.1), (int(eachtss_array[1]) + 1, num+0.1)]
                        (line_xs, line_ys) = zip(*line)
                        ax.add_line(lines.Line2D(line_xs, line_ys,
                                                 solid_capstyle='butt', solid_joinstyle='miter',
                                                 linewidth=int(line_width), alpha=1,
                                                 antialiased=False, color='red'))
                    num = num + 1
                line1 = [(int(eachLine_list[3]) - 1, num), (int(eachLine_list[4]) + 1, num)]
                (line1_xs, line1_ys) = zip(*line1)
                ax.add_line(lines.Line2D(line1_xs, line1_ys,
                                         solid_capstyle='butt', solid_joinstyle='miter',
                                         linewidth=int(line_width), alpha=1,
                                         antialiased=False, color='black'))

                if tName != tName1:
                    minLine = int(min(trStart))
                    maxLine = int(max(trEnd))
                    addtext = (int(transcript_end)-int(transcript_start)) / 30
                    while (minLine < maxLine):
                        ax.scatter(minLine,num, marker=arr, linewidth=0.1,c='blue',s = 15,alpha=0.5)
                        minLine +=addtext
                    trline = [(int(min(trStart))-10, num),
                              (int(max(trEnd))+10, num)]
                    (trline_xs, trline_ys) = zip(*trline)
                    ax.add_line(lines.Line2D(trline_xs, trline_ys,
                                             solid_capstyle='butt', solid_joinstyle='miter',
                                             linewidth=0.5, alpha=1,
                                             antialiased=False, color='black'))
                    num = num + 1
                    trStart.clear()
                    trEnd.clear()

            fig.suptitle('\n\n\nchr' + str(eachLine_list[0]) + ': ' + genename, fontsize=10)
        # print("transcript_start = %s,transcript_end = %s"%(transcript_start,transcript_end))
        fig.savefig(pdf_path, dpi=150)
        # sys.exit(0)

    def read_transcriptNum(self,text_path,geneName):
        total = 0;
        data_num = str(os.getcwd()) + '/' + text_path
        f = open(data_num, 'r')
        allLines_num = f.readlines()
        trNumList = []
        for each_num in allLines_num:
            if not each_num.startswith('#'):
                each_num_list = each_num.split('\t')
                trname = each_num_list[1]
                if trname == geneName:
                    trNumList.append(each_num_list[0])
        for each_num in trNumList:
            # 'gene_name "(.*?)";'
            numlist = re.findall(r"\d+",re.search('/f(.*?)p',each_num).group())
            num = int(numlist[0])
            total = total + num
        f.close()
        return total

    def readHG38(self,text_path,charname,flag,start,end):
        list = []
        stendList = []
        data_num = os.getcwd() + '/' + text_path
        f = open(data_num, 'r')
        allLines_num = f.readlines()
        for each_num in allLines_num:
            each_num_list = each_num.split('\t')
            text = each_num_list[0]
            chr = re.findall('chr(.*?):',text)
            if charname in chr:
                if flag in text:
                    list1 = re.findall('chr'+ charname +':(.*?),', text)
                    list.append(list1)
        for leng in list:
            str = leng[0]
            each_num_list = str.split('..')
            start1 = each_num_list[0]
            end1 = each_num_list[1]
            # print(" %s is %s" % (int(start1) >= int(start),int(end1)<=int(end)))
            if int(start1) >= int(start) and int(end1) <= int(end):
                standend = start1 + ',' +end1
                stendList.append(standend)
        # print(stendList)
        f.close()
        return stendList

def main():
    d = showGene(GTF,GE_ID,TXT)
    # mainloop()

if __name__ == '__main__':
    main()
