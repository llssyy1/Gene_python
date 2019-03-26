import re
import os
import os.path
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from tkinter import *
from tkinter import ttk
import time
import imageio
import shutil
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


class showGene(object):

    def __init__(self,g_id,t_id,exonName,text_path):
       #os.listdir （path）返回path目录    os.getcwd（）返回当前目录
        list1 = os.listdir(os.getcwd())  # 列出文件夹下所有的目录与文件
        list2 = []
        list3 = []
        for i in range(0, len(list1)):
            if list1[i][-4:-1] + list1[i][-1] == '.gtf':
                list2.append(list1[i])
            if list1[i][-4:-1] + list1[i][-1] == '.txt':
                list3.append(list1[i])
        # 选择基因组以及要显示的基因的GUI界面
        self.top = Tk()
        self.genename = StringVar()
        self.genedata = StringVar()
        self.genenum = StringVar()
        self.getinfoFm = Frame(self.top)
        self.data = Label(self.getinfoFm, text="Chooes the gtf").grid(column=0,row=0)  # 添加一个标签，并将其列设置为1，行设置为0
        self.gene = Label(self.getinfoFm, text="       Enter a genename:").grid(column=1, row=0)
        self.choosedata = ttk.Combobox(self.getinfoFm, width=15, textvariable=self.genedata)
        self.choosedata['values'] = [''] + list2
        self.choosedata.current(0)
        self.choosedata.grid(column=0, row=1)
        self.genenameEntry = Entry(self.getinfoFm, textvariable=self.genename, width=12, )
        self.genenameEntry.grid(column=1, row=1)

        self.getinfoFm.pack()
        # 按钮界面
        self.buttonFm = Frame(self.top)
        self.showB = Button(self.buttonFm, text='SHOW', command=lambda:self.show(g_id,t_id,exonName,text_path),
                            activeforeground='white',
                            activebackground='blue').grid(column=0, row=0)
        self.quitB = Button(self.buttonFm, text='QUIT', command=self.top.quit,
                            activeforeground='white',
                            activebackground='blue').grid(column=1, row=0)
        self.buttonFm.pack()

    def show(self,g_id,t_id,exonName,text_path):
        global genetxt_path
        global ax
        parttern = g_id + ' "(.*?)";'
        parttern1 = t_id + ' "(.*?)";'
        traNum = -1
        transcript_num = 0
        gene_exist = 0
        transcript_list = []
        # trNumList = []
        line_width = 5
        if self.genedata.get() == '' or self.genename.get() == '':
            return

        # gtf文件路径，默认为可执行文件（exe）所在文件夹
        data_path = os.getcwd() + '/' + self.genedata.get()
        # data_num = os.getcwd() + '/' + self.genenum.get()
        genename = self.genename.get()
        # 储存各基因结果（txt文本和png图片）的文件夹
        result_path = os.getcwd() + '/My result/' + genename
        # 若路径不存在，创造路径
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        # 临时文件（基因的信息）的路径
        genetxt_path = result_path + '/' + genename + '.txt'

        f = open(data_path, 'r')
        fp = open(genetxt_path, 'w')
        allLines = f.readlines()

        for eachLine in allLines:
            if not eachLine.startswith('#'):
                eachLine_list = eachLine.split('\t')
                name = re.search(parttern, eachLine_list[-1]).group(1)
                if name == genename:
                    gene_exist = 1
                    fp.write(eachLine)
                    if eachLine_list[2] == exonName:
                        transcript_name = re.search(parttern1 , eachLine_list[-1]).group(1)
                        if transcript_name not in transcript_list:
                            transcript_list.append(transcript_name)
                            transcript_num += 1
        f.close()
        fp.close()
        if transcript_num >= 21:
            line_width = 4
        # 读取转录子总数
        lista = []
        for eachnum in transcript_list:
            trLabel = '(' + str(self.read_transcriptNum(text_path, eachnum)) + ')'
            lista.append(trLabel)
        # print(lista)
        # 若gtf文件中不存在相应基因的信息，弹出提示
        if gene_exist == 0:
            self.failedFm = Frame(self.top)
            self.faillb = Label(self.failedFm, text="genename no exist,please check and retry").pack()
            self.failedFm.pack()
            self.top.update()
            time.sleep(3)
            self.genename.set('')
            self.failedFm.destroy()
            shutil.rmtree(result_path)
        # 反之向下执行
        else:
            png_path = self.paint(result_path,genename,parttern1,exonName,transcript_num,transcript_list,lista,line_width,traNum)
            # 显示图片
            self.img_gif = PhotoImage(file=png_path)
            self.label_img = Label(self.top, image=self.img_gif)

            self.clearB = Button(self.top, text='CLEAR', command=self.clear,
                                 activeforeground='white',
                                 activebackground='blue')
            self.clearB.pack()
            self.label_img.pack()

    def paint(self,result_path,genename,parttern1,exonName,transcript_num,transcript_list,lista,line_width,traNum):
        fig = plt.figure(1)
        num = 0
        png_path = result_path + '/' + genename + '.png'
        pdf_path = result_path + '/' + genename + '.pdf'
        # 读取之前保存的包含基因信息的txt文件，若不存在，弹出警告信息
        fp = open(genetxt_path, 'r')
        allLines = fp.readlines()
        # 所有转录本的长短
        startList = []
        endList = []
        # 记录每一个转录本的长短
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
            # 某一个转录本添加长度
            trStart.append(eachLine_list[3])
            trEnd.append(eachLine_list[4])
            # 获取transcript_id的名字，当tname和tname1不相同时，画上横线
            # 当基因转录子只有一个时，设置tName1为空#当基因转录子只有一个时，设置tName1为空
            tName = re.search(parttern1, allLines[traNum]).group(1)
            if (traNum + 1) < len(allLines):
                tName1 = re.search(parttern1, allLines[traNum + 1]).group(1)
            else:
                tName1 = ""
            if eachLine_list[2] == exonName:
                # 判断方向
                if eachLine_list[6] == '+':
                    arr = '4'
                else:
                    arr = '3'
                # 新建区域ax
                ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])
                # 坐标轴标签
                ax.set_xlim(int(transcript_start) - 10, int(transcript_end) + 10)
                ax.set_ylim(-0.5, transcript_num + 1)
                ax.set_xticks(np.linspace(int(transcript_start) - 10, int(transcript_end) + 10, 2))
                ax.set_yticks([0.1] + list(range(1, transcript_num + 1)))
                ax.set_yticklabels(['tss'] + transcript_list)
                ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
                (int(transcript_start) - 10, int(transcript_end) + 10, 2)]])
                # 坐标轴显示
                ax.spines['top'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.get_xaxis().set_tick_params(direction='out')
                ax.tick_params(axis=u'y', which=u'both', length=0)  # 纵坐标刻度线不显示（length=0）
                ax1 = ax.twinx()  # 创建第二个坐标轴
                ax1.set_ylim(-0.5, transcript_num + 1)
                ax1.set_yticks([0.1] + list(range(1, transcript_num + 1)))
                ax1.set_yticklabels([' '] + lista)
                ax1.spines['top'].set_visible(False)
                ax1.spines['left'].set_visible(False)
                ax1.spines['right'].set_visible(False)
                ax1.tick_params(axis=u'y', which=u'both', length=0)  # 纵坐标刻度线不显示（length=0）
                for label in ax1.get_yticklabels():
                    label.set_fontsize(6)
                # 坐标轴字体大小
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(6)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(6)

                # 画外显子---线上的小方块
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
                # if num != 0:

                # 画基因----水平线
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
        # 保存图片
        # print("transcript_start = %s,transcript_end = %s"%(transcript_start,transcript_end))
        fig.savefig(png_path, dpi=140)
        fig.savefig(pdf_path, dpi=150)
        return png_path

    def clear(self):
        global genetxt_path
        self.top.update()
        self.label_img.destroy()
        self.clearB.destroy()
        # 删除保存基因信息的txt文件，只保留图片
        os.remove(genetxt_path)

    # 计算id2id.txt文件的转录子数目
    def read_transcriptNum(self,text_path,geneName):
        total = 0;
        data_num = os.getcwd() + '/' + text_path
        fn = open(data_num, 'r')
        allLines_num = fn.readlines()
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
        fn.close()
        return total

    def readHG38(self,text_path,charname,flag,start,end):
        list = []
        stendList = []
        data_num = os.getcwd() + '/' + text_path
        f = open(data_num, 'r')
        allLines_num = f.readlines()
        for each_num in allLines_num:
            each_num_list = each_num.split('\t')
            text = each_num_list[0]      #筛选留下第一列
            chr = re.findall('chr(.*?):',text)
            if charname in chr:      #通过染色体号筛选
                if flag in text:     #通过基因正负方向筛选
                    list1 = re.findall('chr'+ charname +':(.*?),', text)
                    list.append(list1)   # 留下染色体的长度
        for leng in list:
            str = leng[0]
            each_num_list = str.split('..')
            start1 = each_num_list[0]
            end1 = each_num_list[1]
            # print("输入大小：start = %s  end = %s"% (start1,end1))
            # print(" %s is %s" % (int(start1) >= int(start),int(end1)<=int(end)))
            if int(start1) >= int(start) and int(end1) <= int(end):
                standend = start1 + ',' +end1
                stendList.append(standend)
        print(stendList)
        f.close()
        return stendList

def main():
    d = showGene('gene_id','transcript_id','exon','K510_3rd.id2id.txt')
    # d.read_transcriptNum( 'K510_3rd.id2id.txt','ENST00000372422')
    # d.readHG38('hg38.cage_peak_phase1and2combined_fair_ann.txt.gz.extract.tsv','16','+','29581516','29613725')
    mainloop()


if __name__ == '__main__':
    main()
