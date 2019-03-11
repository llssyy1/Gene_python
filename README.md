# Gene_python
used to Gene(gtf->picture)
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


class showGene(object):

    def __init__(self):
       #os.listdir （path）返回path目录    os.getcwd（）返回当前目录
        list1 = os.listdir(os.getcwd())  # 列出文件夹下所有的目录与文件
        list2 = []
        for i in range(0, len(list1)):
            if list1[i][-4:-1] + list1[i][-1] == '.gtf':
                list2.append(list1[i])
        # 选择基因组以及要显示的基因的GUI界面
        self.top = Tk()
        self.genename = StringVar()
        self.genedata = StringVar()
        self.getinfoFm = Frame(self.top)
        self.data = Label(self.getinfoFm, text="Chooes the gtf of a species").grid(column=0,
                                                                                   row=0)  # 添加一个标签，并将其列设置为1，行设置为0
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
        self.showB = Button(self.buttonFm, text='SHOW', command=self.show,
                            activeforeground='white',
                            activebackground='blue').grid(column=0, row=0)
        self.quitB = Button(self.buttonFm, text='QUIT', command=self.top.quit,
                            activeforeground='white',
                            activebackground='blue').grid(column=1, row=0)
        self.buttonFm.pack()

    def show(self):
        global genetxt_path
        global ax
        traNum = -1
        transcript_num = 0
        gene_exist = 0
        transcript_list = []
        line_width = 5
        if self.genedata.get() == '' or self.genename.get() == '':
            return

        # gtf文件路径，默认为可执行文件（exe）所在文件夹
        data_path = os.getcwd() + '/' + self.genedata.get()
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
            # print(eachLine)
            if not eachLine.startswith('#'):
                eachLine_list = eachLine.split('\t')
                name = re.search('gene_id "(.*?)";', eachLine_list[-1]).group(1)
                if name == genename:
                    gene_exist = 1
                    fp.write(eachLine)
                  #  print(eachLine_list)
                    if eachLine_list[2] == 'exon':
                        #transcript_num += 1
                        transcript_name = re.search('transcript_id "(.*?)";', eachLine_list[-1]).group(1)
                        if transcript_name not in transcript_list:
                            transcript_list.append(transcript_name)
                            transcript_num += 1
        f.close()
        fp.close()
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
            fig = plt.figure(1)
            num = 1
            # tmp_colors = ['lime', 'red', 'blue', 'yellow', 'yellow', 'w']
            # names_tmp_colors = ['gene', 'CDS', 'exon', 'three_prime_utr', 'five_prime_utr', 'stop_codon']
            tmp_colors = ['yellow']
            names_tmp_colors = ['exon']
            colors_legend_name = ['gene']
            color_dict = dict(zip(names_tmp_colors, tmp_colors))
            png_path = result_path + '/' + genename + '.png'
            # 读取之前保存的包含基因信息的txt文件，若不存在，弹出警告信息
            fp = open(genetxt_path, 'r')
            allLines = fp.readlines()
            #所有转录本的长短
            startList = []
            endList = []
            #记录每一个转录本的长短
            trStart = []
            trEnd = []
            for eachLine in allLines:
                eachLine_list = eachLine.split('\t')
                startList.append(eachLine_list[3])
                endList.append(eachLine_list[4])
            # print(min(startList))
            # print(max(endList))
            transcript_start = min(startList)
            transcript_end = max(endList)
            for eachLine in allLines:
                # print(eachLine)
                traNum = traNum+1
                eachLine_list = eachLine.split('\t')
                #某一个转录本添加长度
                trStart.append(eachLine_list[3])
                trEnd.append(eachLine_list[4])
                #获取transcript_id的名字，当tname和tname1不相同时，画上横线
                # 当基因转录子只有一个时，设置tName1为空#当基因转录子只有一个时，设置tName1为空
                tName = re.search('transcript_id "(.*?)";', allLines[traNum]).group(1)
                if (traNum+1) < len(allLines):
                    tName1 = re.search('transcript_id "(.*?)";', allLines[traNum + 1]).group(1)
                else:
                    tName1 = ""
                #print(tName+tName1)
                if eachLine_list[2] == 'exon':
                    # 判断方向
                    if eachLine_list[6] == '+':
                        arr = '->'
                    else:
                        arr = '<-'
                    #新建区域ax
                    #figure的百分比,从figure 10%的位置开始绘制, 宽高是figure的80%
                    #left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
                    ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])
                    #print(ax)
                    arrow = mpatches.FancyArrowPatch(
                        (int(transcript_start), 0.1),
                        (int(transcript_end), 0.1),
                        arrowstyle=arr,
                        mutation_scale=25, lw=1, color='lime', antialiased=True)
                    # 画箭头
                    ax.add_patch(arrow)
                    # 坐标轴标签
                    ax.set_xlim(int(transcript_start)-10, int(transcript_end)+10)
                    # ax.set_xlim(int(eachLine_list[3]), int(eachLine_list[4]))
                    ax.set_ylim(-0.5, transcript_num + 1)
                    ax.set_xticks(np.linspace(int(transcript_start)-10, int(transcript_end)+10, 2))
                    # ax.set_xticks(np.linspace(int(eachLine_list[3]), int(eachLine_list[4]), 5))
                    ax.set_yticks([0.1] + list(range(1, transcript_num + 1)))
                    ax.set_yticklabels(['gene'] + transcript_list)
                    ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
                    (int(transcript_start)-10, int(transcript_end)+10, 2)]])
                    # ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
                    # (int(eachLine_list[3]), int(eachLine_list[4]), 5)]])
                    # 坐标轴显示
                    ax.spines['top'].set_visible(False)
                    ax.spines['left'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    # ax.spines['bottom'].set_visible(False)
                    ax.get_xaxis().tick_bottom()
                    ax.get_yaxis().tick_left()
                    ax.get_xaxis().set_tick_params(direction='out')
                    ax.tick_params(axis=u'y', which=u'both', length=0)  # 纵坐标刻度线不显示（length=0）
                    # ax.tick_params(axis=u'x', which=u'both', length=0)  # 纵坐标刻度线不显示（length=0）
                    # 坐标轴字体大小
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(6)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(6)

                    # 画外显子---线上的小方块
                    line1 = [(int(eachLine_list[3])-1, num), (int(eachLine_list[4])+1, num)]
                    (line1_xs, line1_ys) = zip(*line1)
                    ax.add_line(lines.Line2D(line1_xs, line1_ys,
                                             solid_capstyle='butt', solid_joinstyle='miter',
                                             linewidth=int(line_width), alpha=1,
                                             antialiased=False, color='black'))
                    # if num != 0:

                    #画基因----水平线
                    if tName != tName1:
                       # print(num)
                        trline = [(int(min(trStart)), num ),
                                  (int(max(trEnd)), num)]
                        (trline_xs, trline_ys) = zip(*trline)
                        ax.add_line(lines.Line2D(trline_xs, trline_ys, linewidth=0.2,
                                                 antialiased=False, color='black'))
                        num = num + 1
                        trStart.clear()
                        trEnd.clear()

                fig.suptitle('\n\n\nchr' + str(eachLine_list[0]) + ': ' + genename, fontsize=10)
                # 保存图片
                fig.savefig(png_path, dpi=150)

            # 显示图片
            self.img_gif = PhotoImage(file=png_path)
            self.label_img = Label(self.top, image=self.img_gif)
            self.clearB = Button(self.top, text='CLEAR', command=self.clear,
                                 activeforeground='white',
                                 activebackground='blue')
            self.clearB.pack()
            self.label_img.pack()

    def clear(self):
        global genetxt_path
        self.top.update()
        self.label_img.destroy()
        self.clearB.destroy()
        # 删除保存基因信息的txt文件，只保留图片
        os.remove(genetxt_path)


def main():
    d = showGene()
    mainloop()


if __name__ == '__main__':
    main()
