from matplotlib import pyplot as plt
from typing import Set
from bisect import insort
from time import time 
from datetime import timedelta
import numpy as np
import math, logging, csv, string
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
plt.rcParams.update({'figure.max_open_warning': 0})

color_dict = {
    ('Charged -',): '#A91FEA',
    ('Charged +',): '#45D318',
    ('Hydrophilic',): 'blue',
    ('Hydrophobic',): 'red'
}

name_dict = {
    '1h55': 'HRP',
    '1cf3': 'GOx',
    '2hrw': 'GFP',
    '1yph': r'$ \alpha $'+'-CT'
}

def Graph(proteins: Set[str], options: Set[str], logger, start: int) -> None:
    data = {
    protein:{
    'sizes_outer': {},
    'sizes_inner': {},
    'distances_outer': {},
    'distances_inner': {},
    'touching_outer': {},
    'touching_inner': {}
    } for protein in proteins}
    for protein in proteins:
        for fileType in data[protein]:
            try:
                with open(protein+'/'+fileType+'_'+protein+".csv") as file:
                    for i, row in enumerate(csv.reader(file)):
                        if i == 0:
                            continue
                        if 'inner' in fileType:
                            number = float(''.join(row[-2]))
                            key = tuple(sorted(row[0:-2]))
                        else:
                            number = float(''.join(row[-1]))
                            key = tuple(sorted(row[0:-1]))
                        if not key in data[protein][fileType]:
                            data[protein][fileType][key] = []
                        if 'sizes' in fileType and "use_radius" in options:
                            number = math.sqrt(number/math.pi)
                        data[protein][fileType][key].append(number)
            except:
                raise FileNotFoundError('Could not read: '+fileType+'_'+protein+'.csv')
            for key in data[protein][fileType].keys():
                data[protein][fileType][key] = np.array(data[protein][fileType][key])

    for protein in proteins:
        proteinData = data[protein]
        for fileType, fileData in proteinData.items():
            if not fileType in options:
                continue
            if not fileData:
                continue
            if "remove_outliers" in options:
                for key in fileData.keys():
                    fileData[key] = fileData[key][abs(fileData[key] - np.mean(fileData[key])) <= 2 * np.std(fileData[key])]
            edge = (min((value.min() if value.size != 0 else float("inf")) for value in fileData.values()), max((value.max() if value.size != 0 else 0) for value in fileData.values()))
            total = sum((sum(value) if value.size != 0 else 0) for value in fileData.values())
            bin_num = math.ceil(math.log2(sum(len(value) for value in fileData.values())))+10
            for value in fileData.values():
                for i in value:
                    i = total/i
            
            width = math.ceil(math.sqrt(len(fileData)))
            height = math.ceil(len(fileData)/width)
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            fig2, axes = plt.subplots(height, width, sharex=False, sharey=False)
            axes = axes.flatten() if width*height>1 else (axes,)

            bin_width = (edge[1] - edge[0]) / bin_num

            for i, (key, values) in enumerate(sorted(fileData.items())):
                hist, bin_edges = np.histogram(values, bin_num//len(fileData), range=edge)
                title = key if type(key) == str else ', '.join(key)
                bins = np.delete(bin_edges, len(bin_edges) - 1)
                offset = bin_width * i
                if key in color_dict:
                    color = color_dict[key]
                    bar = ax1.bar(bins + offset, hist, width=bin_width, align='center',
                                  label=title, color=color)
                else:
                    bar = ax1.bar(bins + offset, hist, width=bin_width, align='center',
                                  label=title)
                    color = bar[0].get_facecolor()
                axes[i].hist(values, bins=bin_num, range=edge, color=color)
                axes[i].set_title(title+' '+string.capwords(fileType.replace('_', ' '))+' for '+name_dict.get(protein, protein), fontsize='small')
                axes[i].axvline(np.mean(values), color='red')
                axes[i].set_xlabel('$\AA^2$' if 'sizes' in fileType and 'use_radius' not in options else '$\AA$')
                axes[i].set_ylabel('Frequency')
                if 'log_scale' in options:
                    axes[i].set_xscale('log')

            ax1.set_title(string.capwords(fileType.replace('_', ' '))+' for '+name_dict.get(protein, protein), fontsize='large')
            ax1.legend()
            ax1.set_xlabel('$\AA^2$' if 'sizes' in fileType and 'use_radius' not in options else '$\AA$')
            ax1.set_ylabel('Frequency')
            
            plt.tight_layout()

            if 'save' in options:
                try:
                    mkdir('../images')
                except:
                    pass
                try:
                    fig1.savefig('../images/'+fileType+'_'+protein+'_all.png')
                    if height * width > 1:
                        fig2.savefig('../images/'+fileType+'_'+protein+'_seperate.png')
                except:
                    raise OSError("Could not access images folder.")
            if 'stats' in options:
                with open('../statistics/'+fileType+'_'+protein+'_stats.csv', 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(['']+[', '.join(key) for key in fileData])
                    writer.writerow(['Mean']+[np.mean(fileData[key]) for key in fileData])
                    writer.writerow(['Standard Deviation']+[np.std(fileData[key]) for key in fileData])
        if 'show' in options:
            plt.show()
        plt.close(fig2)
        plt.close(fig1)
    end=time()
    logger.debug('Finished graphing.', extra={'offset': timedelta(seconds=end-start), 'code': 'Main'})