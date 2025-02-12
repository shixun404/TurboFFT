import torch as th
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from random import randint
from utils import load_data, load_data_single
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
if __name__ == '__main__':
    logN = [i for i in range(1, 26, 1)]
    logBS = [i for i in range(0, 29, 1)]

    surface_z = [ 2.5, 1.5, 1, 0.3, 0.2, 0.3]
    zticks = [ [0, 1, 2, 2.5], [0, 0.5, 1, 1.5], [0, 0.5, 1], [0, 0.1, 0.2, 0.3], [0, 0.1, 0.2], [0, 0.1, 0.2, 0.3]]
    zlim_label = [ ['0', '1', '2', '9.7'], ['0', '0.5', '1', '1.55'], ['0', '0.5', '8.1'], ['0', '0.1', '0.2', '0.3'],
                ['0', '0.1', '0.2'], ['0', '0.1', '0.2', '0.3']]
                
    # Creating figure
    # fig = plt.figure()
    # fsz = 32  
    # mksz = 40
    # plt.rcParams['lines.linewidth'] = 0.5
    # plt.rcParams["font.family"] = "Times New Roman"
    # plt.rc('font', size=fsz, weight='bold')
    # plt.rcParams["font.family"] = "Times New Roman"
    # plt.rcParams["hatch.color"] = 'white'
    # plt.rcParams['hatch.linewidth'] = 2.0
    # fig = plt.figure(figsize=(20, 10))
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rc('font', size=30, weight='bold')
    plt.rcParams["hatch.color"] = 'white'
    plt.rcParams['hatch.linewidth'] = 2.0
    fig = plt.figure(figsize=(20, 10))
    name_lists = [['cuFFT_data/Benchmark=1_cuFFT_FP64.csv','TurboFFT_data/Benchmark=1_TurboFFT_FP64.csv'], 
                     ]
    title_list = [['Compute Performance\n     A100 FP64   ', 'Memory Bandwidth\n   A100 FP64  '], 
               ]
    file_name = ['figure11.pdf']
    # c_peak = [19.5, ]
    j = 0
    # for name_list in name_lists:
    #     fig = plt.figure(figsize=(20, 10))
    #     ax1 = fig.add_subplot(1, 2, 1, projection='3d')

    #     # Second 3D subplot
    #     ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    #     for filename in name_list:
    #         data_list, data_tensor = load_data_single('../artifact_data/'+filename)
            


    #         # Creating plot
    #         for i in range(len(data_list)):
    #             data_point = data_list[i]
    #             color = 'r' if data_point[0] == 'turboFFT' else 'g'
    #             marker= 'o' if data_point[0] == 'turboFFT' else '^'
    #             label = None if i  != 0 else data_point[0]
    #             label = 'TurboFFT' if label == 'turboFFT' else label
    #             sct = ax1.scatter(data_point[1], data_point[2], data_point[4] / 1000, c=color, marker=marker, s=mksz, label=label)
    #             sct.set_sizes([mksz * 8])
    #             sct = ax2.scatter(data_point[1], data_point[2], data_point[5] / 1000, c=color, marker=marker, s=mksz)
    #             sct.set_sizes([mksz * 8])
    #             ax1.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, data_point[4] / 1000], color='grey', alpha=0.7)
    #             ax2.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, data_point[5] / 1000], color='grey', alpha=0.7)
    #                 # ax1.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, 2.5], color='grey', alpha=0.7)
    #                 # ax2.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, 1.5], color='grey', alpha=0.7)

    #         xx, yy = np.meshgrid(range(25), range(28))
    #         z = 0 * xx + 0 * yy + surface_z[j * 2]
    #         ax1.plot_surface(xx, yy, z, alpha=0.2, color='grey')
    #         z = 0 * xx + 0 * yy + surface_z[j * 2 + 1]
    #         ax2.plot_surface(xx, yy, z, alpha=0.2, color='grey')
    #         # First 3D subplot with elevation = 30, azimuth = 45
    #         ax1.view_init(elev=20, azim=-135)

    #         # Second 3D subplot with elevation = 60, azimuth = 120
    #         ax2.view_init(elev=20, azim=-135)
    #         # Adding labels and title
    #         ax1.set_xlabel('\nlogN',linespacing=2)
    #         ax1.set_ylabel('\nlogBS',linespacing=2)
    #         # z = 1
    #         ax1.zaxis.set_rotate_label(False)
    #         ax1.set_zlabel('TFLOPS',rotation=90)
    #         # ax1.zaxis.set_label_coords(1, -0.05)
    #         # ax1.text(x=0, y=1, z=1, s="z", color='red', size=8, transform = ax1.transAxes)
    #         ax1.set_zticks(zticks[j * 2])
    #         ax1.set_zticklabels(zlim_label[j * 2],fontsize=fsz - 5, rotation=90)
    #         ax2.set_xlabel('\nlogN',linespacing=2)
    #         ax2.set_ylabel('\nlogBS',linespacing=2)
    #         ax2.set_zlabel('TB/s',rotation=90)
    #         ax2.set_zticks(zticks[j * 2 + 1])
    #         ax2.set_zticklabels(zlim_label[j * 2 + 1],fontsize=fsz - 5,rotation=90)

    #         loc_x = 10
    #         loc_y = 25
    #         # ax1.text(loc_x - 2, loc_y, zticks[j * 2][-1] * 1.2, title_list[j][0] )
    #         # ax2.text(loc_x, loc_y, zticks[j * 2 + 1][-1] * 1.2, title_list[j][1] )
            
    #         ax1.text(loc_x, loc_y, zticks[j * 2][-1] * 0.9, '   Roofline Model' )
            
    #         ax2.text(loc_x, loc_y, zticks[j * 2 + 1][-1] * 0.9, '   Roofline Model' )
            

    #     # t2.set_position([0.5, 2])

    #     # fig.tight_layout()
    #     # fig.subplots_adjust(top=0.95)
    #     top_v = 0.7
    #     plt.subplots_adjust(bottom=0, wspace=0, hspace=-0.3, top=top_v)
    #     handles, labels = ax1.get_legend_handles_labels()
        
    #     lgnd = fig.legend(handles, labels, bbox_to_anchor=(0.67, top_v * 0.85), ncol=2,  framealpha=0)
        

    #     #change the marker size manually for both lines
        
        
    #     # Show plot
    #     # plt.show()
    #     # fig.tight_layout()
        
    #     # fig.subplots_adjust(top=top_v)
    colors = ['g', 'r', 'blue', 'orange']
    markers = ['^', 'o', 's', 'P']
    labels = ['cuFFT', 'turboFFT w/o FT', 'turboFFT w/ FT', 'turboFFT Err. Injection']
    for name_list in name_lists:
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')

        # Second 3D subplot
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        id = 0
        for filename in name_list:
            
            data_list, data_tensor = load_data_single('../artifact_data/'+filename)
            


            # Creating plot
            for i in range(len(data_list)):
                data_point = data_list[i]
                color = colors[id]
                marker= markers[id]
                label = None if i  != 0 else labels[id]
                ax1.scatter(data_point[1], data_point[2], data_point[4] / 1000, c=color, marker=marker, label=label)
                
                ax2.scatter(data_point[1], data_point[2], data_point[5] / 1000, c=color, marker=marker)
               
                ax1.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, data_point[4] / 1000], color='grey', alpha=0.7)
                ax2.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, data_point[5] / 1000], color='grey', alpha=0.7)
                    # ax1.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, 2.5], color='grey', alpha=0.7)
                    # ax2.plot([data_point[1], data_point[1]], [data_point[2], data_point[2]], zs=[0, 1.5], color='grey', alpha=0.7)
            id += 1
        xx, yy = np.meshgrid(range(25), range(28))
        z = 0 * xx + 0 * yy + surface_z[j * 2]
        ax1.plot_surface(xx, yy, z, alpha=0.2, color='grey')
        z = 0 * xx + 0 * yy + surface_z[j * 2 + 1]
        ax2.plot_surface(xx, yy, z, alpha=0.2, color='grey')
        # First 3D subplot with elevation = 30, azimuth = 45
        ax1.view_init(elev=20, azim=-135)

        # Second 3D subplot with elevation = 60, azimuth = 120
        ax2.view_init(elev=20, azim=-135)
        # Adding labels and title
        ax1.set_xlabel('\nlogN',linespacing=2)
        ax1.set_ylabel('\nlogBS',linespacing=2)
        ax1.set_zlabel('TFLOPS')
        ax1.set_zticks(zticks[j * 2])
        ax1.set_zticklabels(zlim_label[j * 2])
        ax2.set_xlabel('\nlogN',linespacing=2)
        ax2.set_ylabel('\nlogBS',linespacing=2)
        ax2.set_zlabel('TB/s')
        ax2.set_zticks(zticks[j * 2 + 1])
        ax2.set_zticklabels(zlim_label[j * 2 + 1])

        loc_x = 10
        loc_y = 25
        ax1.text(loc_x - 2, loc_y, zticks[j * 2][-1] * 1.2, title_list[j][0] )
        ax2.text(loc_x, loc_y, zticks[j * 2 + 1][-1] * 1.2, title_list[j][1] )
        
        ax1.text(loc_x, loc_y, zticks[j * 2][-1] * 0.9, '   Roofline Model' )
        
        ax2.text(loc_x, loc_y, zticks[j * 2 + 1][-1] * 0.9, '   Roofline Model' )

        # t2.set_position([0.5, 2])

        # fig.tight_layout()
        # fig.subplots_adjust(top=0.95)
        plt.subplots_adjust(bottom=0, wspace=0, hspace=-0.3, top=0.9)
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, bbox_to_anchor=(0.67, 0.1), ncol=2)

        fig.savefig(f'../artifact_figures/{file_name[j]}',bbox_inches='tight')
        print(f'./artifact_figures/{file_name[j]}')
        j += 1