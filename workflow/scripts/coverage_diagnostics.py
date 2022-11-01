import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def main(file_list, fname_coverage_all, fname_coverage_plt):

    df = pd.DataFrame(columns=['experiment','cell_line','passage','fname', 'ref', 'coverage'])

    for file in file_list:

        experiment = file.split("/")[-6]
        patient = file.split("/")[-4]
        date = file.split("/")[-3]

        if os.path.isfile(file)==True:
            df_temp = pd.read_csv(file, sep='\t', header=0 )

            column_coverage = df_temp.columns[-1]
            df_temp = df_temp.rename(columns={column_coverage: "coverage"})

            df_temp['experiment']=experiment
            df_temp['cell_line']=patient
            df_temp['passage']=date
            df_temp['fname']=file
            cell_line = patient

            df=pd.concat([df,df_temp])
        else:
            print('Not found: ',coverage_path)

    df.to_csv(fname_coverage_all)

    if cell_line == "parental_stock":
        plt.plot(df['pos'].values, df['coverage'].values)
        plt.plot(1000*np.ones(len(df_temp['pos'].values)), color='red')
        plt.set_xticks(np.arange(0, len(df['pos'].values)+1, 1000))
        plt.xaxis.set_tick_params(rotation=70)
        plt.savefig(fname_coverage_plt,format='pdf')

        plt.plot(df['pos'].values, df['coverage'].values)
        plt.plot(1000*np.ones(len(df_temp['pos'].values)), color='red')
        plt.set_xticks(np.arange(0, len(df['pos'].values)+1, 1000))
        plt.xaxis.set_tick_params(rotation=70)
        plt.set_yscale('log')
        plt.savefig(fname_coverage_plt+'.log.pdf',format='pdf')


    else:
        # plotting for each cell_line and passage separatly.

        passages = ['passage_1', 'passage_5', 'passage_10']
        replicates = ['replicate_a', 'replicate_b', 'replicate_c', 'replicate_d', 'replicate_e']

        fig, ax = plt.subplots(5,3, figsize=(20, 15),constrained_layout=True, sharey=True,  sharex=True)
        fig.suptitle('Coverage along the genome', fontsize=20)


        for idx_rep, rep in enumerate(replicates):
            for idx_vp, vp in enumerate(passages):
                df_temp = df[df['cell_line'] == rep]
                df_temp = df_temp[df_temp['passage']==vp]

                ax[idx_rep][idx_vp].set_xticks(np.arange(0, len(df_temp['pos'].values)+1, 1000))
                ax[idx_rep][idx_vp].xaxis.set_tick_params(rotation=70)
                ax[idx_rep][idx_vp].plot(df_temp['pos'].values, df_temp['coverage'].values)

                ax[idx_rep][idx_vp].set_title(vp)
                ax[idx_rep][idx_vp].set_ylabel(rep)

                ax[idx_rep][idx_vp].plot(1000*np.ones(len(df_temp['pos'].values)), color='red')
                #ax[idx_rep][idx_vp].set_yscale('log')

        plt.savefig(fname_coverage_plt,format='pdf')

        for idx_rep, rep in enumerate(replicates):
            for idx_vp, vp in enumerate(passages):
                df_temp = df[df['cell_line'] == rep]
                df_temp = df_temp[df_temp['passage']==vp]

                ax[idx_rep][idx_vp].set_xticks(np.arange(0, len(df_temp['pos'].values)+1, 1000))
                ax[idx_rep][idx_vp].xaxis.set_tick_params(rotation=70)
                ax[idx_rep][idx_vp].plot(df_temp['pos'].values, df_temp['coverage'].values)

                ax[idx_rep][idx_vp].set_title(vp)
                ax[idx_rep][idx_vp].set_ylabel(rep)

                ax[idx_rep][idx_vp].plot(1000*np.ones(len(df_temp['pos'].values)), color='red')
                ax[idx_rep][idx_vp].set_yscale('log')

        plt.savefig(fname_coverage_plt+'.log.pdf',format='pdf')


if __name__ == "__main__":
    main(
        snakemake.input.fnames,
        snakemake.output.fname_coverage_all,
        snakemake.output.fname_coverage_plt,
    )
