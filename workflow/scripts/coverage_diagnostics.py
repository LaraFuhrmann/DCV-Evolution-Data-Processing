import pandas as pd
import os
import numpy as np

def main(file_list, out_dir):

    df = pd.DataFrame(columns=['experiment','cell_line','passage','fname', 'ref', 'start', 'end', 'coverage'])

    for file in file_list:

        experiment = file.split("/")[-5]
        patient = file.split("/")[-4]
        date = file.split("/")[-3]

        if os.path.isfile(file)==True:
            df_temp = pd.read_csv(file, sep='\t', header=0 )

            column_coverage = df_temp.columns[-1]
            df_temp.rename(columns={column_coverage: "coverage"})

            df_temp['experiment']=experiment
            df_temp['cell_line']=patient
            df_temp['passage']=date

            df=pd.concat([df,df_temp])
        else:
            print('Not found: ',coverage_path)

    df.to_csv(out_dir+'coverage_all.csv')


    # plotting for each cell_line and passage separatly.
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    passages = ['passage_1', 'passage_5', 'passage_10']
    replicates = ['replicate_a', 'replicate_b', 'replicate_c', 'replicate_d', 'replicate_e']

    fig, ax = plt.subplots(5,3, figsize=(20, 15),constrained_layout=True, sharey=True,  sharex=True)
    fig.suptitle('Coverage along the genome', fontsize=20)


    for idx_rep, rep in enumerate(replicates):
        for idx_vp, vp in enumerate(passages):
            df_temp = df[df['cell_line'] == rep]
            df_temp = df_temp[df_temp['passage']==vp]

            ax[idx_rep][idx_vp].set_xticks(np.arange(0, len(x_pos)+1, 1000))
            ax[idx_rep][idx_vp].xaxis.set_tick_params(rotation=70)
            ax[idx_rep][idx_vp].plot(df_temp['pos'].values, df_temp['coverage'].values)

            ax[idx_rep][idx_vp].set_title(vp)
            ax[idx_rep][idx_vp].set_ylabel(rep)

            ax[idx_rep][idx_vp].plot(1000*np.ones(len(x_pos)), color='red')
            #ax[idx_rep][idx_vp].set_yscale('log')

    plt.savefig(out_dir+experiment+'_coverage_distribution_1000coverage.pdf',format='pdf')


if __name__ == "__main__":
    main(
        snakemake.input.fnames,
        snakemake.output.out_dir,
    )
