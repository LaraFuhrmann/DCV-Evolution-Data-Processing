import pandas as pd
import os
import numpy as np

def main(file_list, out_dir):

    df = pd.DataFrame(columns=['experiment','cell_line','passage','fname', 'ref', 'start', 'end', 'coverage'])

    for file in file_list:
        coverage_path = file.split('snvs.vcf')[0] + 'REGION_1/coverage.txt'

        if os.path.isfile(coverage_path)==True:
            df_temp = pd.read_csv(coverage_path, sep='\t', names=['fname', 'ref', 'start', 'end', 'coverage'])
            df_temp['experiment']=3
            df_temp['cell_line']=file.split("/")[-5]
            df_temp['passage']=file.split("/")[-4]
            df=pd.concat([df,df_temp])
        else:
            print('Not found: ',coverage_path)

    df.to_csv(out_dir+'coverage_all.csv')

    # create dataframe where each position is a column
    pos_min = int(np.min(np.array(df['start'])))
    pos_max = int(np.max(np.array(df['end'])))
    col_pos = ['cell_line', 'passage', 'start', 'end'] + [str(pos) for pos in range(pos_min,pos_max)]
    df_pos = pd.DataFrame(columns=col_pos)

    for iter_row, row in df.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        cov = int(row['coverage'])
        tmp_dict = {'cell_line':row['cell_line'],
                    'passage': row['passage'],
                    'start': start,
                    'end': end}
        for pos_tmp in range(start, end):
            tmp_dict.update({str(pos_tmp): cov})

        df_pos = df_pos.append(tmp_dict, ignore_index=True)

    df_pos.to_csv(out_dir+'coverage_position_wise.csv')

    # plotting for each cell_line and passage separatly.
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    for MT in df_pos['cell_line'].unique():

        for passage in sorted(df_pos['passage'].unique()):
            df = df_pos[df_pos['cell_line'] == MT]
            df = df_pos[df_pos['passage'] == passage]

            df = df.sort_values('start')
            if df.shape[0]>0:
                f = plt.figure()
                plt.plot(np.array(df['start']), np.array(df['coverage']))
                plt.hlines(1000,0,10000, color='red')
                f.savefig(out_dir+"/"+MT+"_"+passage+'_coverage.pdf')
                f.clf()
                plt.clf()
                plt.close()

        passages = sorted(df_pos['passage'].unique())

        fig, ax = plt.subplots(3,11, figsize=(20, 5),constrained_layout=True, sharey=True,  sharex=True)
        #fig.suptitle('Coverage of MT2_2 cell line', fontsize=20)

        df_cov_compact_grouped_null = df_pos.groupby(['cell_line', 'passage'], axis=0, dropna=False).mean()
        x_pos = sorted(df_cov_compact_grouped_null['pos'].unique())

        k=0
        for idx_vp, vp in enumerate(passages):
            df_cell_line = df_cov_compact_grouped_null.loc[cell_line, vp]


            if (idx_vp%11 ==0) and (idx_vp>0):
                k+=1

            colum_idx = idx_vp%11
            ax[k][colum_idx].set_xticks(np.arange(0, len(x_pos)+1, 1000))
            ax[k][colum_idx].xaxis.set_tick_params(rotation=70)
            ax[k][colum_idx].plot(x_pos, df_cell_line)
            ax[k][colum_idx].set_title('VP'+str((idx_vp+1)*10))
            ax[k][colum_idx].plot(1000*np.ones(len(x_pos)), color='red')
            ax[k][colum_idx].set_yscale('log')

        plt.savefig( out_dir+"/"+ MT + '_coverage.pdf',format='pdf')


if __name__ == "__main__":
    main(
        snakemake.input.fnames,
        snakemake.output.out_dir,
    )
