import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def append_scores(input_pickle,input_csv,run):
    df_in = pd.read_csv(input_csv)
    df_out = pd.read_pickle(input_pickle)
    #df_out['signal_score'] = -9999
    signal_score_array = np.ones(len(df_out))*-99999.99


    for i in range(len(df_in.index)):

        temp_run = df_in['run_number'][i]
        temp_subrun = df_in['subrun_number'][i]
        temp_event = df_in['event_number'][i]
        signal_score = df_in['signal_score'][i]
        index_array = df_in.query('run_number == {:2f} & subrun_number == {:2f} & event_number == {:2f} '.format(temp_run,temp_subrun,temp_event)).index.values
        #print(index_array, signal_score)
        #df_out['signal_score'][index_array[0]] = signal_score
        signal_score_array[index_array[0]] = signal_score
        #print( signal_score_array[index_array[0]])

    df_out['signal_score'] = signal_score_array
    df_out.to_pickle("{}_CV_high_stats_with_scores.pkl".format(run))


run = "run3"
base_dir = "/home/lmlepin/Desktop/dm_sets/dark_tridents_analysis/{}_samples/".format(run)
input_pkl = base_dir + "{}_CV_high_stats_CNN.pkl".format(run)
input_csv = base_dir + "{}_CV_high_stats_CNN_scores_8441_steps.csv".format(run)

print(input_pkl)
print(input_csv)

append_scores(input_pkl,input_csv,run)