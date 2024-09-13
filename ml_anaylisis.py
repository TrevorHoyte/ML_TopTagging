import functools
import pandas as pd
import numpy as np
import tensorflow as tf
import os
import pathlib
import statistics


def read_data(train_file,test_file,ttbar_file,qcd_file,Trimdata):
  
  df = pd.read_csv(train_file)
  df_test = pd.read_csv(test_file)
  df_ttbar = pd.read_csv(ttbar_file)
  df_qcd = pd.read_csv(qcd_file)
  
  #trim data keep only 100 const
  #list_columns=[]
  if (Trimdata):
    for i in range(301,601):
      label=str(i)
      df.pop(label)
      df_test.pop(label)
      df_ttbar.pop(label)
      df_qcd.pop(label)

  print(df.head())
  print(df.dtypes)

  print("-------------------")

  print(df_test.head())
  print(df_test.dtypes)

  target = df.pop('qcd=0, ttbar=1')
  target_test = df_test.pop('qcd=0, ttbar=1')
  target_tbar = df_ttbar.pop('qcd=0, ttbar=1')
  target_qcd = df_qcd.pop('qcd=0, ttbar=1')

  dataset = tf.data.Dataset.from_tensor_slices((df.values, target.values))
  dataset_test = tf.data.Dataset.from_tensor_slices((df_test.values, target_test.values))
  dataset_ttbar= tf.data.Dataset.from_tensor_slices((df_ttbar.values, target_tbar.values))
  dataset_qcd= tf.data.Dataset.from_tensor_slices((df_qcd.values, target_qcd.values))

  for feat, targ in dataset.take(5):
    print ('Features: {}, Target: {}'.format(feat, targ))



  train_dataset = dataset.shuffle(len(df)).batch(5000)
  test_dataset = dataset_test.shuffle(len(df_test)).batch(100)
  ttbar_dataset = dataset_test.shuffle(len(df_ttbar)).batch(100)
  qcd_dataset = dataset_test.shuffle(len(df_qcd)).batch(100)

  return train_dataset, test_dataset, ttbar_dataset,qcd_dataset

def get_compiled_model():
  model = tf.keras.Sequential([
    tf.keras.layers.Dense(300, activation='relu'),
    tf.keras.layers.Dense(102, activation='relu'),
    tf.keras.layers.Dense(12, activation='relu'),
    tf.keras.layers.Dense(6, activation='relu'),
    tf.keras.layers.Dense(1,activation='sigmoid')

    #tf.keras.layers.Dense(10, activation='relu'),
    #tf.keras.layers.Dense(10, activation='relu'),
    #tf.keras.layers.Dense(1)
  ])

  model.compile(optimizer='adam',
                loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
                metrics=['accuracy'])
  return model




######## plesase change home to location where al the files are############################################################################
home='/home/thoyte/dataset/ml/Results/'
############################################################################################################################################
file_raw=['vect3.csv','vect3_test.csv', 'ttbar_vect3.csv','qcd_vect3.csv']
file_random=['set_train_random.csv','set_test_random.csv','set_ttbar_random.csv','set_qcd_random.csv']
file_set1=['set_train_v3_full.csv','set_test_v3_full.csv','set_ttbar_v3_full.csv','set_qcd_v3_full.csv']
file_set2=['set_train_v3_full_divide.csv','set_test_v3_full_divide.csv','set_ttbar_v3_full_divide.csv','set_qcd_v3_full_divide.csv']
file_set3=['set_train_v3_trans.csv','set_test_v3_trans.csv','set_ttbar_v3_trans.csv','set_qcd_v3_trans.csv']
file_set4=['set_train_v3_trans_divide.csv','set_test_v3_trans_divide.csv','set_ttbar_v3_trans_divide.csv','set_qcd_v3_trans_divide.csv']

allfiles=[file_raw,file_random,file_set1,file_set2,file_set3,file_set4]

f=open('resultsfile.txt','w')


for items in allfiles:
  train_file=pathlib.Path(home+items[0])
  test_file=pathlib.Path(home+items[1])
  ttbar_file=pathlib.Path(home+items[2])
  qcd_file=pathlib.Path(home+items[3])

  data=read_data(train_file,test_file,ttbar_file,qcd_file,Trimdata=True)
  
  train_dataset=data[0]
  test_dataset=data[1]
  ttbar_dataset=data[2]
  qcd_dataset=data[3]


  results=[]
#####################100 runs for each file switch to smaller number if its too long###########
  for i in range(100):
    model = get_compiled_model()
    model.fit(train_dataset, epochs=40)

    test_loss, test_accuracy = model.evaluate(test_dataset)
    ttbar_loss, ttbar_accuracy=model.evaluate(ttbar_dataset)
    qcd_loss, qcd_accuracy=model.evaluate(qcd_dataset)

    print('\n\nTest Loss {}, Test Accuracy {}'.format(test_loss, test_accuracy))
    results.append((test_loss,test_accuracy,ttbar_accuracy,qcd_accuracy))

  #print to file results
  f.write("---------------------------------------------------------------------------- \n")
  f.write("The Train file"+items[0])
  std_loss=[]
  std_acc=[]
  std_tt=[]
  std_qcd=[]
  for i in range(len(results)):
    run=str(results[i])
    f.write(run)
    std_loss.append(results[i][0])
    std_acc.append(results[i][1])
    std_tt.append(results[i][2])
    std_qcd.append(results[i][3])
  
  f.write('Accuracy from testing: {0}  stdev:  {1} '.format(statistics.mean(std_acc),statistics.stdev(std_acc)))
  f.write('Accuracy from ttbar: {0}  stdev:  {1} '.format(statistics.mean(std_tt),statistics.stdev(std_tt)))
  f.write('Accuracy from qcd: {0}  stdev:  {1} '.format(statistics.mean(std_qcd),statistics.stdev(std_qcd)))
  f.write('loss from testing: {0}  stdev:  {1} '.format(statistics.mean(std_loss),statistics.stdev(std_loss)))


  
f.close()
  








